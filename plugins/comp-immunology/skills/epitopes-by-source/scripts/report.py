# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""Epitopes-by-source report.

Pull every IEDB epitope from a given SOURCE ORGANISM and build a self-contained HTML report
(+ CSV + JSON) via the IQ-API. Given a source organism (a pathogen / taxon -- e.g. "EBV",
"SARS-CoV-2", "NCBITaxon:10376"), this fetches all epitopes whose source organism is that taxon
**or any descendant taxon**, using the ontology subtree trick the IQ-API exposes through the
``source_organism_iri_search`` array column::

    source_organism_iri_search=cs.{NCBITaxon:10376}
        -> Human herpesvirus 4 (Epstein Barr virus) AND every EBV strain.

The unit of analysis is an IEDB epitope (``epitope_search.structure_id``, one row per epitope
aggregated across all assays). The report summarises the set (linear vs discontinuous; assay
coverage; MHC class; top antigens / hosts / alleles; length histogram; per-strain counts) and
a sortable table.

Source resolution is forgiving: an NCBITaxon id (``NCBITaxon:10376`` or bare ``10376``), a built-in
abbreviation (EBV, CMV, HIV-1, ...), or a substring of the IEDB organism name.

The HTTP / pagination / counting layer is imported from :mod:`iedb_toolkit.core`; this module keeps
only the report's domain logic.
"""

import argparse
import csv
import html
import json
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone

from core import count, iter_all, fetch_all, fetch_one

EPITOPE_ENDPOINT = "epitope_search"
NCBITAXON_RE = re.compile(r"^(?:NCBITaxon:)?(\d+)$", re.IGNORECASE)

# Endpoints that carry the *scalar* source/host organism columns used for name<->IRI
# resolution. epitope_search has only the list-valued *_names column, so resolution
# (substring scan / IRI->name) must run against the assay endpoints.
RESOLVE_ENDPOINTS = ("tcell_search", "mhc_search")
SOURCE_COLS = ("source_organism_iri", "source_organism_name", "source_organism_iri_search")
HOST_COLS = ("host_organism_iri", "host_organism_name", "host_organism_iri_search")

# Columns pulled for each epitope (everything the report / CSV / breakdowns need).
EPITOPE_SELECT = ",".join([
    "structure_id", "linear_sequence", "linear_sequence_length", "structure_type",
    "parent_source_antigen_names", "parent_source_antigen_iris",
    "source_organism_names", "source_organism_iris",
    "host_organism_names", "host_organism_iris",
    "mhc_allele_names", "mhc_classes", "qualitative_measures", "assay_names",
    "reference_ids", "pubmed_ids", "journal_names",
    "tcell_ids", "bcell_ids", "elution_ids",
])

# Common abbreviation -> species-level NCBITaxon rollup (so the subtree trick captures
# all strains). The resolver re-derives the canonical IEDB name from the API, so a
# wrong IRI here surfaces immediately as a name mismatch. EBV (10376) is verified.
ALIASES = {
    "ebv":               {"iri": "NCBITaxon:10376",   "name": "Human herpesvirus 4 (Epstein Barr virus)"},
    "epstein-barr virus": {"iri": "NCBITaxon:10376",  "name": "Human herpesvirus 4 (Epstein Barr virus)"},
    "epstein barr virus": {"iri": "NCBITaxon:10376",  "name": "Human herpesvirus 4 (Epstein Barr virus)"},
    "cmv":               {"iri": "NCBITaxon:10359",   "name": "Human herpesvirus 5 (Cytomegalovirus)"},
    "hcmv":              {"iri": "NCBITaxon:10359",   "name": "Human herpesvirus 5 (Cytomegalovirus)"},
    "hiv":               {"iri": "NCBITaxon:11676",   "name": "Human immunodeficiency virus 1"},
    "hiv-1":             {"iri": "NCBITaxon:11676",   "name": "Human immunodeficiency virus 1"},
    "hiv1":              {"iri": "NCBITaxon:11676",   "name": "Human immunodeficiency virus 1"},
    "hbv":               {"iri": "NCBITaxon:10407",   "name": "Hepatitis B virus"},
    "hcv":               {"iri": "NCBITaxon:11103",   "name": "Hepatitis C virus"},
    "sars-cov-2":        {"iri": "NCBITaxon:2697049", "name": "Severe acute respiratory syndrome coronavirus 2"},
    "sars2":             {"iri": "NCBITaxon:2697049", "name": "Severe acute respiratory syndrome coronavirus 2"},
    "covid":             {"iri": "NCBITaxon:2697049", "name": "Severe acute respiratory syndrome coronavirus 2"},
    "influenza a":       {"iri": "NCBITaxon:11320",   "name": "Influenza A virus"},
    "influenza":         {"iri": "NCBITaxon:11320",   "name": "Influenza A virus"},
    "flu":               {"iri": "NCBITaxon:11320",   "name": "Influenza A virus"},
    "dengue":            {"iri": "NCBITaxon:12637",   "name": "Dengue virus"},
    "mtb":               {"iri": "NCBITaxon:1773",    "name": "Mycobacterium tuberculosis"},
    "m. tuberculosis":   {"iri": "NCBITaxon:1773",    "name": "Mycobacterium tuberculosis"},
    "tuberculosis":      {"iri": "NCBITaxon:1773",    "name": "Mycobacterium tuberculosis"},
    "plasmodium falciparum": {"iri": "NCBITaxon:5833", "name": "Plasmodium falciparum"},
}

HOST_ALIASES = {
    "human":        {"iri": "NCBITaxon:9606",  "name": "Homo sapiens (human)"},
    "homo sapiens": {"iri": "NCBITaxon:9606",  "name": "Homo sapiens (human)"},
    "mouse":        {"iri": "NCBITaxon:10090", "name": "Mus musculus (mouse)"},
    "mus musculus": {"iri": "NCBITaxon:10090", "name": "Mus musculus (mouse)"},
}


# --------------------------------------------------------------------------- #
# Taxon (source / host organism) resolution
# --------------------------------------------------------------------------- #
def _iri_to_name(cols, iri, timeout):
    """NCBITaxon IRI -> canonical organism name via the scalar columns on the assay endpoints."""
    iri_col, name_col, _arr = cols
    for ep in RESOLVE_ENDPOINTS:
        rec = fetch_one(ep, {iri_col: f"eq.{iri}"}, select=f"{iri_col},{name_col}", timeout=timeout)
        if rec and rec.get(name_col):
            return rec[name_col]
    return None


def _substring_scan(cols, term, timeout, max_rows=4000):
    """Bounded case-insensitive substring scan: returns {iri: name} for matching organisms.

    Pages the assay endpoints via core.iter_all, which uses each endpoint's own primary-key
    cursor (tcell_id / elution_id) -- a unique, sortable keyset column.
    """
    iri_col, name_col, _arr = cols
    cand = {}
    for ep in RESOLVE_ENDPOINTS:
        for row in iter_all(ep, {name_col: f"ilike.*{term}*"},
                            select=f"{iri_col},{name_col}", timeout=timeout,
                            max_rows=max_rows, label="resolve"):
            if row.get(iri_col):
                cand.setdefault(row[iri_col], row.get(name_col))
        if cand:
            break
    return cand


def resolve_taxon(cols, name_or_iri, label, aliases, timeout=90):
    """Resolve an organism name/abbreviation/id to {query, iri, name, resolution}.

    Tiers: (1) NCBITaxon id/CURIE -> use directly; (2) built-in alias; (3) bounded
    case-insensitive substring scan of the organism name (prefer exact / 'X (common)'
    match; fail fast with the candidate list if ambiguous).
    """
    raw = name_or_iri.strip()

    m = NCBITAXON_RE.match(raw)
    if m:
        iri = f"NCBITaxon:{m.group(1)}"
        return {"query": raw, "iri": iri, "name": _iri_to_name(cols, iri, timeout) or iri,
                "resolution": "ncbitaxon-id"}

    key = raw.lower()
    if key in aliases:
        iri = aliases[key]["iri"]
        return {"query": raw, "iri": iri,
                "name": _iri_to_name(cols, iri, timeout) or aliases[key]["name"],
                "resolution": "alias"}

    cand = _substring_scan(cols, raw, timeout)
    if not cand:
        raise ValueError(
            f"No {label} organism in IEDB matching {raw!r}. Try an NCBITaxon id "
            f"(e.g. 'NCBITaxon:10376' or '10376'), a known alias "
            f"({', '.join(sorted(aliases))}), or part of the organism's IEDB name.")
    if len(cand) == 1:
        iri = next(iter(cand))
        return {"query": raw, "iri": iri, "name": cand[iri], "resolution": "substring"}

    q = raw.lower()
    preferred = [i for i, nm in cand.items()
                 if nm and (nm.strip().lower() == q or nm.strip().lower().startswith(q + " ("))]
    if len(preferred) == 1:
        iri = preferred[0]
        return {"query": raw, "iri": iri, "name": cand[iri], "resolution": "substring"}

    shown = "; ".join(f"{nm} [{i}]" for i, nm in list(cand.items())[:8])
    raise ValueError(
        f"{label.capitalize()} organism {raw!r} is ambiguous ({len(cand)} matches): {shown}"
        + (" ..." if len(cand) > 8 else "")
        + ". Re-run with the exact name or its NCBITaxon id.")


# --------------------------------------------------------------------------- #
# Per-epitope classification & breakdowns
# --------------------------------------------------------------------------- #
def _as_items(value):
    """Normalise a column to a list of non-empty strings (list -> elements; scalar -> [scalar])."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v).strip() for v in value if v not in (None, "")]
    s = str(value).strip()
    return [s] if s else []


def _list_len(value):
    return len(value) if isinstance(value, list) else (1 if value not in (None, "") else 0)


def _assay_types(value):
    """assay_names elements are '|'-joined (assay|method); return the distinct tokens."""
    out = set()
    for el in _as_items(value):
        for tok in str(el).split("|"):
            tok = tok.strip()
            if tok:
                out.add(tok)
    return sorted(out)


def _most_specific(names):
    """Pick the most specific source-organism name (prefer one naming a strain; else longest)."""
    names = _as_items(names)
    if not names:
        return None
    strained = [n for n in names if "strain" in n.lower() or "isolate" in n.lower()]
    return max(strained or names, key=len)


def classify_epitope(row):
    """Derive a flat, report-ready record from one epitope_search row."""
    seq = row.get("linear_sequence")
    stype = (row.get("structure_type") or "").strip()
    is_linear = stype.lower() == "linear peptide"

    quals = _as_items(row.get("qualitative_measures"))
    classes = set(_as_items(row.get("mhc_classes")))
    if "I" in classes and "II" in classes:
        class_bucket = "both"
    elif "I" in classes:
        class_bucket = "I"
    elif "II" in classes:
        class_bucket = "II"
    elif classes:
        class_bucket = "other"
    else:
        class_bucket = "none"

    pubmeds = _as_items(row.get("pubmed_ids"))
    n_refs = _list_len(row.get("reference_ids")) or len(pubmeds)

    return {
        "structure_id": row.get("structure_id"),
        "linear_sequence": seq,
        "length": row.get("linear_sequence_length"),
        "structure_type": stype,
        "linear_bucket": "linear" if is_linear else "discontinuous",
        "antigens": _as_items(row.get("parent_source_antigen_names")),
        "source_organisms": _as_items(row.get("source_organism_names")),
        "strain": _most_specific(row.get("source_organism_names")),
        "hosts": _as_items(row.get("host_organism_names")),
        "alleles": _as_items(row.get("mhc_allele_names")),
        "mhc_classes": sorted(classes),
        "class_bucket": class_bucket,
        "qualitative_measures": quals,
        "assay_types": _assay_types(row.get("assay_names")),
        "any_positive": any(q.lower().startswith("positive") for q in quals),
        "any_negative": any(q.lower().startswith("negative") for q in quals),
        "has_tcell": _list_len(row.get("tcell_ids")) > 0,
        "has_bcell": _list_len(row.get("bcell_ids")) > 0,
        "has_mhc_ligand": _list_len(row.get("elution_ids")) > 0,
        "n_refs": n_refs,
        "pubmed_ids": pubmeds,
        "example_pubmed": pubmeds[0] if pubmeds else None,
        "journals": _as_items(row.get("journal_names")),
        "iedb_url": f"https://www.iedb.org/epitope/{row.get('structure_id')}",
    }


def apply_post_filters(epitopes, linear_only=False, positive_only=False, with_mhc=False,
                       mhc_class=None, allele=None):
    """Apply the client-side (list-union / structure-type) filters main exposes as flags.

    Server-side filters (source/host subtree, length range) are pushed into the query; these
    operate on the aggregated per-epitope columns that can't be expressed as a single PostgREST
    predicate. The allele filter keeps an exact allele and its higher-resolution sub-typings
    (e.g. HLA-A*02:01:01) but NOT the broader serotype (HLA-A2).
    """
    if linear_only:
        epitopes = [e for e in epitopes if e["linear_bucket"] == "linear"]
    if positive_only:
        epitopes = [e for e in epitopes if e["any_positive"]]
    if with_mhc:
        epitopes = [e for e in epitopes if e["class_bucket"] != "none"]
    if mhc_class:
        epitopes = [e for e in epitopes if mhc_class in e["mhc_classes"]]
    if allele:
        target = allele.strip().upper()
        prefix = target + ":"
        epitopes = [e for e in epitopes
                    if any(a.upper() == target or a.upper().startswith(prefix) for a in e["alleles"])]
    return epitopes


def _top(counter, n=None):
    return [{"label": k, "count": v} for k, v in counter.most_common(n)]


CLASS_ORDER = ["I", "II", "both", "other", "none"]
CLASS_LABEL = {"I": "MHC class I", "II": "MHC class II", "both": "class I + II",
               "other": "other", "none": "no MHC data"}


def compute_breakdowns(epitopes, top_n=15):
    """Aggregate counters / histograms / top-N tables over the classified epitopes."""
    total = len(epitopes)
    linear = sum(1 for e in epitopes if e["linear_bucket"] == "linear")

    antigen_c, host_c, allele_c, strain_c = Counter(), Counter(), Counter(), Counter()
    length_c = Counter()
    class_c = Counter(e["class_bucket"] for e in epitopes)
    for e in epitopes:
        for a in set(e["antigens"]):
            antigen_c[a] += 1
        for h in set(e["hosts"]):
            host_c[h] += 1
        for al in set(e["alleles"]):
            allele_c[al] += 1
        if e["strain"]:
            strain_c[e["strain"]] += 1
        if e["linear_bucket"] == "linear" and e["length"]:
            length_c[e["length"]] += 1

    return {
        "total": total,
        "linear": linear,
        "discontinuous": total - linear,
        "with_tcell": sum(1 for e in epitopes if e["has_tcell"]),
        "with_bcell": sum(1 for e in epitopes if e["has_bcell"]),
        "with_mhc_ligand": sum(1 for e in epitopes if e["has_mhc_ligand"]),
        "with_positive": sum(1 for e in epitopes if e["any_positive"]),
        "mhc_class_dist": [{"label": CLASS_LABEL[k], "count": class_c[k]}
                           for k in CLASS_ORDER if class_c.get(k)],
        "top_antigens": _top(antigen_c, top_n),
        "top_hosts": _top(host_c, top_n),
        "top_alleles": _top(allele_c, top_n),
        "strains": _top(strain_c),
        "length_hist": [{"label": str(L), "count": c} for L, c in sorted(length_c.items())],
        "n_antigens": len(antigen_c),
        "n_hosts": len(host_c),
        "n_alleles": len(allele_c),
        "n_strains": len(strain_c),
    }


# --------------------------------------------------------------------------- #
# CSV / JSON output
# --------------------------------------------------------------------------- #
CSV_FIELDS = [
    "structure_id", "linear_sequence", "length", "structure_type",
    "source_organism", "parent_antigens", "mhc_alleles", "mhc_classes",
    "host_organisms", "qualitative_measures", "assay_types",
    "has_tcell", "has_bcell", "has_mhc_ligand", "any_positive",
    "n_references", "pubmed_ids", "journals", "iedb_url",
]


def write_csv(epitopes, path):
    def j(xs):
        return "; ".join(str(x) for x in xs)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_FIELDS)
        for e in epitopes:
            w.writerow([
                e["structure_id"], e["linear_sequence"] or "", e["length"] if e["length"] is not None else "",
                e["structure_type"], e["strain"] or "", j(e["antigens"]), j(e["alleles"]),
                j(e["mhc_classes"]), j(e["hosts"]), j(e["qualitative_measures"]), j(e["assay_types"]),
                e["has_tcell"], e["has_bcell"], e["has_mhc_ligand"], e["any_positive"],
                e["n_refs"], j(e["pubmed_ids"]), j(e["journals"]), e["iedb_url"],
            ])


def write_json(payload, path):
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)


# --------------------------------------------------------------------------- #
# HTML report (light theme, self-contained: inline CSS + vanilla-JS sortable table)
# --------------------------------------------------------------------------- #
def _esc(x):
    return html.escape("" if x is None else str(x))


CSS = """
:root{--bg:#fafafa;--card:#ffffff;--ink:#1a1a1a;--mut:#6a6a6a;--line:#d8d8d8;
--hover:#f0f4ff;--link:#1665c1;--bar:#2980b9;--bar-bg:#eef2f7;--pos:#197a19;--chip-bg:#f0f0f0}
*{box-sizing:border-box}
body{margin:0 auto;max-width:1340px;padding:32px 24px 80px;background:var(--bg);color:var(--ink);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
h1{font-size:26px;margin:0 0 4px}
h2{font-size:20px;margin:38px 0 12px;border-bottom:2px solid var(--line);padding-bottom:6px}
h3{font-size:14px;margin:14px 0 6px;color:var(--mut);font-weight:600;text-transform:uppercase;letter-spacing:.03em}
.sub{color:var(--mut);font-size:13px;margin:0 0 20px}
a{color:var(--link);text-decoration:none}a:hover{text-decoration:underline}
.tagq{display:inline-block;font-size:11px;padding:1px 8px;border-radius:10px;background:#e3edf7;color:var(--link);font-weight:600}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(148px,1fr));gap:12px;margin:14px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px}
.card .n{font-size:25px;font-weight:700}.card .l{color:var(--mut);font-size:12px;margin-top:2px}
.card.hl{border-color:var(--link)}.card.hl .n{color:var(--link)}
.meta{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0;font-size:13px;line-height:1.8}
.meta code{background:var(--chip-bg);padding:1px 6px;border-radius:4px;font-size:12px}
.note{background:#fff8e1;border:1px solid #f5e1a4;border-radius:8px;padding:12px 16px;margin:14px 0;font-size:13.5px}
.charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:6px 30px;margin:6px 0 8px}
.barrow{display:flex;align-items:center;gap:10px;margin:3px 0;font-size:13px}
.barlbl{flex:0 0 240px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.bartrack{flex:1;background:var(--bar-bg);border-radius:3px;height:15px;overflow:hidden}
.barfill{display:block;height:15px;border-radius:3px}
.barval{flex:0 0 58px;text-align:right;color:var(--mut);font-variant-numeric:tabular-nums}
.tbl-wrap{max-height:760px;overflow:auto;border:1px solid var(--line);border-radius:8px;margin-top:8px}
table{width:100%;border-collapse:collapse;font-size:13px;background:var(--card)}
th,td{text-align:left;padding:7px 9px;border-bottom:1px solid var(--line);vertical-align:top}
th{position:sticky;top:0;background:#efefef;cursor:pointer;user-select:none;white-space:nowrap;font-weight:600}
th:hover{background:#e5e5e5}tr:hover td{background:var(--hover)}
.mono{font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace}
.tag{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;font-weight:600;white-space:nowrap}
.tag-tcell{background:#e3edf7;color:#1f6fb2}.tag-bcell{background:#fbe6d4;color:#b25a16}
.tag-ligand{background:#ece2f5;color:#6d3f9c}
.tag-classI{background:#fff4e0;color:#b86a00}.tag-classII{background:#ece2f5;color:#6d3f9c}
.tag-other{background:#eee;color:#555}.tag-disc{background:#f0f0f0;color:#888}
.pos{color:var(--pos);font-weight:700}
.more{color:var(--mut);font-size:11px}
.foot{color:var(--mut);font-size:12px;margin-top:40px;border-top:1px solid var(--line);padding-top:14px}
"""

JS = """
function sortTable(t,n,numeric){const tb=t.tBodies[0];const rows=[...tb.rows];
const dir=t.getAttribute('data-sd')==='asc'?-1:1;t.setAttribute('data-sd',dir===1?'asc':'desc');
rows.sort((a,b)=>{let x=a.cells[n].getAttribute('data-v')??a.cells[n].innerText;
let y=b.cells[n].getAttribute('data-v')??b.cells[n].innerText;
if(numeric){x=parseFloat(x)||0;y=parseFloat(y)||0;return (x-y)*dir;}
return x.localeCompare(y)*dir;});rows.forEach(r=>tb.appendChild(r));}
"""

CLASS_BAR_COLORS = {"MHC class I": "#b86a00", "MHC class II": "#6d3f9c",
                    "class I + II": "#2980b9", "other": "#888", "no MHC data": "#bbb"}


def _bar(label, n, maxval, color):
    pct = (n / maxval * 100) if maxval else 0
    return (f'<div class="barrow"><span class="barlbl" title="{_esc(label)}">{_esc(label)}</span>'
            f'<span class="bartrack"><span class="barfill" style="width:{pct:.1f}%;background:{color}"></span></span>'
            f'<span class="barval">{n:,}</span></div>')


def _chart(title, items, color):
    """Render a labelled block of horizontal bars. `color` is a hex string or label->hex fn."""
    if not items:
        return ""
    maxval = max(it["count"] for it in items) or 1
    def col(lbl):
        return color(lbl) if callable(color) else color
    rows = "".join(_bar(it["label"], it["count"], maxval, col(it["label"])) for it in items)
    return f'<div><h3>{_esc(title)}</h3>{rows}</div>'


def _card(n, label, hl=False):
    return (f'<div class="card{" hl" if hl else ""}"><div class="n">{n:,}</div>'
            f'<div class="l">{label}</div></div>')


def _join_trunc(items, k=2):
    items = [i for i in items if i]
    if not items:
        return "&mdash;"
    head = "; ".join(_esc(i) for i in items[:k])
    if len(items) > k:
        head += f' <span class="more">(+{len(items) - k})</span>'
    return head


def _class_tags(classes):
    out = []
    for c in classes:
        if c == "I":
            out.append('<span class="tag tag-classI">I</span>')
        elif c == "II":
            out.append('<span class="tag tag-classII">II</span>')
        else:
            out.append(f'<span class="tag tag-other">{_esc(c)}</span>')
    return " ".join(out)


def _assay_tags(e):
    out = []
    if e["has_tcell"]:
        out.append('<span class="tag tag-tcell">T&#8209;cell</span>')
    if e["has_bcell"]:
        out.append('<span class="tag tag-bcell">B&#8209;cell</span>')
    if e["has_mhc_ligand"]:
        out.append('<span class="tag tag-ligand">MHC&nbsp;ligand</span>')
    return " ".join(out) or "&mdash;"


def render_html(payload, table_cap=1000):
    md = payload["metadata"]
    bd = payload["breakdowns"]
    eps = payload["epitopes"]
    src = md["source"]
    total = bd.get("total", 0)

    parts = [f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>IEDB epitopes &mdash; {_esc(src['name'])}</title><style>{CSS}</style></head>
<body>
<h1>Epitopes from {_esc(src['name'])}</h1>
<p class="sub">source <code>{_esc(src['iri'])}</code>
<span class="tagq">{_esc(src['resolution'])}</span> &middot;
IEDB Query API &middot; endpoint <code>epitope_search</code> &middot;
generated {_esc(md['generated_at'])} &middot; data &copy; IEDB</p>"""]

    # provenance
    qf = "; ".join(f"{k}={v}" for k, v in md["query_filter"].items())
    prov = [f'<b>Query.</b> All epitopes whose source organism is <code>{_esc(src["name"])}</code> '
            f'or any descendant taxon, via the NCBITaxon subtree filter '
            f'<code>{_esc(qf)}</code>.']
    if md.get("host_filter"):
        hf = md["host_filter"]
        prov.append(f'<b>Host filter.</b> Restricted to host <code>{_esc(hf["name"])}</code> '
                    f'(<code>{_esc(hf["iri"])}</code>).')
    pf = md["post_filters"]
    active_pf = [k.replace("_", " ") + (f"={v}" if v not in (True, None) else "")
                 for k, v in pf.items() if v not in (None, False)]
    if active_pf:
        prov.append(f'<b>Filters applied.</b> {_esc(", ".join(active_pf))} '
                    f'(reduced {md["n_fetched"]:,} &rarr; {md["n_after_filters"]:,} epitopes).')
    est = md.get("count_estimated")
    line = f'<b>Size.</b> estimated matches {est:,} &middot; fetched {md["n_fetched"]:,}'
    if md.get("truncated"):
        line += f' &middot; <b>capped</b> by --max-rows={md["max_rows"]:,} (not all matches fetched)'
    prov.append(line + ".")
    parts.append('<div class="meta">' + "<br>".join(prov) + '</div>')

    if total == 0:
        parts.append('<div class="note">No epitopes found for this source organism under the '
                     'current filters. A null result means no matching annotation in IEDB, not '
                     'biological absence.</div>')

    # summary cards
    parts.append('<div class="cards">')
    parts.append(_card(total, "epitopes (total)", hl=True))
    parts.append(_card(bd["linear"], "linear peptides"))
    parts.append(_card(bd["discontinuous"], "discontinuous / other"))
    parts.append(_card(bd["with_positive"], "with a positive assay"))
    parts.append(_card(bd["with_tcell"], "with T-cell assay data"))
    parts.append(_card(bd["with_bcell"], "with B-cell assay data"))
    parts.append(_card(bd["with_mhc_ligand"], "with MHC-ligand data"))
    parts.append(_card(bd["n_strains"], "distinct source taxa"))
    parts.append('</div>')

    # breakdown charts
    if total:
        parts.append('<h2>Breakdowns</h2><div class="charts">')
        parts.append(_chart(f"Top source antigens ({bd['n_antigens']:,} distinct)",
                            bd["top_antigens"], "var(--bar)"))
        parts.append(_chart("MHC class distribution", bd["mhc_class_dist"],
                            lambda lbl: CLASS_BAR_COLORS.get(lbl, "#888")))
        parts.append(_chart("Peptide-length histogram (linear)", bd["length_hist"], "#2980b9"))
        parts.append(_chart(f"Top host organisms ({bd['n_hosts']:,} distinct)",
                            bd["top_hosts"], "#27ae60"))
        parts.append(_chart(f"Top MHC alleles ({bd['n_alleles']:,} distinct)",
                            bd["top_alleles"], "#b86a00"))
        if bd["strains"] and bd["n_strains"] > 1:
            parts.append(_chart(f"Source taxa under the rollup ({bd['n_strains']:,})",
                                bd["strains"][:15], "#7c3aed"))
        parts.append('</div>')

    # main epitope table
    if total:
        eps_sorted = sorted(eps, key=lambda e: (
            0 if e["any_positive"] else 1, -(e["n_refs"] or 0), -(e["length"] or 0),
            e["linear_sequence"] or ""))
        shown = eps_sorted[:table_cap]
        cap_note = (f' &middot; showing top {table_cap:,} of {total:,} (full set in CSV/JSON)'
                    if total > table_cap else '')
        parts.append(f'<h2>Epitopes ({total:,}){cap_note}</h2>')
        if total > table_cap:
            parts.append(f'<div class="note">Table truncated to the first {table_cap:,} epitopes '
                         f'(ranked: positive assay, then #references, then length). '
                         f'All {total:,} epitopes are in the CSV and JSON.</div>')
        cols = [("Epitope", True), ("Sequence", False), ("Len", True), ("Source antigen", False),
                ("Source organism", False), ("MHC", False), ("Host", False), ("Assays", False),
                ("+ve", True), ("Refs", True), ("PubMed", False)]
        head = "".join(
            f'<th onclick="sortTable(this.closest(\'table\'),{i},{str(num).lower()})">{_esc(name)}</th>'
            for i, (name, num) in enumerate(cols))
        parts.append(f'<div class="tbl-wrap"><table data-sd="desc"><thead><tr>{head}</tr></thead><tbody>')
        for e in shown:
            sid = e["structure_id"]
            if e["linear_sequence"]:
                seq_cell = f'<td class="mono" data-v="{_esc(e["linear_sequence"])}">{_esc(e["linear_sequence"])}</td>'
            else:
                seq_cell = '<td data-v=""><span class="tag tag-disc">discontinuous</span></td>'
            length = e["length"]
            mhc_cell = (_class_tags(e["mhc_classes"]) + " " + _join_trunc(e["alleles"], 2)).strip() or "&mdash;"
            assay_code = ("T" if e["has_tcell"] else "") + ("B" if e["has_bcell"] else "") + ("L" if e["has_mhc_ligand"] else "")
            pos = ('<span class="pos">&#10003;</span>', 1) if e["any_positive"] else ("&mdash;", 0)
            pm = e["example_pubmed"]
            pm_cell = (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{_esc(pm)}/" target="_blank">{_esc(pm)}</a>'
                       if pm else "&mdash;")
            parts.append(
                f'<tr>'
                f'<td data-v="{sid}"><a href="{_esc(e["iedb_url"])}" target="_blank">{_esc(sid)}</a></td>'
                f'{seq_cell}'
                f'<td data-v="{length or 0}">{_esc(length) if length is not None else "&mdash;"}</td>'
                f'<td>{_join_trunc(e["antigens"], 1)}</td>'
                f'<td>{_esc(e["strain"]) if e["strain"] else "&mdash;"}</td>'
                f'<td data-v="{_esc("".join(e["mhc_classes"]))}">{mhc_cell}</td>'
                f'<td>{_join_trunc(e["hosts"], 2)}</td>'
                f'<td data-v="{assay_code}">{_assay_tags(e)}</td>'
                f'<td data-v="{pos[1]}">{pos[0]}</td>'
                f'<td data-v="{e["n_refs"]}">{e["n_refs"]:,}</td>'
                f'<td>{pm_cell}</td>'
                f'</tr>')
        parts.append('</tbody></table></div>')

    # methodology + caveats
    parts.append("""
<h2>Methodology &amp; caveats</h2>
<div class="meta">
<b>Data source.</b> IEDB Query API (IQ-API), endpoint <code>epitope_search</code> (one row per
epitope, aggregated across all of that epitope's assays); PostgREST, public, no auth.<br>
<b>Source matching.</b> Epitopes are selected by the NCBITaxon subtree
(<code>source_organism_iri_search=cs.{IRI}</code>), which matches the taxon <i>and every
descendant taxon</i> &mdash; so a species-level query (e.g. EBV) includes all of its strains.<br>
<b>Assay coverage.</b> &ldquo;With T-cell / B-cell / MHC-ligand data&rdquo; counts epitopes that
carry at least one assay of that kind (T-cell assays, B-cell assays, or MHC mass-spec/elution
ligand assays). &ldquo;With a positive assay&rdquo; means &ge;1 qualitative measure begins with
&ldquo;Positive&rdquo; (positive and negative results are both present in IEDB).<br>
<b>Caveats.</b> (1) A null result means no matching annotation in IEDB, not biological absence.
(2) Discontinuous / non-peptidic epitopes have no linear sequence and are excluded from the
length histogram. (3) Counts are per epitope; an epitope contributing the same antigen/host/allele
more than once is counted once for that breakdown. (4) Estimated match counts are PostgreSQL-planner
estimates; the fetched/after-filter counts in the cards are exact.
</div>
""")
    parts.append(f'<div class="foot">Generated by the <code>epitopes-by-source</code> tool '
                 f'&middot; {_esc(md["generated_at"])} &middot; data &copy; IEDB &middot; '
                 f'click a column header to sort</div>')
    parts.append(f'<script>{JS}</script></body></html>')
    return "".join(parts)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="iedb-source",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source",
                    help="Source organism: an abbreviation (EBV, SARS-CoV-2, ...), an NCBITaxon id "
                         "('NCBITaxon:10376' or '10376'), or a substring of the IEDB organism name.")
    ap.add_argument("--out", "-o", default=".", help="Output directory (created if missing).")
    ap.add_argument("--max-rows", type=int, default=None,
                    help="Cap epitopes fetched (logs when it truncates).")
    ap.add_argument("--timeout", type=int, default=90, help="Per-request timeout (s).")
    ap.add_argument("--top-n", type=int, default=15,
                    help="Top-N length for antigen/host/allele breakdowns (default 15).")
    ap.add_argument("--table-cap", type=int, default=1000,
                    help="Max epitope rows rendered in the HTML (full set always in CSV/JSON).")
    ap.add_argument("--host", default=None,
                    help="Restrict to a host organism (name or NCBITaxon id), e.g. 'human'.")
    ap.add_argument("--positive-only", action="store_true",
                    help="Keep only epitopes with >=1 positive assay outcome.")
    ap.add_argument("--mhc-class", choices=["I", "II"], default=None,
                    help="Keep only epitopes tested against the given MHC class.")
    ap.add_argument("--allele", default=None,
                    help="Keep only epitopes restricted by this MHC allele (IEDB name, e.g. "
                         "'HLA-A*02:01'). Matches the allele and its higher-resolution sub-typings "
                         "(e.g. HLA-A*02:01:01); does NOT include the broader serotype (HLA-A2).")
    ap.add_argument("--with-mhc", action="store_true",
                    help="Keep only epitopes with MHC data (drop the 'no MHC data' bucket, "
                         "i.e. antibody/B-cell epitopes with no MHC restriction).")
    ap.add_argument("--linear-only", action="store_true",
                    help="Keep only linear peptide epitopes (drop discontinuous).")
    ap.add_argument("--min-length", type=int, default=None, help="Minimum peptide length (linear).")
    ap.add_argument("--max-length", type=int, default=None, help="Maximum peptide length (linear).")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # 1) resolve source (and optional host)
    src = resolve_taxon(SOURCE_COLS, args.source, "source", ALIASES, timeout=args.timeout)
    sys.stderr.write(f"source: {src['name']}  [{src['iri']}]  (resolved via {src['resolution']})\n")
    host = None
    if args.host:
        host = resolve_taxon(HOST_COLS, args.host, "host", HOST_ALIASES, timeout=args.timeout)
        sys.stderr.write(f"host:   {host['name']}  [{host['iri']}]\n")

    # 2) build server-side filters (subtree on source, optional subtree on host, length range)
    filters = {SOURCE_COLS[2]: f"cs.{{{src['iri']}}}"}
    if host:
        filters[HOST_COLS[2]] = f"cs.{{{host['iri']}}}"
    lo, hi = args.min_length, args.max_length
    if lo is not None and hi is not None:
        filters["and"] = f"(linear_sequence_length.gte.{lo},linear_sequence_length.lte.{hi})"
    elif lo is not None:
        filters["linear_sequence_length"] = f"gte.{lo}"
    elif hi is not None:
        filters["linear_sequence_length"] = f"lte.{hi}"

    # 3) size preflight (advisory)
    try:
        est = count(EPITOPE_ENDPOINT, filters, mode="estimated", timeout=args.timeout)
        sys.stderr.write(f"estimated epitopes: {est:,}\n")
        if args.max_rows and est > args.max_rows:
            sys.stderr.write(f"NOTE: --max-rows={args.max_rows:,} will cap the pull "
                             f"(~{est - args.max_rows:,} matching epitopes NOT fetched)\n")
    except Exception as e:  # advisory only
        est = -1
        sys.stderr.write(f"(size preflight skipped: {e})\n")

    # 4) fetch + classify
    raw_rows = fetch_all(EPITOPE_ENDPOINT, filters, select=EPITOPE_SELECT,
                         timeout=args.timeout, max_rows=args.max_rows, label="epitopes")
    n_fetched = len(raw_rows)
    epitopes = [classify_epitope(r) for r in raw_rows]

    # 5) post-filters (list-union columns / structure type)
    epitopes = apply_post_filters(
        epitopes, linear_only=args.linear_only, positive_only=args.positive_only,
        with_mhc=args.with_mhc, mhc_class=args.mhc_class, allele=args.allele)

    # 6) breakdowns + payload
    breakdowns = compute_breakdowns(epitopes, top_n=args.top_n)
    payload = {
        "metadata": {
            "tool": "epitopes-by-source",
            "generated_at": generated_at,
            "endpoint": EPITOPE_ENDPOINT,
            "source": src,
            "host_filter": host,
            "query_filter": filters,
            "post_filters": {
                "positive_only": args.positive_only,
                "with_mhc": args.with_mhc,
                "mhc_class": args.mhc_class,
                "allele": args.allele,
                "linear_only": args.linear_only,
                "min_length": args.min_length,
                "max_length": args.max_length,
            },
            "count_estimated": est,
            "n_fetched": n_fetched,
            "n_after_filters": len(epitopes),
            "max_rows": args.max_rows,
            "truncated": bool(args.max_rows and n_fetched >= args.max_rows),
        },
        "breakdowns": breakdowns,
        "epitopes": epitopes,
    }

    # 7) write artifacts
    safe = re.sub(r"[^A-Za-z0-9]+", "_", src["iri"]).strip("_") or "source"
    csv_path = os.path.join(args.out, f"epitopes_{safe}.csv")
    json_path = os.path.join(args.out, f"epitopes_{safe}.json")
    html_path = os.path.join(args.out, "report.html")
    write_csv(epitopes, csv_path)
    write_json(payload, json_path)
    with open(html_path, "w") as fh:
        fh.write(render_html(payload, table_cap=args.table_cap))

    sys.stderr.write(
        f"\n{len(epitopes):,} epitopes "
        f"({breakdowns['linear']:,} linear, {breakdowns['with_positive']:,} with a positive assay) "
        f"for {src['name']}.\n"
        f"  {csv_path}\n  {json_path}\n  {html_path}\n")
    return 0
