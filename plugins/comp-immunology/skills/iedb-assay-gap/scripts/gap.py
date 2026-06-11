# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""Assay-gap analysis -- epitopes with a T-cell assay but NO MHC-binding assay for an allele.

For an allele X, this finds epitopes that:
  * HAVE >= 1 T-cell assay restricted by X   (tcell_search)
  * have ZERO MHC-binding assays against X   (mhc_search)

i.e. epitopes shown to be T-cell antigenic in the context of X whose direct binding to X has never
been experimentally measured -- prime candidates for binding validation and a quantifiable gap in
the binding data.

Allele matching is grounded in the MRO (MHC Restriction Ontology) the IQ-API exposes through the
``mhc_allele_iri_search`` array field (each record carries its allele IRI plus every ancestor IRI):

  * exact + higher-resolution sub-alleles  -> mhc_allele_iri_search=cs.{<IRI>}
  * serotype rollup (default on)           -> exact name match on the serotype ancestor
                                              (e.g. HLA-A2, HLA-DR4), derived from the ontology.

Every matched record is tagged with the resolution at which it matched (``exact``, ``sub-allele``,
or ``serotype``) so the rollup is fully auditable.

The HTTP / pagination / ontology layer is imported from :mod:`iedb_toolkit.core`; this module keeps
only the gap's domain logic. ``resolve_allele`` is a thin adapter over ``core.resolve_entity``.
"""

import argparse
import csv
import html
import json
import os
import re
import sys
from datetime import datetime, timezone

import core
from core import fetch_all

# Fields pulled for each T-cell assay row (enough to characterise a candidate).
TCELL_SELECT = ",".join([
    "tcell_id", "structure_id", "linear_sequence", "linear_sequence_length",
    "qualitative_measure", "mhc_allele_name", "mhc_allele_resolution",
    "mhc_class", "source_organism_name", "parent_source_antigen_name",
    "disease_names", "assay_names", "host_organism_name", "pubmed_id",
])
# We only need the structure_id (epitope key) from the binding side; the match
# resolution is known from which query produced the row, so no allele columns.
MHC_SELECT = "structure_id"


# --------------------------------------------------------------------------- #
# Allele / ontology resolution
# --------------------------------------------------------------------------- #
def resolve_allele(allele, timeout):
    """Look up an allele's IRI, MRO ancestor IRIs, class, and derived serotype.

    Thin adapter over :func:`iedb_toolkit.core.resolve_entity` (which already returns iri /
    mhc_class / ancestors / serotype_name / serotype_iri). Ancestors are narrowed to the MRO
    ontology, matching the standalone script's output shape.
    """
    info = core.resolve_entity("allele", allele, timeout=timeout)
    return {
        "allele": allele,
        "iri": info["iri"],
        "mhc_class": info.get("mhc_class"),
        "ancestors": [a for a in info["ancestors"]
                      if isinstance(a.get("iri"), str) and a["iri"].startswith("MRO")],
        "serotype_name": info.get("serotype_name"),
        "serotype_iri": info.get("serotype_iri"),
    }


def build_match_queries(info, mode):
    """Return list of (resolution_label, filter_params) for this allele + mode.

    The exact (iri-contains) and serotype (name-eq) record sets are disjoint:
    a serotype-level record (e.g. HLA-A2) is an *ancestor* of the allele, so its
    iri_search array does NOT contain the allele's IRI.
    """
    queries = [("exact", {"mhc_allele_iri_search": f"cs.{{{info['iri']}}}"})]
    if mode == "serotype" and info["serotype_name"]:
        queries.append(
            ("serotype", {"mhc_allele_name": f"eq.{info['serotype_name']}"}))
    return queries


# --------------------------------------------------------------------------- #
# Core gap computation
# --------------------------------------------------------------------------- #
def collect_side(endpoint, select, info, mode, order_key, timeout):
    """Fetch all rows for an allele's match set, tagging each row's resolution.

    Returns (rows, counts) where each row gets a `_resolution` key:
    'exact', 'sub-allele', or 'serotype'. Paging is delegated to core.fetch_all, which
    keyset-paginates by `order_key` (the endpoint's unique primary key).
    """
    rows = []
    counts = {}
    for label, filt in build_match_queries(info, mode):
        fetched = fetch_all(endpoint, filt, select=select, order_key=order_key,
                            timeout=timeout, label=label)
        counts[label] = len(fetched)
        for r in fetched:
            if label == "serotype":
                r["_resolution"] = "serotype"
            else:  # exact (iri-contains) query; name column present only T-cell side
                nm = r.get("mhc_allele_name")
                r["_resolution"] = "sub-allele" if (nm and nm != info["allele"]) else "exact"
        rows.extend(fetched)
    return rows, counts


def _as_items(value, sep="|"):
    """Normalise an IQ-API field to a list of non-empty strings.

    Some columns are JSON arrays (disease_names), others are sep-joined strings
    (assay_names is '|'-separated). Handle both.
    """
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v).strip() for v in value if v not in (None, "")]
    return [p.strip() for p in str(value).split(sep) if p.strip()]


def aggregate_tcell(rows):
    """Group T-cell assay rows by epitope structure_id."""
    epis = {}
    for r in rows:
        sid = r["structure_id"]
        e = epis.get(sid)
        if e is None:
            e = epis[sid] = {
                "structure_id": sid,
                "linear_sequence": r.get("linear_sequence"),
                "length": r.get("linear_sequence_length"),
                "mhc_class": r.get("mhc_class"),
                "source_organism": r.get("source_organism_name"),
                "parent_antigen": r.get("parent_source_antigen_name"),
                "host_organism": r.get("host_organism_name"),
                "diseases": set(),
                "assays": set(),
                "allele_names": set(),
                "resolutions": set(),
                "n_assays": 0, "n_pos": 0, "n_neg": 0, "n_other": 0,
                "pubmeds": set(),
            }
        e["n_assays"] += 1
        e["resolutions"].add(r["_resolution"])
        if r.get("mhc_allele_name"):
            e["allele_names"].add(r["mhc_allele_name"])
        qm = (r.get("qualitative_measure") or "").lower()
        if qm.startswith("positive"):
            e["n_pos"] += 1
        elif qm.startswith("negative"):
            e["n_neg"] += 1
        else:
            e["n_other"] += 1
        for d in _as_items(r.get("disease_names")):
            e["diseases"].add(d)
        for a in _as_items(r.get("assay_names")):
            e["assays"].add(a)
        if r.get("pubmed_id"):
            e["pubmeds"].add(r["pubmed_id"])
        # keep first non-null peptide/organism if earlier row was null
        for k_row, k_epi in (("linear_sequence", "linear_sequence"),
                             ("source_organism_name", "source_organism"),
                             ("parent_source_antigen_name", "parent_antigen")):
            if not e[k_epi] and r.get(k_row):
                e[k_epi] = r.get(k_row)
    return epis


def binding_structs_by_resolution(rows):
    """Sets of structure_ids that have a binding assay, split by resolution."""
    exact = set()
    sero = set()
    for r in rows:
        sid = r["structure_id"]
        if r["_resolution"] == "serotype":
            sero.add(sid)
        else:
            exact.add(sid)
    return exact, sero


def compute_gap(allele, mode, timeout):
    """Run the full gap analysis for one allele."""
    info = resolve_allele(allele, timeout)
    sys.stderr.write(
        f"  allele={allele}  iri={info['iri']}  class={info['mhc_class']}  "
        f"serotype={info['serotype_name']}\n")

    tcell_rows, tcell_counts = collect_side(
        "tcell_search", TCELL_SELECT, info, mode, "tcell_id", timeout)
    mhc_rows, mhc_counts = collect_side(
        "mhc_search", MHC_SELECT, info, mode, "elution_id", timeout)

    epis = aggregate_tcell(tcell_rows)
    tcell_exact = {s for s, e in epis.items() if e["resolutions"] & {"exact", "sub-allele"}}
    tcell_sero = {s for s, e in epis.items() if "serotype" in e["resolutions"]}
    tcell_all = set(epis.keys())

    mhc_exact, mhc_sero = binding_structs_by_resolution(mhc_rows)
    mhc_all = mhc_exact | mhc_sero

    # Gap sets
    gap_exact = tcell_exact - mhc_exact                # exact resolution only
    gap_inclusive = tcell_all - mhc_all                # exact + serotype rollup

    # Build candidate records (use the inclusive gap as the headline set;
    # mark which ones are also in the exact-only gap).
    candidates = []
    for sid in gap_inclusive:
        e = epis[sid]
        res = sorted(e["resolutions"])
        candidates.append({
            "structure_id": sid,
            "linear_sequence": e["linear_sequence"],
            "length": e["length"],
            "mhc_class": e["mhc_class"],
            "n_tcell_assays": e["n_assays"],
            "n_tcell_positive": e["n_pos"],
            "n_tcell_negative": e["n_neg"],
            "n_tcell_other": e["n_other"],
            "tcell_resolution": "+".join(res),
            "in_exact_gap": sid in gap_exact,
            "source_organism": e["source_organism"],
            "parent_antigen": e["parent_antigen"],
            "host_organism": e["host_organism"],
            "diseases": sorted(e["diseases"]),
            "assays": sorted(e["assays"]),
            "example_pubmed": next(iter(sorted(e["pubmeds"])), None) if e["pubmeds"] else None,
            "iedb_url": f"https://www.iedb.org/epitope/{sid}",
        })
    # Sort: positive T-cell first, then total assays, then peptide
    candidates.sort(key=lambda c: (-c["n_tcell_positive"], -c["n_tcell_assays"],
                                   c["linear_sequence"] or ""))

    return {
        "allele": allele,
        "match_mode": mode,
        "iri": info["iri"],
        "mhc_class": info["mhc_class"],
        "serotype_name": info["serotype_name"],
        "serotype_iri": info["serotype_iri"],
        "ancestors": info["ancestors"],
        "assay_counts": {"tcell": tcell_counts, "mhc": mhc_counts},
        "epitope_counts": {
            "tcell_epitopes_exact": len(tcell_exact),
            "tcell_epitopes_serotype": len(tcell_sero),
            "tcell_epitopes_total": len(tcell_all),
            "binding_epitopes_exact": len(mhc_exact),
            "binding_epitopes_serotype": len(mhc_sero),
            "binding_epitopes_total": len(mhc_all),
            "gap_exact": len(gap_exact),
            "gap_inclusive": len(gap_inclusive),
        },
        "n_candidates": len(candidates),
        "n_candidates_positive": sum(1 for c in candidates if c["n_tcell_positive"] > 0),
        "candidates": candidates,
    }


# --------------------------------------------------------------------------- #
# Output rendering
# --------------------------------------------------------------------------- #
def write_csv(result, path):
    fields = ["structure_id", "linear_sequence", "length", "mhc_class",
              "n_tcell_assays", "n_tcell_positive", "n_tcell_negative",
              "tcell_resolution", "in_exact_gap", "source_organism",
              "parent_antigen", "host_organism", "diseases", "assays",
              "example_pubmed", "iedb_url"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for c in result["candidates"]:
            w.writerow([
                c["structure_id"], c["linear_sequence"], c["length"], c["mhc_class"],
                c["n_tcell_assays"], c["n_tcell_positive"], c["n_tcell_negative"],
                c["tcell_resolution"], c["in_exact_gap"], c["source_organism"],
                c["parent_antigen"], c["host_organism"],
                "; ".join(c["diseases"]), "; ".join(c["assays"]),
                c["example_pubmed"], c["iedb_url"],
            ])


def _esc(x):
    return html.escape("" if x is None else str(x))


def render_html(results, mode, generated_at, table_cap=500):
    css = """
:root{--bg:#0f1419;--card:#1a2230;--ink:#e6edf3;--mut:#8b98a9;--acc:#4fd1c5;
--pos:#3fb950;--warn:#d29922;--line:#2d3748;--chip:#243044}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1200px;margin:0 auto;padding:32px 20px 80px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:20px;margin:38px 0 12px;
border-bottom:1px solid var(--line);padding-bottom:8px}
h3{font-size:16px;margin:22px 0 8px;color:var(--acc)}
.sub{color:var(--mut);font-size:13px;margin:0 0 24px}
a{color:var(--acc);text-decoration:none}a:hover{text-decoration:underline}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:14px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .n{font-size:26px;font-weight:700}.card .l{color:var(--mut);font-size:12px;margin-top:2px}
.card.hl{border-color:var(--acc)}.card.hl .n{color:var(--acc)}
.meta{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:12px 0;font-size:13.5px}
.meta b{color:var(--ink)}.meta code{background:var(--chip);padding:1px 6px;border-radius:4px;color:var(--acc)}
table{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}
th,td{text-align:left;padding:7px 9px;border-bottom:1px solid var(--line);vertical-align:top}
th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none;white-space:nowrap}
th:hover{color:var(--acc)}tr:hover td{background:#161d29}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.chip{display:inline-block;background:var(--chip);border-radius:4px;padding:1px 6px;font-size:11px;margin:1px}
.pos{color:var(--pos);font-weight:600}.warn{color:var(--warn)}
.tag{font-size:10.5px;text-transform:uppercase;letter-spacing:.04em;padding:1px 6px;border-radius:4px}
.tag.ex{background:#15331f;color:var(--pos)}.tag.se{background:#3a2f12;color:var(--warn)}
.note{background:#11202b;border-left:3px solid var(--acc);padding:12px 16px;border-radius:6px;margin:14px 0;font-size:13.5px;color:#bcd}
.tbl-wrap{max-height:640px;overflow:auto;border:1px solid var(--line);border-radius:10px}
.foot{color:var(--mut);font-size:12px;margin-top:40px;border-top:1px solid var(--line);padding-top:14px}
"""
    js = """
function sortTable(t,n,numeric){const tb=t.tBodies[0];const rows=[...tb.rows];
const dir=t.getAttribute('data-sd')==='asc'?-1:1;t.setAttribute('data-sd',dir===1?'asc':'desc');
rows.sort((a,b)=>{let x=a.cells[n].getAttribute('data-v')??a.cells[n].innerText;
let y=b.cells[n].getAttribute('data-v')??b.cells[n].innerText;
if(numeric){x=parseFloat(x)||0;y=parseFloat(y)||0;return (x-y)*dir;}
return x.localeCompare(y)*dir;});rows.forEach(r=>tb.appendChild(r));}
"""
    parts = [f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>IEDB T-cell-without-binding gap report</title><style>{css}</style></head>
<body><div class="wrap">
<h1>Epitopes with T-cell assays but no MHC-binding assay</h1>
<p class="sub">IEDB Query API (IQ-API) &middot; match mode: <b>{_esc(mode)}</b> &middot;
generated {_esc(generated_at)}</p>
<div class="note"><b>What this is.</b> For each MHC allele below, these are epitopes that have at
least one <b>T-cell assay restricted by that allele</b> but <b>zero MHC-binding assays</b>
measured against the same allele &mdash; epitopes proven T-cell antigenic in the allele's context
whose direct binding has never been tested. A &ldquo;gap&rdquo; means <b>no recorded binding
assay</b>, not proven non-binding.</div>
"""]

    # cross-allele overview
    parts.append('<h2>Overview</h2><div class="tbl-wrap"><table data-sd="desc">'
                 '<thead><tr>'
                 '<th onclick="sortTable(this.closest(\'table\'),0,false)">Allele</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),1,false)">Class</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),2,false)">Serotype</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),3,true)">T-cell epitopes</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),4,true)">Binding-tested</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),5,true)">Gap (incl. serotype)</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),6,true)">Gap (exact only)</th>'
                 '<th onclick="sortTable(this.closest(\'table\'),7,true)">Gap w/ +ve T-cell</th>'
                 '</tr></thead><tbody>')
    for r in results:
        ec = r["epitope_counts"]
        parts.append(
            f'<tr><td class="mono">{_esc(r["allele"])}</td>'
            f'<td>{_esc(r["mhc_class"])}</td>'
            f'<td class="mono">{_esc(r["serotype_name"] or "&mdash;")}</td>'
            f'<td data-v="{ec["tcell_epitopes_total"]}">{ec["tcell_epitopes_total"]:,}</td>'
            f'<td data-v="{ec["binding_epitopes_total"]}">{ec["binding_epitopes_total"]:,}</td>'
            f'<td data-v="{ec["gap_inclusive"]}"><b class="warn">{ec["gap_inclusive"]:,}</b></td>'
            f'<td data-v="{ec["gap_exact"]}">{ec["gap_exact"]:,}</td>'
            f'<td data-v="{r["n_candidates_positive"]}"><b class="pos">{r["n_candidates_positive"]:,}</b></td>'
            f'</tr>')
    parts.append('</tbody></table></div>')

    # per-allele sections
    for r in results:
        ec = r["epitope_counts"]
        ac = r["assay_counts"]
        anc = ", ".join(f'{a["name"]}' for a in r["ancestors"] if a["name"]) or "&mdash;"
        parts.append(f'<h2>{_esc(r["allele"])}</h2>')
        parts.append(
            '<div class="cards">'
            f'<div class="card"><div class="n">{ec["tcell_epitopes_total"]:,}</div>'
            f'<div class="l">epitopes w/ T-cell assay</div></div>'
            f'<div class="card"><div class="n">{ec["binding_epitopes_total"]:,}</div>'
            f'<div class="l">epitopes binding-tested</div></div>'
            f'<div class="card hl"><div class="n">{ec["gap_inclusive"]:,}</div>'
            f'<div class="l">GAP &mdash; T-cell but no binding</div></div>'
            f'<div class="card"><div class="n">{r["n_candidates_positive"]:,}</div>'
            f'<div class="l">of which have a <b>positive</b> T-cell assay</div></div>'
            '</div>')
        parts.append(
            '<div class="meta">'
            f'<b>Allele IRI</b> <code>{_esc(r["iri"])}</code> &nbsp; '
            f'<b>Class</b> {_esc(r["mhc_class"])} &nbsp; '
            f'<b>Serotype ancestor</b> <code>{_esc(r["serotype_name"] or "none")}</code> '
            f'(<code>{_esc(r["serotype_iri"] or "&mdash;")}</code>)<br>'
            f'<b>Match resolution split</b> &mdash; T-cell epitopes: '
            f'{ec["tcell_epitopes_exact"]:,} exact/sub-allele, '
            f'{ec["tcell_epitopes_serotype"]:,} serotype-level; '
            f'binding-tested: {ec["binding_epitopes_exact"]:,} exact/sub-allele, '
            f'{ec["binding_epitopes_serotype"]:,} serotype-level.<br>'
            f'<b>Assay rows scanned</b> &mdash; T-cell: '
            f'{sum(ac["tcell"].values()):,} &nbsp; MHC-binding: {sum(ac["mhc"].values()):,}.<br>'
            f'<b>Ontology ancestors</b> &mdash; {_esc(anc)}.'
            '</div>')

        cands = r["candidates"]
        shown = cands[:table_cap]
        cap_note = (f' &middot; showing top {table_cap:,} of {len(cands):,} '
                    f'(full set in CSV/JSON)') if len(cands) > table_cap else ''
        parts.append(f'<h3>Candidate epitopes ({len(cands):,}){cap_note}</h3>')
        parts.append('<div class="tbl-wrap"><table data-sd="desc"><thead><tr>'
                     '<th onclick="sortTable(this.closest(\'table\'),0,false)">Peptide</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),1,true)">Len</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),2,true)">T-cell assays</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),3,true)">+ve</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),4,true)">&minus;ve</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),5,false)">Restriction res.</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),6,false)">Source organism</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),7,false)">Parent antigen</th>'
                     '<th onclick="sortTable(this.closest(\'table\'),8,false)">Diseases</th>'
                     '<th>IEDB</th></tr></thead><tbody>')
        for c in shown:
            res_html = "".join(
                f'<span class="tag {"ex" if t!="serotype" else "se"}">{_esc(t)}</span> '
                for t in c["tcell_resolution"].split("+"))
            dz = ", ".join(c["diseases"][:3]) + (" &hellip;" if len(c["diseases"]) > 3 else "")
            pos_cls = ' class="pos"' if c["n_tcell_positive"] > 0 else ''
            parts.append(
                f'<tr>'
                f'<td class="mono">{_esc(c["linear_sequence"])}</td>'
                f'<td data-v="{c["length"] or 0}">{_esc(c["length"])}</td>'
                f'<td data-v="{c["n_tcell_assays"]}">{c["n_tcell_assays"]}</td>'
                f'<td data-v="{c["n_tcell_positive"]}"{pos_cls}>{c["n_tcell_positive"]}</td>'
                f'<td data-v="{c["n_tcell_negative"]}">{c["n_tcell_negative"]}</td>'
                f'<td>{res_html}</td>'
                f'<td>{_esc(c["source_organism"])}</td>'
                f'<td>{_esc(c["parent_antigen"])}</td>'
                f'<td>{_esc(dz)}</td>'
                f'<td><a href="{_esc(c["iedb_url"])}" target="_blank">#{c["structure_id"]}</a></td>'
                f'</tr>')
        parts.append('</tbody></table></div>')

    # methodology + caveats
    parts.append("""
<h2>Methodology &amp; caveats</h2>
<div class="meta">
<b>Data source.</b> IEDB Query API (IQ-API), endpoints <code>tcell_search</code>,
<code>mhc_search</code>; PostgREST. Public, no auth.<br><br>
<b>Definition.</b> For allele X, an epitope is a candidate if it has &ge;1 T-cell assay
whose MHC restriction matches X and 0 MHC-binding assays whose tested allele matches X
(&ldquo;same-allele&rdquo; gap).<br><br>
<b>Allele matching.</b> Grounded in the MHC Restriction Ontology (MRO) exposed via the
<code>mhc_allele_iri_search</code> array. <span class="tag ex">exact</span>/<span class="tag ex">sub-allele</span>
matches use <code>mhc_allele_iri_search=cs.{IRI}</code> (the allele and any higher-resolution typing).
<span class="tag se">serotype</span> rollup additionally matches records annotated at the serotype
ancestor (e.g. HLA-A2, HLA-DR4) derived from the ontology and shown per allele above. Each candidate
is tagged with the resolution(s) at which its T-cell evidence exists. Set <code>--match-mode exact</code>
to disable the serotype rollup.<br><br>
<b>Caveats.</b> (1) Serotype detection is best-effort; the chosen serotype is printed for audit.
(2) &ldquo;Has a T-cell assay&rdquo; counts positive and negative assays (broken out); epitopes with a
positive T-cell assay are the highest-interest candidates. (3) Non-peptidic structures are out of scope.
(4) Restriction-annotation resolution varies across records. (5) A gap reflects absence of a
<i>recorded</i> binding assay in IEDB, not experimental evidence of non-binding.
</div>
""")
    parts.append(f'<div class="foot">Generated by the <code>iedb-assay-gap</code> tool '
                 f'&middot; {_esc(generated_at)} &middot; data &copy; IEDB</div>')
    parts.append(f'<script>{js}</script></div></body></html>')
    return "".join(parts)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="iedb-gap",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--allele", action="append", required=True,
                    help="MHC allele IEDB name, e.g. 'HLA-A*02:01'. Repeatable.")
    ap.add_argument("--match-mode", choices=["serotype", "exact"], default="serotype",
                    help="serotype: include serotype-level rollup (default). "
                         "exact: only the allele + higher-resolution sub-alleles.")
    ap.add_argument("--out", "-o", default=".", help="Output directory (default: cwd).")
    ap.add_argument("--timeout", type=int, default=90, help="Per-request timeout (s).")
    ap.add_argument("--table-cap", type=int, default=500,
                    help="Max candidate rows rendered per allele in HTML "
                         "(full set always written to CSV/JSON).")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    results = []
    for allele in args.allele:
        sys.stderr.write(f"\n=== {allele} (mode={args.match_mode}) ===\n")
        res = compute_gap(allele, args.match_mode, args.timeout)
        results.append(res)
        safe = re.sub(r"[^A-Za-z0-9]+", "_", allele).strip("_")
        with open(os.path.join(args.out, f"gap_{safe}.json"), "w") as fh:
            json.dump(res, fh, indent=2, default=str)
        write_csv(res, os.path.join(args.out, f"gap_{safe}.csv"))
        ec = res["epitope_counts"]
        sys.stderr.write(
            f"  -> gap (incl. serotype): {ec['gap_inclusive']:,} epitopes "
            f"({res['n_candidates_positive']:,} with a positive T-cell assay); "
            f"gap (exact only): {ec['gap_exact']:,}\n")

    html_out = render_html(results, args.match_mode, generated_at, args.table_cap)
    report_path = os.path.join(args.out, "report.html")
    with open(report_path, "w") as fh:
        fh.write(html_out)
    sys.stderr.write(f"\nReport written: {report_path}\n")
    return 0
