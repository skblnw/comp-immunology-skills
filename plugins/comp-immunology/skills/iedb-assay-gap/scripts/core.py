# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""Shared IEDB Query API (IQ-API) client.

The IQ-API (https://query-api.iedb.org) is a public, no-auth PostgREST front end over the IEDB
Postgres tables. Each table is an endpoint; queries use PostgREST filter syntax; results are JSON
(default) or CSV (``Accept: text/csv``).

This module packages the efficient-access patterns so any ad-hoc IEDB pull is fast and correct:

  * KEYSET (cursor) pagination -- O(1) per page, unlike OFFSET which rescans skipped rows and
    degrades at high offsets. A single response is capped at 10,000 rows, so paging is required.
  * count=estimated  -- the cheap PostgreSQL-planner count for gauging result size before a pull
    (count=exact is accurate but slow over hundreds of thousands of rows).
  * The ONTOLOGY SUBTREE trick -- every search endpoint carries ``<entity>_iri_search`` array
    columns holding an annotation's IRI plus ALL its ancestor IRIs in the relevant ontology.
    ``<field>=cs.{<IRI>}`` therefore matches a node and every descendant (a serotype/genus/disease
    rollup), generalised here across MHC allele, organism, disease, antigen, and assay.
  * Column projection (``select=``), streaming output, and retry on 5xx/timeout.

Use as a library::

    import core
    n = core.count("tcell_search", core.subtree_filter("host_organism", "NCBITaxon:9606"))
    a = core.resolve_entity("allele", "HLA-A*02:01")        # -> iri + ancestors + serotype
    for row in core.iter_all("tcell_search", {"mhc_allele_name": "eq.HLA-A2"},
                             select="structure_id,linear_sequence"):
        ...

Stdlib only. No third-party dependencies. Read-only.

Docs: https://query-api.iedb.org/docs/swagger/   (PostgREST)
"""

import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error

__all__ = [
    "BASE_URL", "PAGE_SIZE", "CURSOR_KEYS", "ENDPOINTS", "ONTOLOGY", "SEROTYPE_RE",
    "request", "count", "iter_all", "fetch_all", "fetch_one",
    "subtree_filter", "resolve_entity", "resolve_iri_name", "mhc_serotype", "write_rows",
]

BASE_URL = "https://query-api.iedb.org"   # bare root; /api/v1/ also resolves
PAGE_SIZE = 10000                          # IQ-API hard cap per response (verified)
MAX_RETRIES = 4
RETRY_BACKOFF = 3.0                         # seconds * attempt

# Serotype names: HLA-A2 / HLA-A24 / HLA-B7 (class I) or HLA-DR4 / HLA-DQ8 / HLA-DP4
# (class II) -- locus letters directly followed by a number, NO chain designator. Gene/locus
# nodes (HLA-A, HLA-DRB1) are deliberately NOT matched. Allele-specific (MHC) convenience.
SEROTYPE_RE = re.compile(r"^HLA-(?:[A-CEFG]\d+|D[PQR]\d+)$")

# Integer/string cursor (primary-key) column per main search endpoint, for keyset pagination.
CURSOR_KEYS = {
    "tcell_search": "tcell_id",
    "mhc_search": "elution_id",
    "bcell_search": "bcell_id",
    "tcr_search": "receptor_group_id",
    "bcr_search": "receptor_group_id",
    "epitope_search": "structure_id",
    "reference_search": "reference_id",
    "antigen_search": "parent_source_antigen_iri",   # string cursor (lexicographic keyset)
}

# Endpoint grain, for the `endpoints` subcommand / discoverability.
ENDPOINTS = {
    "tcell_search": "one row per T-cell assay",
    "mhc_search": "one row per MHC ligand / binding assay",
    "bcell_search": "one row per B-cell assay",
    "tcr_search": "one row per TCR (receptor group)",
    "bcr_search": "one row per BCR (receptor group)",
    "epitope_search": "one row per epitope (aggregated across assays)",
    "antigen_search": "one row per source antigen",
    "reference_search": "one row per reference / publication",
}

# ontology kind -> (subtree_array_field, scalar_iri_field|None, scalar_name_field)
# scalar IRI present  -> name<->IRI resolution supported.
# scalar IRI is None  -> subtree FILTERING works, but reverse IRI->name is unsupported.
ONTOLOGY = {
    "allele":          ("mhc_allele_iri_search",            "mhc_allele_iri",             "mhc_allele_name"),
    "host_organism":   ("host_organism_iri_search",         "host_organism_iri",          "host_organism_name"),
    "source_organism": ("source_organism_iri_search",       "source_organism_iri",        "source_organism_name"),
    "antigen":         ("parent_source_antigen_iri_search", "parent_source_antigen_iri",  "parent_source_antigen_name"),
    "disease":         ("disease_iri_search",               None,                         "disease_names"),
    "assay":           ("assay_iri_search",                 None,                         "assay_names"),
}


# --------------------------------------------------------------------------- #
# Low-level HTTP
# --------------------------------------------------------------------------- #
def _request(url, headers, timeout):
    """Single HTTP GET with retries; returns (response_bytes, headers_dict)."""
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read(), dict(resp.headers)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            code = getattr(exc, "code", None)
            last_err = exc
            # 4xx (except transient 408/429) are not worth retrying
            if code is not None and code < 500 and code not in (408, 429):
                raise
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF * attempt)
    raise RuntimeError(f"GET failed after {MAX_RETRIES} attempts: {url}\n{last_err}")


def _build_url(endpoint, params):
    # safe= keeps PostgREST syntax readable & avoids double-encoding: * : { } . , ( )
    qs = urllib.parse.urlencode(params, safe="*:{}.,()")
    return f"{BASE_URL}/{endpoint}?{qs}"


def _header(headers, name):
    name = name.lower()
    for k, v in headers.items():
        if k.lower() == name:
            return v
    return None


def request(endpoint, params, timeout=90, accept="application/json"):
    """Raw GET against an endpoint. Returns (body_bytes, headers_dict)."""
    return _request(_build_url(endpoint, params), {"Accept": accept}, timeout)


def _cursor_key(endpoint, order_key=None):
    if order_key:
        return order_key
    if endpoint in CURSOR_KEYS:
        return CURSOR_KEYS[endpoint]
    raise ValueError(
        f"No known cursor key for endpoint {endpoint!r}; pass order_key=<unique sortable column>.")


# --------------------------------------------------------------------------- #
# Counting & pagination
# --------------------------------------------------------------------------- #
def count(endpoint, filters=None, mode="estimated", timeout=90):
    """Total matching rows via the Content-Range header. mode: 'estimated' (cheap) or 'exact'.

    Returns the integer total, or -1 if the server did not report one.
    """
    params = dict(filters or {})
    ck = CURSOR_KEYS.get(endpoint)
    if ck:
        params.setdefault("select", ck)   # minimise payload
    params["limit"] = 1
    _body, hdrs = _request(
        _build_url(endpoint, params),
        {"Accept": "application/json", "Prefer": f"count={mode}"},
        timeout)
    cr = _header(hdrs, "Content-Range") or ""
    tail = cr.split("/")[-1].strip() if "/" in cr else ""
    return int(tail) if tail.isdigit() else -1


def iter_all(endpoint, filters=None, select=None, order_key=None,
             timeout=90, page_size=PAGE_SIZE, max_rows=None, label=""):
    """Yield every matching row using keyset (cursor) pagination.

    O(1) per page: order by a unique sortable key and advance with `<key>=gt.<last>`. Coverage is
    complete and duplicate-free. If `select` is given but omits the cursor column, the cursor is
    added for paging and stripped from yielded rows. Stops on a short page or `max_rows`.
    """
    order_key = _cursor_key(endpoint, order_key)
    user_cols = [c.strip() for c in (select or "").split(",") if c.strip()]
    drop_cursor = bool(user_cols) and order_key not in user_cols
    sel_cols = list(user_cols)
    if drop_cursor:
        sel_cols.append(order_key)
    select_str = ",".join(sel_cols) if sel_cols else None

    cursor = None
    yielded = 0
    while True:
        params = dict(filters or {})
        if select_str:
            params["select"] = select_str
        params["order"] = order_key
        params["limit"] = page_size
        if cursor is not None:
            params[order_key] = f"gt.{cursor}"
        body, _hdrs = _request(_build_url(endpoint, params), {"Accept": "application/json"}, timeout)
        page = json.loads(body)
        if not page:
            break
        last = page[-1][order_key]            # capture before any pop
        for row in page:
            if drop_cursor:
                row.pop(order_key, None)
            yield row
            yielded += 1
            if max_rows is not None and yielded >= max_rows:
                sys.stderr.write("\n")
                return
        cursor = last
        sys.stderr.write(f"    [{label or endpoint}] {yielded:,} rows\r")
        sys.stderr.flush()
        if len(page) < page_size:
            break
    sys.stderr.write("\n")


def fetch_all(endpoint, filters=None, select=None, order_key=None, timeout=90, max_rows=None, label=""):
    """Materialise iter_all() into a list."""
    return list(iter_all(endpoint, filters, select, order_key,
                         timeout=timeout, max_rows=max_rows, label=label))


def fetch_one(endpoint, filters=None, select=None, timeout=90):
    params = dict(filters or {})
    if select:
        params["select"] = select
    params["limit"] = 1
    body, _hdrs = _request(_build_url(endpoint, params), {"Accept": "application/json"}, timeout)
    rows = json.loads(body)
    return rows[0] if rows else None


# --------------------------------------------------------------------------- #
# Ontology helpers (generalised MRO/NCBITaxon/DOID... subtree matching)
# --------------------------------------------------------------------------- #
def subtree_filter(kind, iri):
    """Filter dict matching all records annotated at `iri` OR any descendant of it.

    Works for every ontology `kind` (the array column always carries ancestor IRIs).
    e.g. subtree_filter("host_organism", "NCBITaxon:9606") -> {"host_organism_iri_search": "cs.{NCBITaxon:9606}"}
    """
    if kind not in ONTOLOGY:
        raise ValueError(f"Unknown ontology kind {kind!r}; choose from {sorted(ONTOLOGY)}")
    array_field = ONTOLOGY[kind][0]
    return {array_field: f"cs.{{{iri}}}"}


_iri_name_cache = {}


def resolve_iri_name(kind, iri, timeout=90, endpoints=("mhc_search", "tcell_search")):
    """Map an ontology IRI to the name used for records annotated at that node.

    Returns None for kinds without a scalar IRI column (disease, assay) or if not found.
    """
    array_field, iri_field, name_field = ONTOLOGY[kind]
    if iri_field is None:
        return None
    key = (kind, iri)
    if key in _iri_name_cache:
        return _iri_name_cache[key]
    name = None
    for ep in endpoints:
        rec = fetch_one(ep, {iri_field: f"eq.{iri}"},
                        select=f"{iri_field},{name_field}", timeout=timeout)
        if rec and rec.get(name_field):
            name = rec[name_field]
            break
    _iri_name_cache[key] = name
    return name


def resolve_entity(kind, name, timeout=90, endpoints=("mhc_search", "tcell_search")):
    """Resolve a named entity to its IRI + ancestor chain (and, for alleles, class + serotype).

    Only supported for kinds with a scalar IRI column (allele, host_organism, source_organism,
    antigen). Ancestors sharing the entity's own CURIE prefix are name-resolved; foreign-ontology
    ancestors (e.g. GO terms above an allele) are listed with name=None and not queried.
    """
    array_field, iri_field, name_field = ONTOLOGY[kind]
    if iri_field is None:
        raise ValueError(
            f"kind {kind!r} has no scalar IRI column, so a name cannot be resolved to an IRI. "
            f"Once you know the IRI, use subtree_filter({kind!r}, <iri>) to filter.")
    sel = f"{iri_field},{array_field},{name_field}"
    if kind == "allele":
        sel += ",mhc_class"

    # 1) exact match on the canonical name.
    rec = None
    matched_name = name
    for ep in endpoints:
        rec = fetch_one(ep, {name_field: f"eq.{name}"}, select=sel, timeout=timeout)
        if rec:
            break

    # 2) forgiving fallback: case-insensitive substring (free-text names often carry a suffix,
    #    e.g. organism 'Homo sapiens (human)'). Skipped for alleles -- their IEDB names are
    #    canonical/exact, and an ilike scan over the huge mhc table would be needlessly slow.
    if not rec:
        if kind == "allele":
            raise ValueError(
                f"allele {name!r} not found. Use the exact IEDB allele name, e.g. 'HLA-A*02:01'.")
        cand = {}   # iri -> canonical name
        for ep in endpoints:
            for row in iter_all(ep, {name_field: f"ilike.*{name}*"},
                                select=f"{iri_field},{name_field}", timeout=timeout, max_rows=2000):
                if row.get(iri_field):
                    cand.setdefault(row[iri_field], row.get(name_field))
            if cand:
                break
        if not cand:
            raise ValueError(
                f"{kind} {name!r} not found (exact or substring) in {endpoints}. Check the IEDB "
                f"name; `schema {endpoints[0]}` shows columns, or count with an ilike filter.")
        if len(cand) == 1:
            only_iri = next(iter(cand))
        else:
            # prefer an exact (case-insensitive) name, or the canonical 'name (common)' form,
            # so 'Mus musculus' -> 'Mus musculus (mouse)'; otherwise refuse as ambiguous.
            q = name.strip().lower()
            preferred = [iri for iri, nm in cand.items()
                         if nm and (nm.strip().lower() == q or nm.strip().lower().startswith(q + " ("))]
            if len(preferred) == 1:
                only_iri = preferred[0]
            else:
                shown = "; ".join(f"{nm} [{iri}]" for iri, nm in list(cand.items())[:8])
                raise ValueError(
                    f"{kind} {name!r} is ambiguous ({len(cand)} matches): {shown}"
                    + (" ..." if len(cand) > 8 else "") + ". Re-run with the exact name.")
        matched_name = cand[only_iri]
        for ep in endpoints:
            rec = fetch_one(ep, {iri_field: f"eq.{only_iri}"}, select=sel, timeout=timeout)
            if rec:
                break

    iri = rec.get(iri_field)
    prefix = iri.split(":")[0] if iri else None
    ancestors = []
    for anc in (rec.get(array_field) or []):
        if not isinstance(anc, str) or anc == iri:
            continue
        # resolve names only within the same ontology (efficiency + clean output)
        anc_name = resolve_iri_name(kind, anc, timeout=timeout) if (prefix and anc.startswith(prefix)) else None
        ancestors.append({"iri": anc, "name": anc_name})

    out = {"kind": kind, "name": matched_name, "query": name, "iri": iri, "ancestors": ancestors}
    if matched_name != name:
        out["note"] = f"matched {matched_name!r} for query {name!r}"
    if kind == "allele":
        out["mhc_class"] = rec.get("mhc_class")
        s_name, s_iri = mhc_serotype(ancestors)
        out["serotype_name"], out["serotype_iri"] = s_name, s_iri
    return out


def mhc_serotype(ancestors):
    """(name, iri) of the first ancestor whose name matches the HLA serotype pattern, else (None, None)."""
    for a in ancestors:
        nm = a.get("name")
        if nm and SEROTYPE_RE.match(nm):
            return nm, a.get("iri")
    return None, None


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _flat(v):
    if v is None:
        return ""
    if isinstance(v, list):
        return "; ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, default=str)
    return v


def write_rows(rows, fmt, out_fh):
    """Stream rows to a file handle as CSV (header from the first row) or a JSON array."""
    if fmt == "csv":
        writer = None
        for row in rows:
            if writer is None:
                writer = csv.DictWriter(out_fh, fieldnames=list(row.keys()), extrasaction="ignore")
                writer.writeheader()
            writer.writerow({k: _flat(v) for k, v in row.items()})
    elif fmt == "json":
        out_fh.write("[")
        first = True
        for row in rows:
            out_fh.write(("" if first else ",") + "\n  " + json.dumps(row, default=str))
            first = False
        out_fh.write("\n]\n" if not first else "]\n")
    else:
        raise ValueError(f"Unknown format {fmt!r} (use 'csv' or 'json')")
