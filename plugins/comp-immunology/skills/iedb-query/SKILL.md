---
name: iedb-query
description: Efficient general-purpose client and CLI for the IEDB Query API (IQ-API, https://query-api.iedb.org) — the public PostgREST interface to the Immune Epitope Database. Use whenever a user wants to pull, count, filter, or export IEDB data (T-cell assays, MHC-binding/ligand assays, B-cell assays, TCR/BCR receptors, epitopes, antigens, references) for any MHC allele, organism, disease, antigen, or assay type, or wants to query the IEDB efficiently / avoid slow or incomplete pulls / paginate a large IEDB result / count how many records match. Provides keyset (cursor) pagination (the IQ-API caps responses at 10k rows), cheap size estimation (count=estimated), CSV/JSON streaming export, and an ontology-subtree matcher: every search endpoint carries `<entity>_iri_search` arrays of an annotation's ancestor IRIs, so a single `cs.{IRI}` filter rolls a query up an ontology — generalized across MHC allele (MRO; e.g. HLA-A2 serotype), organism (NCBITaxon), disease (DOID), antigen, and assay. Importable as a Python module (`import iedb`) and runnable as a CLI (`python iedb.py endpoints|schema|count|query|resolve`). Stdlib-only, read-only, no auth. For the specific "T-cell-assay-but-no-MHC-binding gap" analysis use the separate iedb-assay-gap skill instead.
---

# iedb-query

A thin, efficient toolkit over the IEDB Query API (IQ-API), a public no-auth **PostgREST**
front end at `https://query-api.iedb.org`. Use it to count, query, paginate, and export IEDB
data correctly and fast, and to roll queries up an ontology (serotype / genus / disease family).

The script is **both** an importable library and a CLI: `scripts/iedb.py`.

## When to use

- Pull / filter / export IEDB records for an allele, organism, disease, antigen, or assay type.
- Count how many records match a filter (without waiting on a slow exact count).
- Paginate a large IEDB result set completely and efficiently.
- Resolve an entity name to its ontology IRI + ancestors, then query a subtree (rollup).

For the dedicated **"T-cell positive but no MHC-binding assay" gap** analysis with an HTML report,
use the sibling `iedb-assay-gap` skill instead — this skill is the general query layer.

## CLI

```bash
S=~/.claude/plugins/comp-immunology/skills/iedb-query/scripts/iedb.py

python "$S" endpoints                       # list endpoints, grain, cursor keys, ontology kinds
python "$S" schema tcell_search             # show one endpoint's columns + a sample row
python "$S" count  tcell_search --where 'mhc_allele_name=eq.HLA-A2'          # estimated (fast)
python "$S" count  tcell_search --where 'structure_id=eq.18719' --exact      # exact

# query: auto keyset-paginated + streamed; prints estimated size first; CSV or JSON
python "$S" query  tcell_search \
    --where 'host_organism_iri_search=cs.{NCBITaxon:9606}' \
    --where 'qualitative_measure=like.Positive*' \
    --select structure_id,linear_sequence,mhc_allele_name \
    --format csv --out human_pos.csv

# resolve a name -> IRI + ancestors (+ serotype for alleles); prints ready-to-use filter hints
python "$S" resolve allele 'HLA-A*02:01'
python "$S" resolve host_organism 'Homo sapiens'      # forgiving: matches 'Homo sapiens (human)'
```

`--where` takes raw PostgREST `col=op.value` tokens (repeatable, AND-combined), so the full
operator set is available with no bespoke syntax: `eq neq gt gte lt lte like ilike in cs ov not`
and `or=(a.eq.1,b.eq.2)`. `query` always prints an estimated row count to stderr first and never
truncates silently (`--max-rows` caps are logged).

## Library

```python
import iedb
n  = iedb.count("mhc_search", iedb.subtree_filter("allele", "MRO:0001007"))   # estimated
a  = iedb.resolve_entity("allele", "HLA-A*02:01")     # {iri, ancestors, mhc_class, serotype_*}
rows = iedb.fetch_all("tcell_search", {"mhc_allele_name": "eq.HLA-A2"},
                      select="structure_id,linear_sequence")          # keyset-paginated list
for row in iedb.iter_all("epitope_search", {"structure_id": "gt.0"}): ...   # streaming generator
with open("out.csv","w",newline="") as fh:
    iedb.write_rows(iedb.iter_all("reference_search", {"reference_id":"lt.1000020"}), "csv", fh)
```

Key functions: `count`, `iter_all`/`fetch_all`, `fetch_one`, `subtree_filter`, `resolve_entity`,
`resolve_iri_name`, `mhc_serotype`, `write_rows`. See `references/iq_api_reference.md`.

## The ontology subtree trick (the efficient way to roll up)

Every search endpoint carries `<entity>_iri_search` array columns holding an annotation's IRI
**plus all its ancestor IRIs** in the relevant ontology. So:

- `mhc_allele_iri_search=cs.{MRO:0001007}` → `HLA-A*02:01` **and** every higher-resolution sub-allele.
- `host_organism_iri_search=cs.{NCBITaxon:9606}` → human **and** every descendant taxon.
- `disease_iri_search=cs.{DOID:9352}` → type 2 diabetes **and** its subtypes.

`resolve_entity(kind, name)` finds the IRI + ancestor chain so you know what to roll up to (and,
for alleles, derives the HLA serotype, e.g. `HLA-A2` = `MRO:0001530`). Supported kinds with
name↔IRI resolution: `allele`, `host_organism`, `source_organism`, `antigen`. `disease` and
`assay` support subtree **filtering** but not reverse IRI→name lookup (no scalar IRI column).

## Efficiency rules (baked in)

1. Responses cap at **10,000 rows** → always paginate. Use **keyset** (`order=<pk>` +
   `<pk>=gt.<last>`), never OFFSET (O(offset) rescans, degrades badly). `iter_all` does this.
2. `count=estimated` to gauge size (cheap); `count=exact` only when you need the precise total.
3. `--select` / `select=` to project only needed columns (smaller payloads).
4. Native CSV via `Accept: text/csv`; this skill streams pages so memory stays flat.
5. Retries on 5xx/timeout with backoff are built in.

## Gotchas

- URL-encode in raw curl: `*`→`%2A`, `:`→`%3A`, `{`→`%7B`, `}`→`%7D` (the client handles this).
- `assay_names` is a `|`-separated string, not a JSON array; `disease_names` is an array.
- `antigen_search`'s cursor (`parent_source_antigen_iri`) is a string (lexicographic keyset).
- `resolve_entity` matches **allele** names exactly (canonical, e.g. `HLA-A*02:01`); for
  organism/antigen names it falls back to a bounded case-insensitive substring scan, preferring
  an exact or `name (common)` match (so `Mus musculus` → `Mus musculus (mouse)`) and listing
  candidates only when genuinely ambiguous.
- A "no record" result means no matching annotation in IEDB, not biological absence.
- Resolving a **real** entity is fast (the lookup stops at the first match), but resolving a
  *misspelled* name is slow (~tens of seconds): confirming a name's absence requires scanning the
  unindexed name column. The common path (valid names) is fast; typos just take a moment to fail.
