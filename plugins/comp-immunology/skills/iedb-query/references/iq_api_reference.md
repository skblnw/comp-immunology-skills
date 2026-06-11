# IEDB Query API (IQ-API) — reference for the `iedb-query` skill

Overview: <https://help.iedb.org/hc/en-us/articles/4402872882189>
Swagger: <https://query-api.iedb.org/docs/swagger/> · Use-cases: <https://github.com/IEDB/IQ-API-use-cases>

The IQ-API is a public, no-auth **PostgREST** front end over the IEDB Postgres tables.
Base URL `https://query-api.iedb.org` (the bare root; `https://query-api.iedb.org/api/v1/...`
also resolves). Each table/view is an endpoint; filters use PostgREST syntax; results are JSON
(default) or CSV (`Accept: text/csv`). TSV is **not** supported.

## Endpoints

39 endpoints total. The main **search** endpoints (one row per the stated grain) and their
integer/string **cursor key** (unique sortable column, used for keyset pagination):

| Endpoint | Grain | Cursor key |
|---|---|---|
| `tcell_search` | one T-cell assay | `tcell_id` |
| `mhc_search` | one MHC ligand / binding assay | `elution_id` |
| `bcell_search` | one B-cell assay | `bcell_id` |
| `tcr_search` | one TCR (receptor group) | `receptor_group_id` |
| `bcr_search` | one BCR (receptor group) | `receptor_group_id` |
| `epitope_search` | one epitope (aggregated across assays) | `structure_id` |
| `antigen_search` | one source antigen | `parent_source_antigen_iri` (string) |
| `reference_search` | one reference / publication | `reference_id` |

`structure_id` is the integer epitope key shared across tables; web page
`https://www.iedb.org/epitope/<structure_id>`. Also present: `*_export` views, `*_to_*` link
tables (e.g. `epitope_to_mhc`, `tcr_to_epitope`), `curie_map` (CURIE-prefix → URL template, for
building external links — **not** an IRI→label resolver), `api_metrics`.

## The ontology subtree trick (the core efficiency pattern)

Every search endpoint exposes the **same nine `*_iri_search` array columns**, each holding an
annotation's own IRI **plus all its ancestor IRIs** in the relevant ontology:

```
host_organism_iri_search   disease_iri_search             mhc_allele_iri_search
assay_iri_search           parent_source_antigen_iri_search
source_organism_iri_search non_peptidic_molecule_iri_search
r_object_source_molecule_iri_search  r_object_source_organism_iri_search
```

So `<field>=cs.{<IRI>}` (array **contains**) matches a node **and every descendant** — an
ontology rollup in one filter:

```bash
# human host + every descendant taxon (NCBITaxon)
curl -s 'https://query-api.iedb.org/tcell_search?host_organism_iri_search=cs.{NCBITaxon:9606}&select=structure_id&limit=1' -H 'Prefer: count=estimated' -D - -o /dev/null | grep -i content-range
# HLA-A*02:01 + higher-resolution sub-alleles (MRO)
curl -s 'https://query-api.iedb.org/mhc_search?mhc_allele_iri_search=cs.{MRO:0001007}&select=structure_id&limit=1' -H 'Prefer: count=estimated' -D - -o /dev/null | grep -i content-range
```

### `ONTOLOGY` registry (verified columns)

| kind | subtree array (filter) | scalar IRI | scalar name | resolve name↔IRI? |
|---|---|---|---|---|
| `allele` | `mhc_allele_iri_search` | `mhc_allele_iri` | `mhc_allele_name` | yes |
| `host_organism` | `host_organism_iri_search` | `host_organism_iri` | `host_organism_name` | yes |
| `source_organism` | `source_organism_iri_search` | `source_organism_iri` | `source_organism_name` | yes |
| `antigen` | `parent_source_antigen_iri_search` | `parent_source_antigen_iri` | `parent_source_antigen_name` | yes |
| `disease` | `disease_iri_search` | *(none)* | `disease_names` | filter only |
| `assay` | `assay_iri_search` | *(none)* | `assay_names` | filter only |

`disease` and `assay` have no scalar IRI column, so subtree **filtering** works but reverse
IRI→name lookup does not. `resolve_entity` works for the four scalar-pair kinds; for alleles it
also derives the HLA **serotype** ancestor via `SEROTYPE_RE = ^HLA-(?:[A-CEFG]\d+|D[PQR]\d+)$`
(matches `HLA-A2`/`HLA-DR4`, not loci `HLA-A`/genes `HLA-DRB1`). Verified:
`HLA-A*02:01 → HLA-A2 (MRO:0001530)`, `HLA-A*24:02 → HLA-A24 (MRO:0001533)`,
`HLA-DRB1*04:01 → HLA-DR4 (MRO:0001646)`.

## PostgREST operators

| Filter | Meaning | Example |
|---|---|---|
| `eq` `neq` | (in)equality | `mhc_allele_name=eq.HLA-A2` |
| `gt` `gte` `lt` `lte` | comparisons / keyset cursor | `tcell_id=gt.42` |
| `like` `ilike` | wildcard (`*`=`%`), case (in)sensitive | `mhc_allele_name=like.HLA-A*` |
| `in` | in list | `structure_id=in.(1,2,3)` |
| `cs` | array **contains** (the subtree trick) | `disease_iri_search=cs.{DOID:9352}` |
| `cd` `ov` | array contained-in / overlaps | `host_organism_iri_search=ov.{NCBITaxon:9606,NCBITaxon:10090}` |
| `not` | negate any op | `qualitative_measure=not.like.Negative*` |
| `or` `and` | logical groups | `or=(mhc_allele_name.eq.HLA-A2,mhc_allele_name.eq.HLA-A24)` |
| `select` | project columns | `select=structure_id,linear_sequence` |
| `order` | sort | `order=tcell_id` |

No full-text-search (`fts`/`plfts`) operators are exposed; use `like`/`ilike` on text columns.
Resource embedding works for `*_summary` fields and `*_to_*` link tables only.

URL-encoding (raw curl): `*`→`%2A`, `:`→`%3A`, `{`→`%7B`, `}`→`%7D`. The client uses
`urlencode(..., safe="*:{}.,()")` so these stay readable.

## Pagination & counting

- A single response is **capped at 10,000 rows** (verified: `limit=20000` still returns 9999 in
  `Content-Range`). You must paginate.
- **Keyset (cursor) pagination** — `order=<cursor>` + `<cursor>=gt.<last_value>` — is O(1) per
  page (the DB seeks via index). **Prefer this over `offset=`**, which rescans all skipped rows
  and degrades badly at high offsets. `iter_all`/`fetch_all` implement keyset.
- **Counting**: `Prefer: count=exact` returns the precise total after `/` in the `Content-Range`
  header but is slow on large tables; `Prefer: count=estimated` returns the PostgreSQL-planner
  estimate ~40–57% faster — use it to gauge size before a pull. `count()` defaults to estimated.
- `Range: <start>-<end>` headers also work as an alternative to `limit`/`offset`.

## Efficiency checklist

1. Paginate with keyset (never OFFSET). 2. Gauge size with `count=estimated`. 3. Project with
`select=`. 4. Stream CSV (`Accept: text/csv`) for big exports. 5. Roll up with `cs.{IRI}` on a
`*_iri_search` array instead of enumerating names. 6. Retry 5xx/timeout with backoff.

## Field gotchas

- `assay_names` is a **`|`-separated string** (e.g. `'IFNg release|ELISPOT'`), not a JSON array.
- `disease_names` **is** a JSON array. (`_as_items()` normalizes both.)
- `mhc_allele_name` can be a serotype (`HLA-A2`), a locus/molecule (`mouse`, `H2-K`), or a full
  allele (`HLA-A*02:01`) depending on the record's annotation resolution.
- `qualitative_measure`: `Positive` / `Positive-High` / `Positive-Intermediate` / `Positive-Low`
  / `Negative` — match positives with `like.Positive*`.

## Worked examples (curl / library / CLI)

```bash
# Count human-host T-cell assays (estimated)
curl -s 'https://query-api.iedb.org/tcell_search?host_organism_iri_search=cs.{NCBITaxon:9606}&select=tcell_id&limit=1' \
     -H 'Prefer: count=estimated' -D - -o /dev/null | grep -i content-range
```
```python
import iedb
iedb.count("tcell_search", iedb.subtree_filter("host_organism", "NCBITaxon:9606"))
iedb.resolve_entity("allele", "HLA-A*02:01")          # -> iri, ancestors, serotype
iedb.fetch_all("epitope_search", {"structure_id": "in.(18719,140621)"})
```
```bash
S=~/.claude/plugins/comp-immunology/skills/iedb-query/scripts/iedb.py
python "$S" count tcell_search --where 'host_organism_iri_search=cs.{NCBITaxon:9606}'
python "$S" query mhc_search   --where 'mhc_allele_iri_search=cs.{MRO:0001007}' \
        --select structure_id --format csv --out a0201_ligands.csv
python "$S" resolve allele 'HLA-A*02:01'
```
