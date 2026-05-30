# IEDB Query API (IQ-API) — notes for the assay-gap skill

Official overview: <https://help.iedb.org/hc/en-us/articles/4402872882189>
Interactive reference (Swagger): <https://query-api.iedb.org/docs/swagger/>
Use-cases repo: <https://github.com/IEDB/IQ-API-use-cases>

The IQ-API is a **PostgREST** front end over the IEDB Postgres tables. Public, no
auth, no account. Base URL: `https://query-api.iedb.org`. Each table is an
endpoint; queries use PostgREST filter syntax. Results are JSON (default) or TSV
(`Accept: text/csv`).

## Endpoints used here

| Endpoint | Grain | Key fields used |
|---|---|---|
| `tcell_search` | one row per **T-cell assay** | `structure_id`, `linear_sequence`, `qualitative_measure`, `mhc_allele_name`, `mhc_allele_resolution`, `mhc_allele_iri`, `mhc_allele_iri_search`, `mhc_class`, `source_organism_name`, `parent_source_antigen_name`, `disease_names`, `assay_names` |
| `mhc_search` | one row per **MHC ligand / binding assay** | `structure_id`, `mhc_allele_name`, `mhc_allele_resolution`, `mhc_allele_iri`, `mhc_allele_iri_search` |
| `epitope_search` | one row per **epitope** (aggregated) | not used for per-allele gap — its `mhc_allele_names` list is a union across all assays, so it can't isolate per-allele binding |

`structure_id` is the integer epitope-structure key shared across all tables. Web
page: `https://www.iedb.org/epitope/<structure_id>`.

Other endpoints (not used): `bcell_search`, `tcr_search`, `bcr_search`,
`antigen_search`, `reference_search`, plus link tables (`epitope_to_mhc`, …).

## PostgREST operators (relevant subset)

| Filter | Meaning | Example |
|---|---|---|
| `eq` | equals | `mhc_allele_name=eq.HLA-A*02:01` |
| `in` | in list | `mhc_allele_name=in.("HLA-A2","HLA-A*02:01")` |
| `like`/`ilike` | wildcard (`*` = `%`) | `mhc_allele_name=like.HLA-A*` |
| `cs` | array **contains** | `mhc_allele_iri_search=cs.{MRO:0001007}` |
| `cd` | array contained-in | |
| `select` | choose columns | `select=structure_id,mhc_allele_name` |
| `order` | sort | `order=tcell_id` |

URL-encoding gotchas: `*` → `%2A`, `:` → `%3A`, array braces `{}` and IRIs may be
sent literally. The script uses `urllib.parse.urlencode(..., safe="*:{}.,()")`.

## Pagination

A single response is capped at **10,000 rows**. To page:

- `Prefer: count=exact` makes the response include `Content-Range: 0-9999/39248`
  — the number after `/` is the exact total.
- Then loop with `offset=` / `limit=` (or `Range` headers) ordered by a stable key
  (`order=tcell_id`, `order=elution_id`) until `offset >= total`.

## The MRO ontology trick (core of this skill)

Every record carries `mhc_allele_iri` (its own allele node) **and**
`mhc_allele_iri_search` — an array of that node's IRI plus **all ancestor IRIs**
in the MHC Restriction Ontology (MRO), e.g. for `HLA-A*02:01`:

```
mhc_allele_iri        = MRO:0001007
mhc_allele_iri_search = [GO:0032991, GO:0042611, MRO:0000011, MRO:0001005 (HLA-A),
                         MRO:0001454 (human), MRO:0001530 (HLA-A2 serotype),
                         MRO:0001007 (self), ...]
```

Consequences exploited here:

- **Subtree-down match** (allele + higher-resolution typings):
  `mhc_allele_iri_search=cs.{<allele IRI>}`. A child (e.g. `HLA-A*02:01:01`) has the
  allele's IRI among its ancestors, so it matches; an ancestor (serotype) does not.
- **Serotype rollup** (go up to the serotype, then take records annotated *at* that
  node): find the serotype ancestor by resolving each ancestor IRI to a name and
  matching the serotype pattern `^HLA-(?:[A-CEFG]\d+|D[PQR]\d+)$` (locus letters +
  number, no chain letter — so `HLA-A2`/`HLA-DR4` match but loci `HLA-A`/genes
  `HLA-DRB1` do not). Then match `mhc_allele_name=eq.<serotype>`.
- Exact (`cs` on allele IRI) and serotype (`eq` on serotype name) record sets are
  **disjoint** → safe to union without double-counting.

Verified serotype derivations: `HLA-A*02:01 → HLA-A2 (MRO:0001530)`,
`HLA-A*24:02 → HLA-A24 (MRO:0001533)`, `HLA-DRB1*04:01 → HLA-DR4 (MRO:0001646)`.

## Worked queries

```bash
# exact total of T-cell assays restricted by HLA-A*02:01 (+ sub-alleles)
curl -s 'https://query-api.iedb.org/tcell_search?mhc_allele_iri_search=cs.{MRO:0001007}&select=structure_id&limit=1' \
     -H 'Prefer: count=exact' -D - -o /dev/null | grep -i content-range
# -> Content-Range: 0-0/39248

# serotype-level T-cell rows (annotated exactly "HLA-A2")
curl -s 'https://query-api.iedb.org/tcell_search?mhc_allele_name=eq.HLA-A2&select=structure_id&limit=1' \
     -H 'Prefer: count=exact' -D - -o /dev/null | grep -i content-range

# binding assays for the same allele
curl -s 'https://query-api.iedb.org/mhc_search?mhc_allele_iri_search=cs.{MRO:0001007}&select=structure_id&limit=1' \
     -H 'Prefer: count=exact' -D - -o /dev/null | grep -i content-range
```

The gap = distinct `structure_id` on the T-cell side minus distinct `structure_id`
on the binding side, over the matched allele set.
