# IQ-API notes for `epitopes-by-source`

The IEDB Query API (IQ-API, `https://query-api.iedb.org`) is a public, no-auth
**PostgREST** front end over the IEDB Postgres tables. Each table is an endpoint;
filters use PostgREST `col=op.value` syntax; results are JSON (or CSV via
`Accept: text/csv`). Docs: `https://query-api.iedb.org/docs/swagger/`.

## Endpoint used: `epitope_search`

- **Grain:** one row per epitope, **aggregated across all of that epitope's assays**.
- **Primary key / keyset cursor:** `structure_id` (integer, ascending). Web page:
  `https://www.iedb.org/epitope/<structure_id>`.
- Most informative columns are **list-valued** (an epitope can have many assays, hosts,
  MHC restrictions, references): `parent_source_antigen_names`, `source_organism_names`,
  `host_organism_names`, `mhc_allele_names`, `mhc_classes`, `qualitative_measures`,
  `assay_names`, `reference_ids`, `pubmed_ids`, `journal_names`, `tcell_ids`,
  `bcell_ids`, `elution_ids`. Scalars: `linear_sequence`, `linear_sequence_length`,
  `structure_type`.
- `assay_names` list elements are themselves `|`-joined (`"<assay>|<method>"`).
- **`epitope_search` has no scalar `source_organism_iri` / `source_organism_name`** —
  only the list `source_organism_names` and the array `source_organism_iri_search`.
  So name↔IRI resolution must run against `tcell_search` / `mhc_search`, which DO carry
  the scalar columns. (That's what `resolve_taxon` / `_iri_to_name` / `_substring_scan`
  in the script do.)

## The ontology subtree trick (how source filtering works)

Every search endpoint carries `<entity>_iri_search` array columns holding an
annotation's IRI **plus all its ancestor IRIs**. So:

```
source_organism_iri_search=cs.{NCBITaxon:10376}
   -> Human herpesvirus 4 (Epstein Barr virus) AND every descendant taxon (all strains).
host_organism_iri_search=cs.{NCBITaxon:9606}
   -> Homo sapiens AND every descendant taxon.
```

`cs.{IRI}` is the PostgREST "array contains" operator. A species-level rollup therefore
sweeps up all strains/isolates with one filter. Verified: EBV `NCBITaxon:10376` →
**6,134 epitopes** (12 distinct source taxa, incl. `strain B95-8` = `NCBITaxon:10377`).

## Efficiency rules (baked into the script)

1. Responses cap at **10,000 rows** → always paginate. Use **keyset** (`order=structure_id`
   + `structure_id=gt.<last>`), never OFFSET (O(offset) rescans). `iter_all` does this.
2. `Prefer: count=estimated` (cheap planner estimate) for the size preflight; the report's
   card counts are the exact fetched/after-filter numbers.
3. `select=` projects only the ~20 columns the report needs (smaller payloads).
4. Retries on 5xx / 408 / 429 / timeout with backoff are built in.

## Gotchas

- URL-encode in raw curl: `*`→`%2A`, `:`→`%3A`, `{`→`%7B`, `}`→`%7D`. The client passes
  `safe="*:{}.,()"` to `urlencode` so PostgREST syntax stays readable.
- A two-bound length filter uses PostgREST's `and=(...)`:
  `and=(linear_sequence_length.gte.8,linear_sequence_length.lte.11)`.
- The list-union columns (`mhc_classes`, `qualitative_measures`) describe an epitope's
  whole assay history; `--mhc-class` / `--positive-only` are applied **after** fetch so
  pre/post-filter counts are exact and explainable.
- A "no record" result = no matching annotation in IEDB, not biological absence.
- Resolving a valid organism is fast (stops at first match); confirming a **misspelled**
  name is slow (scans the unindexed name column before failing).

## Common NCBITaxon rollups (the alias map)

| Alias | NCBITaxon | Organism |
|---|---|---|
| ebv | 10376 | Human herpesvirus 4 (Epstein Barr virus) |
| cmv / hcmv | 10359 | Human herpesvirus 5 (Cytomegalovirus) |
| hiv / hiv-1 | 11676 | Human immunodeficiency virus 1 |
| hbv | 10407 | Hepatitis B virus |
| hcv | 11103 | Hepatitis C virus |
| sars-cov-2 / sars2 / covid | 2697049 | SARS coronavirus 2 |
| influenza a / flu | 11320 | Influenza A virus |
| dengue | 12637 | Dengue virus |
| mtb / m. tuberculosis | 1773 | Mycobacterium tuberculosis |
| plasmodium falciparum | 5833 | Plasmodium falciparum |

The resolver re-derives the canonical IEDB name from the API for whatever IRI it uses, so
a wrong alias IRI surfaces immediately as a name mismatch. Only EBV (10376) is verified
against the live API; the others are standard NCBITaxon species ids and should be
spot-checked when first used.
