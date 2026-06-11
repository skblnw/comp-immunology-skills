---
name: epitopes-by-source
description: Pull EVERY epitope in the IEDB whose SOURCE ORGANISM is a given pathogen / taxon and build a self-contained HTML report (plus CSV + JSON). Use whenever a user wants the epitopes from a source organism — e.g. "epitopes from EBV", "all IEDB epitopes for SARS-CoV-2", "show me every dengue epitope", "IEDB epitopes by source organism", "what epitopes does IEDB have for M. tuberculosis", "give me a report of CMV epitopes". Resolves common abbreviations (EBV, CMV, HIV-1, HBV, HCV, SARS-CoV-2, influenza A, dengue, M. tuberculosis, P. falciparum), bare NCBITaxon ids ('NCBITaxon:10376' or '10376'), or a substring of the organism's IEDB name, and rolls the query up the NCBITaxon subtree so ALL strains/descendants are included. Computes breakdowns — linear vs discontinuous, T-cell / B-cell / MHC-ligand assay coverage, positive-assay count, MHC class distribution, top source antigens / host organisms / MHC alleles, peptide-length histogram, and per-strain counts — and renders a sortable epitope table. Optional filters: host organism, positive-only, MHC class, linear-only, peptide-length range. Queries the IEDB Query API (IQ-API) epitope_search endpoint. Stdlib-only, read-only, no dependencies.
---

# Epitopes by source organism

Given a **source organism** (a pathogen or any taxon), return **every epitope in the IEDB that comes from it**, summarised as a self-contained HTML report plus a CSV and JSON of the full set. This is the "what does IEDB know about pathogen X" view: how many epitopes, what proteins they come from, which MHC classes/alleles restrict them, how many have positive assays, and the peptide-length profile.

## When to use this skill

Invoke whenever the user names a pathogen/organism and wants its epitopes:

- "Fetch the epitopes for EBV / Epstein-Barr virus."
- "Make a report of all IEDB epitopes from SARS-CoV-2."
- "How many dengue epitopes are in IEDB, and from which proteins?"
- "Give me the human-restricted, MHC class I, positive epitopes for influenza A."
- "List every epitope whose source organism is *Mycobacterium tuberculosis*."

The unit of analysis is an **IEDB epitope** (`epitope_search.structure_id`, one row per epitope aggregated across all of its assays), **not** a peptide string and **not** an individual assay. For the "T-cell positive but no MHC-binding assay" gap analysis, use the sibling `iedb-assay-gap` skill; for general ad-hoc IQ-API querying, use `iedb-query`.

## The unit of analysis

The IQ-API `epitope_search` endpoint returns one row per epitope. "Source organism" is the taxon the epitope is derived from. Matching is done by the **NCBITaxon subtree** the IQ-API exposes through the `source_organism_iri_search` array column (every record carries its taxon IRI plus all ancestor IRIs):

```
source_organism_iri_search=cs.{NCBITaxon:10376}
   -> Human herpesvirus 4 (Epstein Barr virus)  AND every EBV strain/descendant.
```

So a species-level query automatically includes all of that species' strains.

## Source resolution (3-tier, forgiving)

The `source` argument is resolved in three tiers:

1. **NCBITaxon id** — `NCBITaxon:10376` or bare `10376` → used directly (canonical name fetched from the API for display).
2. **Built-in alias** — common abbreviations map to a species-level rollup: `ebv`, `cmv`/`hcmv`, `hiv`/`hiv-1`, `hbv`, `hcv`, `sars-cov-2`/`sars2`/`covid`, `influenza a`/`flu`, `dengue`, `mtb`/`m. tuberculosis`, `plasmodium falciparum`, …
3. **Substring scan** — a bounded case-insensitive scan of the IEDB organism-name column (`tcell_search`/`mhc_search`). Prefers an exact / `"Name (common)"` match; if a name is genuinely **ambiguous** it fails fast and lists the candidate taxa + IRIs.

Note: "EBV" is **not** a substring of the IEDB name "Human herpesvirus 4 (Epstein Barr virus)" — that's exactly why the alias map exists. For organisms without an alias, pass the NCBITaxon id or part of the full name.

## How to invoke the bundled script

Stdlib-only. Run from the skill directory:

```bash
python scripts/epitopes_by_source.py EBV -o out/
python scripts/epitopes_by_source.py "SARS-CoV-2" --positive-only --mhc-class I -o out/
python scripts/epitopes_by_source.py dengue --linear-only --min-length 8 --max-length 11 --host human -o out/
python scripts/epitopes_by_source.py NCBITaxon:1773 -o out/        # M. tuberculosis by id
```

### Options

| Flag | Default | Purpose |
|---|---|---|
| `source` | (required) | Source organism: abbreviation, NCBITaxon id (`NCBITaxon:10376`/`10376`), or a name substring. |
| `--out`, `-o` | `.` | Output directory (created if missing). |
| `--max-rows` | (none) | Cap epitopes fetched (logs + flags the report when it truncates). |
| `--timeout` | `90` | Per-request timeout (s). |
| `--top-n` | `15` | Top-N length for the antigen / host / allele breakdowns. |
| `--table-cap` | `1000` | Max epitope rows rendered in the HTML. The **full** set is always in CSV/JSON. |
| `--host` | (none) | Restrict to a host organism (name or NCBITaxon id), e.g. `human`. *(server-side subtree)* |
| `--positive-only` | off | Keep only epitopes with ≥1 positive assay outcome. *(post-filter)* |
| `--with-mhc` | off | Keep only epitopes with MHC data — drops the "no MHC data" bucket (antibody/B-cell epitopes with no MHC restriction). *(post-filter)* |
| `--mhc-class {I,II}` | (none) | Keep only epitopes tested against the given MHC class. *(post-filter)* |
| `--allele NAME` | (none) | Keep only epitopes restricted by a specific MHC allele (e.g. `'HLA-A*02:01'`). Matches the allele + higher-resolution sub-typings (`HLA-A*02:01:01`); excludes the broader serotype (`HLA-A2`). *(post-filter)* |
| `--linear-only` | off | Keep only linear peptide epitopes (drop discontinuous). *(post-filter)* |
| `--min-length` / `--max-length` | (none) | Peptide-length bounds for linear epitopes. *(server-side)* |

Source and `--host` filter on the server via the NCBITaxon subtree (cheap); length bounds are server-side numeric filters. The list-union filters (`--positive-only`, `--mhc-class`) and `--linear-only` are applied after fetch so the report can state exact pre- vs post-filter counts.

### Pipeline

1. **Resolve** the source organism (and `--host` if given) to an NCBITaxon IRI.
2. **Count** (estimated) for a size preflight + the `--max-rows` warning.
3. **Fetch** all matching epitopes from `epitope_search` via keyset pagination (`order=structure_id` + `structure_id=gt.<last>`; pages past the 10k/response cap automatically).
4. **Classify** each epitope (linear/discontinuous, assay-type presence, positivity, MHC class, antigen, strain, host, alleles, refs).
5. **Post-filter** (positive / MHC class / linear) and **compute breakdowns**.
6. **Render** `report.html` + write `epitopes_<IRI>.csv` / `.json`.

## Output schema

Everything is written inside `--out`:

```
<out>/
├── report.html                    # self-contained: cards + bar-chart breakdowns + sortable table
├── epitopes_NCBITaxon_10376.csv   # one row per epitope (full set)
└── epitopes_NCBITaxon_10376.json  # {metadata, breakdowns, epitopes[]}
```

**`epitopes_<IRI>.csv` columns:**
`structure_id, linear_sequence, length, structure_type, source_organism, parent_antigens, mhc_alleles, mhc_classes, host_organisms, qualitative_measures, assay_types, has_tcell, has_bcell, has_mhc_ligand, any_positive, n_references, pubmed_ids, journals, iedb_url`

- `source_organism` — the most specific source taxon name (e.g. the strain) for that epitope.
- `assay_types` — the distinct assay/method tokens recorded for the epitope (IEDB stores them `|`-joined).
- `iedb_url` — `https://www.iedb.org/epitope/<structure_id>`.

The HTML renders the top `--table-cap` epitopes (ranked by positive assay → #references → length; cap stated when hit); CSV/JSON hold the full set.

## Breakdowns

`total`; **linear vs discontinuous**; **with T-cell / B-cell / MHC-ligand assay data** (presence of `tcell_ids` / `bcell_ids` / `elution_ids`); **with ≥1 positive assay**; **MHC class distribution** (I / II / both / other / no-MHC-data); **top-N source antigens**; **top-N host organisms**; **top-N MHC alleles**; **peptide-length histogram** (linear only); **distinct source taxa** under the rollup with per-taxon counts. All counts are per epitope (a value appearing twice in one epitope is counted once for that breakdown).

## Examples

```bash
# EBV (alias -> NCBITaxon:10376), full report:
python scripts/epitopes_by_source.py EBV -o ebv_report/

# SARS-CoV-2 positive class I epitopes only:
python scripts/epitopes_by_source.py "SARS-CoV-2" --positive-only --mhc-class I -o cov2/

# Dengue linear 8-11mers assayed in human host:
python scripts/epitopes_by_source.py dengue --linear-only --min-length 8 --max-length 11 --host human -o dengue/

# EBV 9-mers restricted by HLA-A*02:01 (a focused MHC-I candidate list):
python scripts/epitopes_by_source.py EBV --allele "HLA-A*02:01" --min-length 9 --max-length 9 -o ebv_A0201_9mer/
```

## Workflow guidance for the assistant

1. Always run the script — don't hand-craft IQ-API URLs.
2. If the organism isn't a built-in alias and resolution is ambiguous, the script prints the candidate taxa; pick the right NCBITaxon id and re-run.
3. After the run, point the user at `report.html` and lead with the headline: total epitopes, the positive-assay subset, the dominant source antigens, and the MHC class split.
4. Don't paste the HTML into chat — quote the headline numbers and a couple of standout epitopes/antigens.

## Gotchas

- "EBV" (and other abbreviations) are not substrings of the IEDB organism name — that's what the alias map is for. Unknown organisms: pass the NCBITaxon id or a name fragment.
- A null/zero result means no matching annotation in IEDB, **not** biological absence.
- Discontinuous / non-peptidic epitopes have no `linear_sequence`; they're shown with a "discontinuous" tag and excluded from the length histogram (still counted in the total).
- Large rollups (e.g. SARS-CoV-2, influenza A) can have very many epitopes; keyset pagination pages through automatically. Use `--max-rows` to cap (the report flags truncation).
- Estimated match counts are PostgreSQL-planner estimates; the fetched / after-filter counts shown in the report cards are exact.
- Resolving a real organism is fast; confirming a **misspelled** name is slow (it scans the unindexed name column before failing).

## Extended API patterns

For the IQ-API endpoint/operator reference and the ontology subtree approach used here, see `references/iq_api_notes.md`.
