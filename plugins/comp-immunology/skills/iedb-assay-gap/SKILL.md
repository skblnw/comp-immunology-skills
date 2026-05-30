---
name: iedb-assay-gap
description: For a given MHC allele, find epitopes in the IEDB that HAVE T-cell assays but have NO MHC-binding assay against that allele — the "binding-data gap": peptides proven T-cell antigenic in the allele's context whose direct MHC binding has never been experimentally measured. Use whenever a user gives one or more MHC alleles (e.g. HLA-A*02:01, HLA-DRB1*04:01) and asks which epitopes are T-cell positive / T-cell tested but lack binding data, wants binding-validation candidates, wants to quantify the gap between T-cell and MHC-binding coverage for an allele, or asks to compare T-cell vs MHC assay coverage. Queries the IEDB Query API (IQ-API) tcell_search + mhc_search endpoints; allele matching is grounded in the MHC Restriction Ontology (MRO) so a serotype rollup (e.g. HLA-A2 for HLA-A*02:01) is included and every match is tagged by resolution. Outputs a self-contained HTML report plus per-allele CSV and JSON. Stdlib-only, no dependencies.
---

# IEDB assay gap

Given an MHC allele, return the epitopes that have **T-cell assays but no MHC-binding assay** for that allele. These are peptides shown to be T-cell antigenic in the allele's context whose binding to the allele has never been directly measured in IEDB — high-value candidates for binding validation, and a quantifiable hole in the binding data.

## When to use this skill

Invoke whenever the user supplies one or more MHC allele names and asks any of:

- "Which epitopes have a T-cell assay for HLA-A\*02:01 but no MHC-binding assay?"
- "Find binding-validation candidates for HLA-DRB1\*04:01."
- "How big is the gap between T-cell and binding coverage for this allele?"
- "Give me peptides that are T-cell positive on HLA-A\*24:02 but were never tested for binding."
- "Compare T-cell vs MHC assay coverage for these alleles."

The unit of analysis is an **IEDB epitope (`structure_id`)**, not a peptide string. The allele must be a real IEDB allele name (e.g. `HLA-A*02:01`, `HLA-DRB1*04:01`, `H2-Kb`).

## The definition (read before reporting results)

For allele **X**, an epitope is a **candidate** when it has:

- **≥ 1 T-cell assay** whose MHC restriction matches X (`tcell_search`), AND
- **0 MHC-binding assays** whose tested allele matches X (`mhc_search`).

This is the **same-allele** gap. "No binding assay" means no *recorded* binding assay in IEDB — **not** experimental evidence of non-binding.

## Allele matching (ontology-grounded, auditable)

Matching is driven by the **MHC Restriction Ontology (MRO)**, which the IQ-API exposes through the `mhc_allele_iri_search` array on every record (the record's own allele IRI plus all ancestor IRIs).

| Resolution | How it matches | Example for `HLA-A*02:01` |
|---|---|---|
| `exact` | record's allele name == X | `HLA-A*02:01` |
| `sub-allele` | `mhc_allele_iri_search=cs.{<X IRI>}`, higher-resolution typings | `HLA-A*02:01:01` |
| `serotype` | exact-name match on the serotype ancestor derived from MRO | `HLA-A2` |

The serotype ancestor is derived from the ontology and **printed/shown for audit** (e.g. `HLA-A*02:01 → HLA-A2`, `HLA-DRB1*04:01 → HLA-DR4`). Serotype rollup is **on by default** (`--match-mode serotype`); pass `--match-mode exact` to restrict to the allele + sub-alleles only. The exact and serotype record sets are disjoint (a serotype node is an *ancestor* of the allele, so its IRI array never contains the allele's IRI), so there is no double-counting. Every candidate is tagged with the resolution(s) at which its T-cell evidence exists, and the report reports both an **exact-only** and a **serotype-inclusive** gap count.

## How to invoke the bundled script

Stdlib-only. Run from the skill directory:

```bash
python scripts/iedb_assay_gap.py --allele "HLA-A*02:01" -o out/
python scripts/iedb_assay_gap.py \
    --allele "HLA-A*02:01" --allele "HLA-A*24:02" --allele "HLA-DRB1*04:01" \
    --out report_dir/
```

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--allele` | (required, repeatable) | IEDB allele name, e.g. `HLA-A*02:01`. Pass multiple times for several alleles in one report. |
| `--match-mode {serotype,exact}` | `serotype` | `serotype` includes the serotype rollup (default). `exact` matches only the allele + higher-resolution sub-alleles. |
| `--out` | `.` | Output directory (created if missing). |
| `--timeout` | `90` | Per-request timeout in seconds. |
| `--table-cap` | `500` | Max candidate rows rendered per allele in the HTML. The **full** candidate set is always written to CSV and JSON. |

### Pipeline (per allele)

1. **Resolve the allele** — one `mhc_search`/`tcell_search` lookup gets its IRI, MRO ancestor IRIs, and MHC class; ancestor IRIs are resolved to names and the serotype is detected.
2. **T-cell side** — page through `tcell_search` for the allele's match set (exact + sub-allele via IRI-contains; serotype via name match), tagging each row's resolution. Aggregate to epitopes with assay counts (positive / negative / other), source organism, parent antigen, diseases, assays.
3. **Binding side** — page through `mhc_search` for the same match set; collect the set of `structure_id`s that have a binding assay.
4. **Gap** — candidate epitopes = T-cell epitopes − binding-tested epitopes. Computed both exact-only and serotype-inclusive.
5. **Render** — `report.html` plus `gap_<allele>.csv` / `gap_<allele>.json`.

Pagination uses PostgREST `Range`/`Content-Range` with `Prefer: count=exact`; responses are capped at 10k rows/request so large alleles (e.g. HLA-A\*02:01 has >200k binding assays) page through in chunks. Transient 5xx/timeouts are retried.

## Output schema

Everything is written inside `--out`:

```
<out>/
├── report.html          # self-contained: overview + per-allele cards/tables + methodology
├── gap_HLA_A_02_01.json # full result for one allele (counts, ontology info, all candidates)
└── gap_HLA_A_02_01.csv  # one row per candidate epitope
```

**`gap_<allele>.csv` columns:**
`structure_id, linear_sequence, length, mhc_class, n_tcell_assays, n_tcell_positive, n_tcell_negative, tcell_resolution, in_exact_gap, source_organism, parent_antigen, host_organism, diseases, assays, example_pubmed, iedb_url`

- `tcell_resolution` — `exact`, `sub-allele`, `serotype`, or a `+`-joined combination, recording at what resolution the epitope's T-cell evidence exists.
- `in_exact_gap` — `True` if the epitope is in the gap even under exact-only matching (strongest candidates).
- `iedb_url` — `https://www.iedb.org/epitope/<structure_id>`.

Candidates are sorted by positive T-cell assays, then total assays. The HTML renders the top `--table-cap` per allele (cap is stated in the report when hit); CSV/JSON hold the full set.

## Workflow guidance for the assistant

1. Always run the script — don't hand-craft IQ-API URLs.
2. Confirm the **exact IEDB allele name** first. `HLA-A2` (serotype), `HLA-A*02:01` (allele), and `HLA-A*02:01:01` (sub-allele) are different nodes. If the user gives a fuzzy name, the script's resolver fails fast with a hint.
3. After the run, point the user at `report.html`. Lead with the **gap count** and the **subset that has a positive T-cell assay** (the highest-interest candidates).
4. Note which gap you're quoting: serotype-inclusive (headline) vs exact-only (`in_exact_gap`). State the serotype that was rolled up so the rollup is transparent.
5. The candidates with `n_tcell_positive > 0` are the ones worth flagging for binding validation; pure-negative-T-cell candidates are weaker.
6. Don't dump the HTML into chat — quote the headline numbers and a couple of standout epitopes.

## Gotchas

- A "gap" is the absence of a *recorded* binding assay, not proof the peptide doesn't bind. Frame it that way.
- "Has a T-cell assay" includes negative assays; always cite the positive subset alongside the total.
- Serotype detection is best-effort (regex on the MRO ancestor names). The chosen serotype is shown per allele; if an exotic allele has no detectable serotype, the run falls back to exact-only matching for that allele and the report shows serotype = none.
- Large class I alleles (HLA-A\*02:01, A\*24:02) have hundreds of thousands of binding assay rows — the binding-side fetch dominates runtime (a few minutes). The script fetches only `structure_id`/allele on that side to keep payloads small.
- Non-peptidic structures are out of scope.

## Extended API patterns

For the full IQ-API endpoint/operator reference and the ontology approach used here, see `references/iq_api_notes.md`.
