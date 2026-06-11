# comp-immunology

Computational immunology toolkit for Claude Code — featurize MHC / TCR / immunogenicity inputs with protein language model embeddings that respect parent-protein context.

## Skills

### `epitope-parent-featurize`

Given a CSV of epitope-to-parent mappings (MHC peptide ↔ source UniProt protein with occurrence range) plus a directory of parent FASTAs, compute **two** ESM C 300M embeddings per row: the epitope on its own (`direct`), and the same residues as they appear inside a window of their parent (`parent-context`). Outputs a manifest joining each input row to its `.npz` pair with slice offsets, so downstream MHC / TCR / immunogenicity models can train on either or both views.

**Triggers:** "Embed epitopes and their parents with ESM", "Featurize epitope-parent pairs", "Compute direct and contextual epitope features", "Build epitope features for MHC/TCR modeling with parent context", "9-mer embeddings sliced from full proteins"

**Output:** `direct/<epitope>.npz` (one per unique peptide), `parent/<acc>__w<ws>_<we>.npz` (one per unique parent window), `manifest.csv`, `summary.json`, `report.md`. Slice verification (`parent[start-1:end] == epitope`) is enforced by default — silent off-by-one annotations are dropped, not embedded.

### `epitopes-by-source`

Pull **every epitope in the IEDB from a given source organism** (a pathogen / taxon) and build a self-contained HTML report. Source resolution is forgiving — a common abbreviation (`EBV`, `CMV`, `HIV-1`, `SARS-CoV-2`, `dengue`, `M. tuberculosis`, …), a bare NCBITaxon id (`NCBITaxon:10376` / `10376`), or a substring of the IEDB organism name — and the query rolls up the **NCBITaxon subtree** (`source_organism_iri_search=cs.{IRI}`) so all strains/descendants are included. The general query layer is `iedb-query`; this is the specific "all epitopes for pathogen X" report.

**Triggers:** "Fetch the epitopes for EBV", "All IEDB epitopes from SARS-CoV-2", "How many dengue epitopes are in IEDB and from which proteins?", "Report of every epitope whose source organism is *M. tuberculosis*", "IEDB epitopes by source organism"

**Output:** a self-contained `report.html` (summary cards, pure-CSS breakdown bars for top source antigens / MHC class / peptide length / host / allele / strain, and a sortable epitope table) plus `epitopes_<IRI>.csv` / `.json` (full set, one row per epitope). Breakdowns cover linear vs discontinuous, T-cell / B-cell / MHC-ligand assay coverage, positive-assay count, MHC class distribution, and per-strain counts. Optional filters: `--host`, `--positive-only`, `--mhc-class`, `--linear-only`, length range. Stdlib-only; queries the IEDB IQ-API (`epitope_search`) with keyset pagination.

### `iedb-assay-gap`

For a given MHC allele, find epitopes in the IEDB that **have T-cell assays but no MHC-binding assay** against that allele — the "binding-data gap": peptides proven T-cell antigenic in the allele's context whose direct MHC binding has never been experimentally measured. Allele matching is grounded in the **MHC Restriction Ontology (MRO)** that the IQ-API exposes via the `mhc_allele_iri_search` array, so a serotype rollup (e.g. `HLA-A2` for `HLA-A*02:01`, `HLA-DR4` for `HLA-DRB1*04:01`) is included and every match is tagged by resolution (`exact` / `sub-allele` / `serotype`).

**Triggers:** "Which epitopes have a T-cell assay for HLA-A\*02:01 but no MHC-binding assay?", "Find binding-validation candidates for an allele", "Quantify the gap between T-cell and MHC-binding coverage", "Peptides that are T-cell positive on HLA-A\*24:02 but never tested for binding"

**Output:** a self-contained `report.html` (cross-allele overview, per-allele cards + sortable candidate tables, methodology/caveats) plus `gap_<allele>.csv` / `.json` (one row per candidate epitope, with positive/negative T-cell breakdown, restriction resolution, source organism, disease, and IEDB link). Stdlib-only; queries the IEDB IQ-API (`tcell_search` + `mhc_search`). Pages with keyset (cursor) pagination so large alleles (HLA-A\*02:01 has >200k binding assays) stay fast.

### `iedb-query`

A general-purpose, efficient client **and** CLI for the IEDB Query API (IQ-API) — the public PostgREST interface to the Immune Epitope Database. Use it to count, filter, paginate, and export records from any IEDB endpoint (T-cell / MHC / B-cell assays, TCR/BCR receptors, epitopes, antigens, references) and to roll a query up an ontology. Its core trick: every search endpoint carries `<entity>_iri_search` arrays of an annotation's ancestor IRIs, so one `cs.{IRI}` filter matches a node **and all descendants** — generalized across MHC allele (MRO, incl. serotype rollup), organism (NCBITaxon), disease (DOID), antigen, and assay. Built-in keyset (cursor) pagination (the API caps responses at 10k rows), cheap `count=estimated` size gauging, retry/backoff, and CSV/JSON streaming. Importable as `import iedb` or runnable as `python iedb.py endpoints|schema|count|query|resolve`. Stdlib-only, read-only, no auth.

**Triggers:** "Pull/count/export IEDB records for an allele/organism/disease/antigen", "Query the IEDB efficiently", "Paginate a large IEDB result set", "How many IEDB T-cell assays are there for human host?", "Resolve an allele/organism to its ontology IRI and roll the query up to the serotype/genus"

**Output:** streamed CSV or JSON of the matched records (to file or stdout), an estimated row count printed before any pull, and for `resolve` an IRI + ancestor chain (+ HLA serotype for alleles) with ready-to-paste subtree-filter hints. Companion to `iedb-assay-gap`, which is the specific "T-cell-but-no-binding gap" analysis; `iedb-query` is the general query layer.

## Cross-plugin dependency

This skill shells out to `esm-featurize` (which actually runs ESM C). `esm-featurize` lives in the [`structural-bioinfo`](https://github.com/skblnw/structural-bioinfo-skills) plugin, not here, so install both:

```
/plugin marketplace add skblnw/structural-bioinfo-skills
/plugin install structural-bioinfo

/plugin marketplace add skblnw/comp-immunology-skills
/plugin install comp-immunology
```

The driver auto-discovers `esm-featurize` from `~/.claude/plugins/*/skills/esm-featurize/scripts/featurize.py` — no configuration needed once both plugins are installed.

## Installation

### Claude Code

```bash
cc --plugin-dir ~/.claude/plugins/comp-immunology
```

### OpenCode (via OpenPackage)

```bash
command -v opkg >/dev/null || npm install -g opkg
opkg install ~/.claude/plugins/comp-immunology --platforms opencode
```

## Requirements

- Python 3.8+ with NumPy (the manifest driver is stdlib + NumPy)
- A conda env named `esmc` with the `esm` package installed (required by `esm-featurize`):
  ```bash
  conda create -n esmc python=3.11 -y
  conda activate esmc
  pip install esm numpy httpx
  ```
- The `structural-bioinfo` plugin (provides `esm-featurize`)

## License

MIT
