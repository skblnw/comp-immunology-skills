# comp-immunology

Computational immunology toolkit for Claude Code — featurize MHC / TCR / immunogenicity inputs with protein language model embeddings that respect parent-protein context.

## Skills

### `epitope-parent-featurize`

Given a CSV of epitope-to-parent mappings (MHC peptide ↔ source UniProt protein with occurrence range) plus a directory of parent FASTAs, compute **two** ESM C 300M embeddings per row: the epitope on its own (`direct`), and the same residues as they appear inside a window of their parent (`parent-context`). Outputs a manifest joining each input row to its `.npz` pair with slice offsets, so downstream MHC / TCR / immunogenicity models can train on either or both views.

**Triggers:** "Embed epitopes and their parents with ESM", "Featurize epitope-parent pairs", "Compute direct and contextual epitope features", "Build epitope features for MHC/TCR modeling with parent context", "9-mer embeddings sliced from full proteins"

**Output:** `direct/<epitope>.npz` (one per unique peptide), `parent/<acc>__w<ws>_<we>.npz` (one per unique parent window), `manifest.csv`, `summary.json`, `report.md`. Slice verification (`parent[start-1:end] == epitope`) is enforced by default — silent off-by-one annotations are dropped, not embedded.

### `iedb-assay-gap`

For a given MHC allele, find epitopes in the IEDB that **have T-cell assays but no MHC-binding assay** against that allele — the "binding-data gap": peptides proven T-cell antigenic in the allele's context whose direct MHC binding has never been experimentally measured. Allele matching is grounded in the **MHC Restriction Ontology (MRO)** that the IQ-API exposes via the `mhc_allele_iri_search` array, so a serotype rollup (e.g. `HLA-A2` for `HLA-A*02:01`, `HLA-DR4` for `HLA-DRB1*04:01`) is included and every match is tagged by resolution (`exact` / `sub-allele` / `serotype`).

**Triggers:** "Which epitopes have a T-cell assay for HLA-A\*02:01 but no MHC-binding assay?", "Find binding-validation candidates for an allele", "Quantify the gap between T-cell and MHC-binding coverage", "Peptides that are T-cell positive on HLA-A\*24:02 but never tested for binding"

**Output:** a self-contained `report.html` (cross-allele overview, per-allele cards + sortable candidate tables, methodology/caveats) plus `gap_<allele>.csv` / `.json` (one row per candidate epitope, with positive/negative T-cell breakdown, restriction resolution, source organism, disease, and IEDB link). Stdlib-only; queries the IEDB IQ-API (`tcell_search` + `mhc_search`). Pages with keyset (cursor) pagination so large alleles (HLA-A\*02:01 has >200k binding assays) stay fast.

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
