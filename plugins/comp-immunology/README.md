# comp-immunology

Computational immunology toolkit for Claude Code — featurize MHC / TCR / immunogenicity inputs with protein language model embeddings that respect parent-protein context.

## Skills

### `epitope-parent-featurize`

Given a CSV of epitope-to-parent mappings (MHC peptide ↔ source UniProt protein with occurrence range) plus a directory of parent FASTAs, compute **two** ESM C 300M embeddings per row: the epitope on its own (`direct`), and the same residues as they appear inside a window of their parent (`parent-context`). Outputs a manifest joining each input row to its `.npz` pair with slice offsets, so downstream MHC / TCR / immunogenicity models can train on either or both views.

**Triggers:** "Embed epitopes and their parents with ESM", "Featurize epitope-parent pairs", "Compute direct and contextual epitope features", "Build epitope features for MHC/TCR modeling with parent context", "9-mer embeddings sliced from full proteins"

**Output:** `direct/<epitope>.npz` (one per unique peptide), `parent/<acc>__w<ws>_<we>.npz` (one per unique parent window), `manifest.csv`, `summary.json`, `report.md`. Slice verification (`parent[start-1:end] == epitope`) is enforced by default — silent off-by-one annotations are dropped, not embedded.

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
