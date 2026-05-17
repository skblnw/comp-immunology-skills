# Comp Immunology Skills

Claude Code skills for computational immunology — featurize MHC / TCR / immunogenicity inputs with protein language model embeddings that respect parent-protein context. Stdlib Python + NumPy on top of ESM C — no API keys, no Forge.

## Quickstart

In Claude Code:

```
/plugin marketplace add skblnw/comp-immunology-skills
/plugin install comp-immunology
```

This plugin's skills shell out to `esm-featurize`, which lives in the `structural-bioinfo` plugin. Install that one alongside:

```
/plugin marketplace add skblnw/structural-bioinfo-skills
/plugin install structural-bioinfo
```

Then just ask Claude in natural language — the right skill triggers on intent.

## Why these exist

Epitope-immunogenicity, MHC-binding, and TCR-specificity models all live or die on how the peptide is represented. The two embeddings most teams end up wanting are:

- the peptide on its own (context-free, how MHC-binding predictors see it), and
- the same residues embedded as part of the full parent protein (context-aware, what ESM C actually adds).

Computing both, deduping unique peptides and parent windows, and joining them back to the input CSV with slice offsets is fiddly enough to keep getting redone in notebooks. These skills are that workflow, cleaned up and exposed as Claude Code triggers.

## Skills

<!-- BEGIN SKILLS -- auto-generated from the plugin README; do not edit by hand -->
- **[epitope-parent-featurize](plugins/comp-immunology/skills/epitope-parent-featurize/SKILL.md)** — Compute two ESM C 300M embeddings per epitope row: the peptide on its own (`direct`) and the same residues inside a window of its parent protein (`parent-context`). Outputs deduped `.npz` files plus a manifest joining each input row to both views with slice offsets. Slice verification (`parent[start-1:end] == epitope`) is enforced to catch silent off-by-one annotations.
<!-- END SKILLS -->

## Typical pipeline

```bash
python plugins/comp-immunology/skills/epitope-parent-featurize/scripts/featurize_pairs.py \
    --csv path/to/epitopes.csv \
    --parents-dir path/to/parents/ \
    --out path/to/output/
```

Read `output/manifest.csv` first — every input row points to its `direct/<epitope>.npz` and `parent/<acc>__w<ws>_<we>.npz` plus the slice offsets into the parent embedding.

## Requirements

- Python 3.8+ with NumPy (the manifest driver is stdlib + NumPy)
- The `structural-bioinfo` plugin (provides `esm-featurize`, which actually runs ESM C)
- A conda env named `esmc` with the `esm` package installed:

  ```bash
  conda create -n esmc python=3.11 -y
  conda activate esmc
  pip install esm numpy httpx
  ```

## Updating

This repo is a downstream mirror of the local plugin at `~/.claude/plugins/comp-immunology/`. Sync is manual today — rsync the plugin into `plugins/comp-immunology/`, then commit and push.

## License

MIT — see [LICENSE](LICENSE).
