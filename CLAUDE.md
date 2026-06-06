# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **Claude Code plugin marketplace** shipping a single plugin, `comp-immunology`, with two computational-immunology skills. There is no application to run and no build step — the deliverable is the skill definitions (`SKILL.md`) and their bundled Python scripts. End users install it with `/plugin marketplace add skblnw/comp-immunology-skills` then `/plugin install comp-immunology`.

This repo is a **downstream mirror** of a local plugin developed at `~/.claude/plugins/comp-immunology/`. Sync is manual: rsync the plugin into `plugins/comp-immunology/`, commit, push. Do not assume changes here propagate back upstream.

## Layout

```
.claude-plugin/marketplace.json          # marketplace manifest → points at plugins/comp-immunology
plugins/comp-immunology/
├── .claude-plugin/plugin.json           # plugin manifest (name, version, keywords)
├── README.md                            # source of truth for the auto-generated SKILLS block
└── skills/<skill>/
    ├── SKILL.md                         # YAML frontmatter (name + description trigger phrases) + body
    ├── scripts/*.py                     # the actual implementation
    ├── references/*.md                  # extended docs the skill body links to (optional)
    └── evals/evals.json                 # smoke-test cases (optional)
```

The `## Skills` block in the root `README.md` is auto-generated from the plugin README (see the `BEGIN SKILLS`/`END SKILLS` markers) — edit skill descriptions in the plugin README, not by hand in the root README.

## The two skills

- **`epitope-parent-featurize`** — `featurize_pairs.py`. Takes a CSV of epitope→parent mappings plus a dir of parent FASTAs; emits two ESM C 300M embeddings per row (the peptide alone = `direct`, and the same residues inside a window of their parent = `parent-context`), joined by `manifest.csv` with 0-based slice offsets. It does **not** run the model itself — it shells out to `featurize.py` from the `esm-featurize` skill (see cross-plugin dependency below). Pure stdlib + NumPy.
- **`iedb-assay-gap`** — `iedb_assay_gap.py`. For an MHC allele, finds epitopes with a T-cell assay but no MHC-binding assay (the "binding-data gap") by querying the IEDB IQ-API. **Stdlib only, zero dependencies.**

## Cross-plugin dependency (epitope-parent-featurize only)

`featurize_pairs.py` never embeds anything itself — it calls the `featurize.py` from the **`esm-featurize`** skill, which lives in the separate [`structural-bioinfo`](https://github.com/skblnw/structural-bioinfo-skills) plugin. `discover_featurize_script()` auto-locates it in this order: (1) bundled next to the script, (2) `~/.claude/skills/esm-featurize/scripts/featurize.py` (legacy), (3) `~/.claude/plugins/*/skills/esm-featurize/scripts/featurize.py`. `esm-featurize` requires a conda env named `esmc` (`conda create -n esmc python=3.11 -y && pip install esm numpy httpx`); `featurize.py` re-execs itself into that env. The driver here runs under any Python.

## Running / testing the scripts

No test runner, linter, or CI is configured — validation is by running the scripts directly.

```bash
# epitope-parent-featurize (needs the esmc env + esm-featurize installed)
python plugins/comp-immunology/skills/epitope-parent-featurize/scripts/featurize_pairs.py \
    --csv epitopes.csv --parents-dir parents/ --out out/
# smoke test without the full model: sample N rows, or reuse npz + rebuild manifest
#   --n 5            sample 5 rows after filtering (reproducible via --seed)
#   --skip-featurize reuse existing npz, only rebuild manifest+report

# iedb-assay-gap (stdlib only; hits the live IQ-API, runs a few minutes for big alleles)
python plugins/comp-immunology/skills/iedb-assay-gap/scripts/iedb_assay_gap.py \
    --allele "HLA-A*02:01" -o out/
```

`iedb-assay-gap` ships `evals/evals.json` — declarative cases (expected allele IRI, serotype rollup, gap invariants) validated against the live API. The file's `how_to_run` field has the exact command.

## Conventions that matter here

- **Stdlib-first.** Scripts depend only on the Python standard library (plus NumPy for `featurize_pairs.py`). Don't add third-party deps — it's a hard design constraint that keeps the skills installable without an environment.
- **Self-contained outputs.** Every run writes a machine-readable artifact (`summary.json` / `gap_*.json`), a CSV, and a **self-contained `report.html`** with inlined CSS and no JS or external assets. Match this pattern when adding output.
- **Correctness checks are load-bearing, not optional.** `featurize_pairs.py` verifies `parent[start-1:end] == epitope` before writing the manifest and drops mismatches by default (silent off-by-one slices are the insidious failure mode this guards against). `iedb_assay_gap.py` grounds allele matching in the MRO ontology (`mhc_allele_iri_search` ancestor IRIs) so exact / sub-allele / serotype resolutions are disjoint and auditable. Preserve these guarantees; don't weaken them for convenience.
- **Dedup before embedding.** `featurize_pairs.py` embeds each unique epitope and each unique `(parent_acc, window_start, window_end)` once, then references them from many manifest rows — a single parent can host hundreds of epitopes.
- **A "gap" means no *recorded* assay, not proof of non-binding** — frame `iedb-assay-gap` results that way, and always cite the positive-T-cell subset alongside the total.
- **SKILL.md descriptions are the trigger surface.** The `description` frontmatter is dense with natural-language trigger phrases on purpose — that's how the skill gets matched to user intent. Keep it that way when editing.

## License

MIT.
