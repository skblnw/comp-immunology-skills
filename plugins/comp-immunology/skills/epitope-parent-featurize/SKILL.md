---
name: epitope-parent-featurize
description: Use this skill when the user has a CSV of epitope-to-parent mappings (typically MHC peptide ↔ source protein) and wants ESM C 300M embeddings for BOTH the epitope on its own AND the same epitope as it appears in its parent's context, with a manifest joining them. Triggers on phrases like "embed epitopes and their parents with ESM", "compute direct and contextual epitope features", "featurize epitope-parent pairs", "ESM embeddings for peptides plus parent slices", "vectorise epitopes in context", "build epitope features for MHC/TCR modeling with parent context", "9-mer embeddings sliced from full proteins", "generate per-residue context vectors for epitopes from UniProt sequences". Also use whenever someone supplies a CSV with epitope sequences, parent UniProt IDs, and occurrence ranges plus a directory of parent FASTAs, since that input shape only makes sense for this dual-embedding workflow. This skill builds on top of the `esm-featurize` skill (which it shells out to) — do not re-invoke esm-featurize separately when this skill applies; it handles parent windowing, slice bookkeeping, and the manifest that esm-featurize alone does not provide.
---

# epitope-parent-featurize

## Purpose

Many MHC/TCR / epitope-immunogenicity workflows need two views of every
candidate peptide:

- **direct**: the peptide on its own, fed to ESM C as a standalone
  sequence — useful when downstream models treat peptides as
  context-free strings (as MHC-binding predictors usually do).
- **parent-context**: the same residues, but embedded as part of the full
  source protein (or a window of it). Downstream models that want to
  exploit ESM's context — local secondary structure, conservation, the
  rest of the antigen — need this.

These two embeddings can differ substantially for the same peptide; the
caller often wants both side by side. This skill computes them in one
shot from a CSV plus a directory of parent FASTAs and writes a manifest
that maps each input row to its two `.npz` files plus the slice offsets
into the parent embedding.

## When to use

Use this skill when the user supplies (or describes):

- a **CSV** with at least one row per candidate epitope, containing
  - the **epitope sequence** (peptide string),
  - a **parent UniProt accession** for that epitope, and
  - an **occurrence range** like `251-259` (1-based inclusive) telling
    you where in the parent the epitope sits;
- a **directory of parent FASTAs** named `<ACCESSION>.fasta` (one
  sequence per file);
- the request is "give me embeddings for both the peptide and its
  parent" or "I want context-aware features for these epitopes".

If the user only wants direct embeddings of a flat list of sequences and
has no parent context, use the simpler `esm-featurize` skill instead.

## Prerequisites

This skill calls the `featurize.py` shipped by the `esm-featurize`
skill, which in turn requires a conda env named `esmc` with the `esm`
package installed:

```bash
conda create -n esmc python=3.11 -y
conda activate esmc
pip install esm numpy httpx
```

The driver script in this skill is pure-stdlib + NumPy; it can be
invoked with any Python and the underlying featurize.py auto-reexecs
into the `esmc` env.

## Quick start

```bash
python ~/.claude/plugins/comp-immunology/skills/epitope-parent-featurize/scripts/featurize_pairs.py \
    --csv path/to/epitopes.csv \
    --parents-dir path/to/parents/ \
    --out path/to/output/
```

That's the minimum. The output directory will contain:

```
output/
├── direct/<epitope>.npz                # one per unique epitope sequence
├── parent/<acc>__w<ws>_<we>.npz        # one per unique parent window
├── fasta/{epitopes,parents}.fasta      # inputs to the featurize skill
├── featurize_direct.log,
├── featurize_parent.log                # logs from each model run
├── manifest.csv                        # one row per input row, joins it
│                                       #  all together (read this first)
├── summary.json
└── report.md
```

### Loading a pair from the manifest

```python
import csv, numpy as np

manifest = list(csv.DictReader(open("output/manifest.csv")))
row = manifest[0]
# Direct embedding of the epitope on its own:
direct = np.load(row["direct_npz"])
direct_per_res = direct["per_residue"]          # (L, 960), L = len(epitope)
direct_mean    = direct["mean_pooled"]          # (960,)

# Parent-context embedding sliced at the occurrence positions:
parent = np.load(row["parent_npz"])
ss = int(row["slice_start_0based"])
se = int(row["slice_end_0based"])
ctx_per_res = parent["per_residue"][ss:se]      # (L, 960), same length
ctx_mean    = ctx_per_res.mean(axis=0)          # (960,)
```

The `parent_npz` is the embedding of a **window of the parent** (default
480 residues centred on the occurrence) rather than the whole protein.
This keeps memory and compute bounded while preserving ~235 aa of
flanking context on each side of a 9-mer — far more than any plausible
local-context window. Parents shorter than the window size are embedded
in full.

## CSV columns

The script auto-detects columns by trying common name variants:

| field | candidates (case-sensitive) |
|---|---|
| epitope sequence | `epitope`, `peptide`, `sequence`, `epi`, `epitope_seq` |
| parent UniProt accession | `uniprot_acc`, `uniprot`, `uniprot_id`, `accession`, `parent_uniprot`, `parent_uniprot_id`, `parent_id`, `acc` |
| occurrence range | `occurrences`, `occurrence`, `position`, `positions`, `range`, `epitope_position` |

Override with `--col-epitope`, `--col-uniprot`, `--col-occurrences` if
the user's CSV uses different names. Extra columns are ignored.

The occurrence value must be a single 1-based inclusive range like
`251-259`. Multi-range values (`"164-172, 170-178"`) are rejected —
they're rare enough that special-casing them is not worth it.

## Row filtering

Rows are dropped (with reason counts in `summary.json`) when:

- the epitope contains non-canonical AAs (`X`, `U`, `B`, `Z`, `*`, `-`),
- the parent FASTA isn't in `--parents-dir` (`<acc>.fasta` or `<acc>.fa`),
- the parent sequence contains non-canonical AAs,
- the occurrence doesn't parse, or its length doesn't match the epitope,
- the occurrence is out of bounds for the parent,
- `parent[start-1:end] != epitope` — i.e. the position annotation
  disagrees with the sequence. This catches off-by-one errors in the
  input CSV and prevents silently comparing unrelated tensors. Pass
  `--position-check warn` to keep such rows with a logged warning
  instead of dropping them.

The skill always verifies the slice equals the epitope before writing
the manifest — this is the single most useful sanity check for an
epitope-context workflow and trying to skip it has burned past callers.

## Common flags

| flag | default | meaning |
|---|---|---|
| `--window N` | `480` | window length for parents longer than N. Smaller = more local context only; larger = more flanking residues but ESM C's stated max is 2048. 480 gives ~235 aa of flank, which covers any reasonable local-context effect for a 9-mer. |
| `--n K` | `0` (all) | sample K rows after filtering (for smoke tests). Reproducible via `--seed`. |
| `--sequential` | off | run direct + parent embedding sequentially (default: in parallel). Sequential is useful when MPS memory is tight. |
| `--skip-featurize` | off | reuse existing npz outputs; only rebuild the manifest + report. Use this after changing the manifest format or fixing a CSV bug — don't waste model time. |
| `--position-check {drop,warn}` | `drop` | what to do when `parent[start-1:end] != epitope`. |
| `--featurize-script` | auto | override the path to `esm-featurize`'s featurize.py. Auto-discovery checks, in order: (1) `scripts/featurize.py` alongside this skill, (2) `~/.claude/skills/esm-featurize/scripts/featurize.py` (legacy standalone), (3) `~/.claude/plugins/*/skills/esm-featurize/scripts/featurize.py` (any installed plugin — currently `structural-bioinfo`). |

## Output sizes

For a CSV of N input rows that produces E unique epitope sequences and W
unique parent windows after dedup:

- `direct/` contains E files; each is small (length × 960 floats).
- `parent/` contains W files; each is up to `window × 960` floats
  (~1.9 MB for a 480-residue window).
- The model runs once per unique epitope and once per unique window, not
  once per CSV row. Duplicate epitopes or windows are embedded once and
  referenced N times from `manifest.csv`.

## Verification

After a run, the caller should:

1. Open `manifest.csv` and pick any row.
2. Load `parent_npz`, slice `per_residue[slice_start:slice_end]`, and
   confirm `parent["sequence"][slice_start:slice_end] == epitope`.
3. Load `direct_npz` and confirm `per_residue.shape == (len(epitope), 960)`.
4. Check `summary.json` for the drop counts — high `position_mismatch`
   counts usually mean the parent FASTAs were downloaded from a
   different UniProt release than the CSV was generated against.

## Notes on design choices

- **Why a window, not the full protein?** ESM C's stated max is 2048
  residues; embedding an 8,797-aa polyprotein wastes most of the
  signal for a 9-mer and risks model edge-effects. A 480-residue window
  centred on the occurrence keeps every "context" embedding the same
  length so the downstream comparison is apples-to-apples.
- **Why dedupe by `(parent_acc, window_start, window_end)`?** A single
  parent often hosts dozens of epitopes (max observed: 217 from a
  single Ebola Nucleoprotein). Without dedup we'd embed the same
  window many times. With dedup, each unique window is embedded once
  and the manifest joins multiple rows to it.
- **Why verify the slice?** Position annotations in epitope databases
  drift relative to UniProt sequence versions. Silent off-by-ones
  produce embeddings that *look* fine but are sliced at the wrong
  residues — a uniquely insidious failure mode.
