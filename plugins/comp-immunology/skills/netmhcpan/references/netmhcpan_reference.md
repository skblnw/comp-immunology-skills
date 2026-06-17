# NetMHCpan reference

Detailed reference for the local NetMHCpan installs. Read the section you need.

- [Output columns](#output-columns)
- [Score interpretation](#score-interpretation)
- [Allele naming & listing](#allele-naming--listing)
- [Useful CLI flags](#useful-cli-flags)
- [Version differences (4.1b vs 4.2c)](#version-differences)
- [What the model was trained on](#training-data-background)

## Output columns

NetMHCpan emits two layouts. The space-aligned **`.out`** (stdout) and the tab-separated
**`.xls`** (`-xls -xlsfile`). Prefer `.xls` for parsing.

**`.out` columns (per row, allele already on the row):**

| Column | Meaning |
|---|---|
| `Pos` | start position of the peptide in the input (1-based in 4.2, 0-based in 4.1) |
| `MHC` | allele, printed as `HLA-A*02:01` (note the `*`; the parser strips it) |
| `Peptide` | the input peptide |
| `Core` | predicted 9mer binding core |
| `Of Gp Gl Ip Il` | alignment bookkeeping: offset, gap position/length, insertion position/length (usually all 0 for plain 9mers) |
| `Icore` | the residues actually presented (interaction core) — differs from Peptide only for bulged/longer ligands |
| `Identity` | the FASTA/sequence name, or `PEPLIST` for `-p` input |
| `Score_EL` | eluted-ligand score (presentation likelihood, 0–1) |
| `%Rank_EL` | EL score ranked vs random peptides for this allele — **use this for calls** |
| `Score_BA` `%Rank_BA` `Aff(nM)` | binding-affinity score, its %Rank, and predicted IC50 in nM (only with `-BA`) |
| `BindLevel` | `<= SB` (strong) or `<= WB` (weak), blank otherwise |

**`.xls` columns:** same quantities, but **wide** — `Pos Peptide ID`, then a repeated block
of `core icore EL_score EL_rank [BA_score BA_rank]` **per allele**, then trailing `Ave NB`
(mean score and number of binding alleles across the requested alleles, per peptide). The
allele names sit in the row above the header (the row pandas skips). 4.1 uses hyphenated
labels (`EL-score`, `EL_Rank`); 4.2 uses underscores (`EL_score`, `EL_rank`).

## Score interpretation

- **Always rank on %Rank, not raw score.** Raw EL/BA scores are not comparable between
  alleles; %Rank is normalized per allele against a large random-peptide pool, so the same
  threshold means the same thing everywhere.
- **EL %Rank default thresholds:** SB (strong binder) ≤ 0.5%, WB (weak binder) ≤ 2.0%.
  Change the run with `-rth`/`-rlt`; if you do, pass matching `--sb/--wb` to the parser.
- **Affinity (`Aff(nM)`):** lower = tighter. Conventional cutoffs ≤ 50 nM strong, ≤ 500 nM weak.
- **EL vs BA:** EL models actual presentation (MS eluted ligands; the whole processing
  pathway) and is the better presentation predictor. BA is pure binding strength — useful when
  affinity per se matters. They can disagree; report whichever the question asks for, and say which.
- **core vs icore:** for ordinary 9mers both equal the peptide. For 8/10/11+mers, `icore` is
  what's presented; `core` is the aligned 9mer used internally.

## Allele naming & listing

- Format: `HLA-A02:01` (locus letter, 2-digit allele group, colon, 2-digit protein). No `*` on input.
- Multiple alleles: comma-separated, no spaces — `-a HLA-A02:01,HLA-B07:02`.
- Non-human supported too: mouse `H-2-Kb`, bovine `BoLA-*`, swine `SLA-*`, macaque `Mamu-*`, etc.
- List everything an install supports: `netMHCpan -listMHC` (12k+ named alleles; the model is
  pan-specific so any allele with a known pseudosequence can be predicted).
- User-defined MHC by sequence: `-hlaseq mhc.fasta` (skip the named-allele list entirely).

## Useful CLI flags

| Flag | Effect |
|---|---|
| `-p FILE` | input is a peptide list (one per line) |
| `-f FILE` | input is FASTA (default mode); window length set by `-l` |
| `-l 8,9,10,11` | peptide lengths to scan from FASTA (default 9) |
| `-a A,B,C` | alleles (comma-separated) |
| `-BA` | add binding-affinity prediction |
| `-xls -xlsfile OUT` | write the tab-separated table (parse this) |
| `-hlaseq FILE` | predict for a user MHC sequence |
| `-rth -rlt` | strong / weak %Rank thresholds (defaults 0.5 / 2.0) |
| `-listMHC` | print supported alleles |
| `-s` | sort output by score |
| `-xlsfile` path must be writable; 4.1 also needs `gawk` on PATH or this file is empty |

## Version differences

| | 4.2c (`netMHCpan`) | 4.1b (`netMHCpan-4.1`) |
|---|---|---|
| Binary | native arm64 | Intel x86_64 (Rosetta 2) |
| `.xls` line 1 | `#<command>` → `skiprows=2` | allele-group row → `skiprows=1` |
| `.xls` score labels | `EL_score`, `EL_rank`, `BA_score`, `BA_rank` | `EL-score`, `EL_Rank`, `BA-score`, `BA_Rank` |
| `Pos` base | 1 | 0 |
| `-xls` needs gawk | no | **yes** |
| Training data | adds CEDAR cancer epitopes + updated IEDB | NAR-2020 release |

EL scores are correlated but not identical across versions (the EL net was retrained), so
results from the two are not interchangeable — never pool them in one benchmark.

## Training data background

NetMHCpan-4.x is pan-specific: each MHC is reduced to a 34-residue pseudosequence (the
groove-lining polymorphic positions), so it generalizes to any allele with a known sequence.
It co-trains on binding affinity (BA) and mass-spec eluted-ligand (EL) data; EL negatives are
random decoy peptides. The NetMHCpan-4.2 training/eval corpus (≈17.3M EL records, ≈215K BA
records, 5-fold CV, IEDB + CEDAR benchmarks) is unpacked locally under
`immunology/netmhc_data/`, with a full teardown in `immunology/netmhc_data/netmhc_dataset_report.html`.
Consult that when the goal is to compare a new method against NetMHCpan or to understand its
coverage and weaknesses (allele imbalance, decoy label noise, 9mer bias, presentation-only).
