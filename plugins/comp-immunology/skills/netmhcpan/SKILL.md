---
name: netmhcpan
description: Run the two locally-installed NetMHCpan MHC class I predictors correctly and parse their output — reach for this for ANY peptide–MHC-I binding or antigen-presentation prediction, even when NetMHCpan isn't named. NetMHCpan-4.2c (native arm64) is the command `netMHCpan`; NetMHCpan-4.1b (Intel via Rosetta) is `netMHCpan-4.1`. Use whenever the user wants to predict which peptides bind or are presented by an HLA / MHC class I allele (human HLA-A/B/C/E/G or non-human, e.g. mouse H-2, BoLA), get eluted-ligand (EL) or binding-affinity (BA) scores, %Rank or strong/weak (SB/WB) calls, pick or rank neoantigen / epitope candidates, or scan a protein FASTA for binders — e.g. "run netmhcpan on these peptides", "which of these peptides are presented by HLA-A*02:01", "predict MHC class I binders for the patient's HLA type A*11:01/B*15:01", "what's the EL %rank of GILGFVFTL", "rank these neoantigens by presentation", "scan this spike FASTA for HLA-A24:02 binders". Also use when choosing between or comparing the two NetMHCpan versions, when a run errors or writes an empty / 0-byte .xls, or when reading / melting a NetMHCpan .xls or .out into a tidy table. Encodes the version pitfalls (which command is which, Rosetta, the gawk dependency, .xls skiprows 1-vs-2, never mixing 4.1 and 4.2 scores) and bundles a version-aware output parser plus an install doctor. Predicts presentation/binding only — NOT class II (use NetMHCIIpan), NOT T-cell immunogenicity / TCR-recognition scoring, NOT pulling epitopes or assays from the IEDB database (use iedb-query), and NOT pMHC structure lookup.
---

# NetMHCpan (local 4.2c + 4.1b)

Two NetMHCpan installs live on this Mac. They predict, for an MHC class I allele and a
peptide, how likely the peptide is **presented** (eluted-ligand, EL) and optionally how
tightly it **binds** (binding affinity, BA). This skill is about driving them correctly and
turning their output into clean data — the traps are all in the version differences, not in
the biology.

| | Command | Version | Runtime | Use it for |
|---|---|---|---|---|
| **Default** | `netMHCpan` | 4.2c | native arm64 | New predictions — newest EL network, includes CEDAR cancer-epitope training |
| **Legacy** | `netMHCpan-4.1` | 4.1b | Intel x86_64 via **Rosetta 2** | Reproducing older 4.1-era results, or 4.1-vs-4.2 comparisons |

Both are exposed on `PATH` via `~/.local/bin`. Default to **4.2c (`netMHCpan`)** unless the
user is explicitly reproducing or comparing against 4.1.

## Verify the installs first

If anything seems off (a fresh machine, a run that errors, an empty `.xls`), run the doctor —
it checks both commands, prerequisites, and a live prediction, and pinpoints what's wrong:

```bash
bash scripts/check_install.sh
```

## Quick start

```bash
# A peptide list (one peptide per line), single allele, EL only (4.2c)
netMHCpan -p peptides.pep -a HLA-A02:01

# EL + binding affinity, several alleles, structured .xls output for parsing
netMHCpan -p peptides.pep -BA -a HLA-A01:01,HLA-A02:01,HLA-B07:02 -xls -xlsfile out.xls

# Scan a protein FASTA for 9mers (default length) — predict & rank all sub-peptides
netMHCpan -f antigen.fasta -a HLA-A02:01 -l 9

# Legacy 4.1b — identical flags, just the other command
netMHCpan-4.1 -p peptides.pep -BA -a HLA-A02:01 -xls -xlsfile out_v41.xls
```

Then parse (handles either version, either format, and melts multi-allele tables):

```bash
python scripts/parse_netmhcpan.py out.xls --out tidy.csv
```

## Input modes

- `-p FILE` — peptide list, one sequence per line (the common case). Add `-inptype 1` only if NetMHCpan misreads it.
- `-f FILE` — FASTA; NetMHCpan slides a window of length `-l` over each protein. `-l 8,9,10,11` for multiple lengths (9 is default and dominant for class I).
- `-a ALLELES` — comma-separated, **`HLA-A02:01` format** (locus letter, no `*`, colon before the 2-digit field). `H-2-Kb`, `BoLA-...`, etc. for non-human. List supported alleles with `-listMHC`.
- `-BA` — also predict binding affinity (adds BA columns + `Aff(nM)`). Without it you get EL only.
- `-hlaseq FILE.fasta` — predict for a user-supplied MHC sequence instead of a named allele.
- `-context` — context-aware EL (requires context-formatted input); rarely needed for ad-hoc scoring.
- `-xls -xlsfile OUT.xls` — write the parseable tab-separated table. **Always prefer `-xls` for downstream parsing** over capturing stdout.

## Reading the output

Score every peptide on **%Rank**, not the raw score — raw EL/BA scores are not comparable
across alleles, but %Rank is normalized per allele (against a large pool of random peptides),
so a 0.5% rank means the same thing for HLA-A*02:01 and HLA-B*07:02.

- **EL_score / %Rank_EL** — presentation likelihood. This is the primary readout. Default bind levels: **SB (strong) ≤ 0.5%**, **WB (weak) ≤ 2.0%**.
- **BA_score / %Rank_BA / Aff(nM)** — binding affinity (only with `-BA`). `Aff(nM)`: lower = tighter; ~≤50 nM strong, ≤500 nM weak, by convention.
- **icore** — the peptide actually presented (after any predicted bulge/deletion); **core** — the predicted 9mer binding core. For ordinary 9mers both equal the input.
- A peptide with low %Rank_EL is "likely presented"; combine with BA if you care about affinity specifically.

Full column glossary, allele-naming rules, and the training-data background are in
`references/netmhcpan_reference.md` — read it when a column is unfamiliar or when the user
asks what the model was trained on / how to compare against it.

## Critical gotchas (this is why the skill exists)

These are version-specific and cost real time when hit blind:

1. **Two different commands.** `netMHCpan` = 4.2c, `netMHCpan-4.1` = 4.1b. Picking the wrong one silently gives different numbers.

2. **Never mix 4.1 and 4.2 scores in one analysis.** The EL network was retrained between versions, so EL_score / %Rank_EL shift (≈0.98 correlation but enough boundary call-flips to matter; BA tracks much tighter). Pick one version for a given experiment and stay on it.

3. **4.1 `-xls` needs `gawk`.** The 4.1b binary shells out to `gawk` to write the `.xls`; without it the file is **0 bytes while stdout still looks fine** — a silent failure. `brew install gawk` fixes it. (4.2c does not need gawk.)

4. **`.xls` parsing differs by version** — the #1 reason a parse comes out garbled:
   - **4.2c**: line 1 is a `#<command>` line → pandas `skiprows=2`. Headers use underscores (`EL_score`, `EL_rank`, `BA_score`, `BA_rank`). Position column starts at **1**.
   - **4.1b**: no command line → pandas `skiprows=1`. Headers use hyphens/mixed case (`EL-score`, `EL_Rank`, `BA-score`, `BA_Rank`). Position column starts at **0**.
   Don't hand-roll this — `scripts/parse_netmhcpan.py` auto-detects the version and normalizes columns.

5. **Multi-allele `.xls` is "wide".** Each requested allele gets its own block of columns, with the allele name in the skipped allele-group row above its block. The parser melts this into tidy one-row-per-(peptide, allele) records; doing it by hand is where mistakes happen.

6. **4.1 runs under Rosetta** (Intel binary on Apple Silicon). It works transparently, but Rosetta 2 must be installed (`softwareupdate --install-rosetta`); the doctor checks this.

## Parsing — use the bundled script

`scripts/parse_netmhcpan.py` is both a CLI and an importable library. It auto-detects 4.1
vs 4.2 and `.xls` vs `.out`, melts multi-allele tables, normalizes the column names across
versions, and derives the SB/WB bind level from %Rank_EL when it isn't in the file.

```bash
# CLI — tidy CSV (default) or JSON; accepts multiple files
python scripts/parse_netmhcpan.py out.xls
python scripts/parse_netmhcpan.py run_v41.xls run_v42.xls --format json --out all.json
python scripts/parse_netmhcpan.py run.xls --sb 0.5 --wb 2.0      # custom thresholds if you used -rth/-rlt
```

```python
from parse_netmhcpan import parse_file          # add scripts/ to sys.path
records = parse_file("out.xls")                  # list[dict]: pos, peptide, allele, el_score, el_rank, ba_*, bind_level, version, ...
# df = to_dataframe(records)                     # optional, needs pandas
```

Tidy columns: `pos, peptide, id, allele, core, icore, el_score, el_rank, ba_score, ba_rank, aff_nM, bind_level, version, source_file`.

## Scope

MHC **class I** presentation/binding only (the local install). Class II (HLA-DR/DQ/DP) needs
NetMHCIIpan, which is not covered here. This skill predicts presentation/binding — not T-cell
immunogenicity.
