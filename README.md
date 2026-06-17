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
- **[epitope-parent-featurize](plugins/comp-immunology/skills/epitope-parent-featurize/SKILL.md)** — Given a CSV of epitope-to-parent mappings (MHC peptide ↔ source UniProt protein with occurrence range) plus a directory of parent FASTAs, compute **two** ESM C 300M embeddings per row: the epitope on its own (`direct`), and the same residues as they appear inside a window of their parent (`parent-context`). Outputs a manifest joining each input row to its `.npz` pair with slice offsets, so downstream MHC / TCR / immunogenicity models can train on either or both views.
- **[epitopes-by-source](plugins/comp-immunology/skills/epitopes-by-source/SKILL.md)** — Pull **every epitope in the IEDB from a given source organism** (a pathogen / taxon) and build a self-contained HTML report. Source resolution is forgiving — a common abbreviation (`EBV`, `CMV`, `HIV-1`, `SARS-CoV-2`, `dengue`, `M. tuberculosis`, …), a bare NCBITaxon id (`NCBITaxon:10376` / `10376`), or a substring of the IEDB organism name — and the query rolls up the **NCBITaxon subtree** (`source_organism_iri_search=cs.{IRI}`) so all strains/descendants are included. The general query layer is `iedb-query`; this is the specific "all epitopes for pathogen X" report.
- **[iedb-assay-gap](plugins/comp-immunology/skills/iedb-assay-gap/SKILL.md)** — For a given MHC allele, find epitopes in the IEDB that **have T-cell assays but no MHC-binding assay** against that allele — the "binding-data gap": peptides proven T-cell antigenic in the allele's context whose direct MHC binding has never been experimentally measured. Allele matching is grounded in the **MHC Restriction Ontology (MRO)** that the IQ-API exposes via the `mhc_allele_iri_search` array, so a serotype rollup (e.g. `HLA-A2` for `HLA-A*02:01`, `HLA-DR4` for `HLA-DRB1*04:01`) is included and every match is tagged by resolution (`exact` / `sub-allele` / `serotype`).
- **[iedb-query](plugins/comp-immunology/skills/iedb-query/SKILL.md)** — A general-purpose, efficient client **and** CLI for the IEDB Query API (IQ-API) — the public PostgREST interface to the Immune Epitope Database. Use it to count, filter, paginate, and export records from any IEDB endpoint (T-cell / MHC / B-cell assays, TCR/BCR receptors, epitopes, antigens, references) and to roll a query up an ontology. Its core trick: every search endpoint carries `<entity>_iri_search` arrays of an annotation's ancestor IRIs, so one `cs.{IRI}` filter matches a node **and all descendants** — generalized across MHC allele (MRO, incl. serotype rollup), organism (NCBITaxon), disease (DOID), antigen, and assay. Built-in keyset (cursor) pagination (the API caps responses at 10k rows), cheap `count=estimated` size gauging, retry/backoff, and CSV/JSON streaming. Importable as `import iedb` or runnable as `python iedb.py endpoints|schema|count|query|resolve`. Stdlib-only, read-only, no auth.
- **[netmhcpan](plugins/comp-immunology/skills/netmhcpan/SKILL.md)** — Run the two locally-installed **NetMHCpan** MHC **class I** predictors correctly and parse their output. NetMHCpan-4.2c (native arm64) is the bare command `netMHCpan`; NetMHCpan-4.1b (Intel via Rosetta 2) is `netMHCpan-4.1`. Predicts peptide **presentation** (eluted-ligand, EL) and optional **binding affinity** (BA) for any class I allele — from a peptide list or a protein FASTA — with per-allele `%Rank` and strong/weak (SB/WB) calls. The point of the skill is the version pitfalls that otherwise cost real time: which command is which, the Rosetta requirement, the `gawk` dependency that makes 4.1 silently write a **0-byte `.xls`**, the `.xls` `skiprows` difference (4.2 = 2, 4.1 = 1) and hyphen-vs-underscore column labels, and the rule to **never mix 4.1 and 4.2 EL scores** (the EL network was retrained). Bundles a version-aware parser that auto-detects version + format and melts wide multi-allele `.xls` into tidy long records, plus an install doctor. Class I only — for class II use NetMHCIIpan.
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
