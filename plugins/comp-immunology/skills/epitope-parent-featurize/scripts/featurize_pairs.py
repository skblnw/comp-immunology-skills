#!/usr/bin/env python3
"""
Compute ESM C 300M embeddings for (epitope, parent) pairs.

For each row of an input CSV (columns: epitope, parent UniProt accession,
occurrence range like "251-259"), emit:

  A. a direct embedding of the epitope sequence (fed to the model in
     isolation), and
  B. a parent-context embedding from a window of the parent protein
     centred on the occurrence positions, plus the 0-based slice offsets
     so the caller can extract the per-residue context vectors for the
     epitope.

Outputs land in <out>/:
    direct/<epitope>.npz            — one per unique epitope sequence
    parent/<acc>__w<ws>_<we>.npz    — one per unique parent window
    manifest.csv                    — one row per input row, joining each
                                       pair to its npz files + slice offsets
    fasta/epitopes.fasta            — input to the featurize skill (epi)
    fasta/parents.fasta             — input to the featurize skill (parent)
    featurize_direct.log,
    featurize_parent.log            — logs from the underlying featurize.py

The script does no analysis of its own — it produces the embeddings and
the bookkeeping. Downstream analysis (similarity, classification, etc.)
is the caller's responsibility.

The script auto-discovers the `featurize.py` shipped with the
`esm-featurize` skill; it does not re-implement embedding logic. Pass
`--featurize-script` to override.

CSV column auto-detection
-------------------------
The script looks for the epitope, UniProt accession, and occurrence
columns by trying common name variants. Override with --col-epitope,
--col-uniprot, --col-occurrences if your CSV uses other names.

Occurrence format
-----------------
The occurrence value must be a single 1-based inclusive range like
``251-259``. Multi-range values ("164-172, 170-178") are rejected. The
inferred slice is verified against the parent sequence — if the parent
residues at the named positions don't match the epitope, the row is
dropped (set --strict-position-check=warn to keep them with a logged
warning).
"""

import argparse
import csv
import datetime
import html
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

CANON = set("ACDEFGHIKLMNPQRSTVWY")

EPITOPE_COL_CANDIDATES = ("epitope", "peptide", "sequence", "epi", "epitope_seq")
UNIPROT_COL_CANDIDATES = ("uniprot_acc", "uniprot", "uniprot_id", "accession",
                          "parent_uniprot", "parent_uniprot_id",
                          "parent_id", "acc")
OCC_COL_CANDIDATES = ("occurrences", "occurrence", "position", "positions",
                      "range", "epitope_position")

DEFAULT_WINDOW = 480


# ── helpers ───────────────────────────────────────────────────────────
def pick_column(fieldnames: list[str], candidates: tuple[str, ...],
                override: str | None, kind: str) -> str:
    if override:
        if override not in fieldnames:
            raise SystemExit(f"--col-{kind} {override!r} not in CSV "
                             f"columns: {fieldnames}")
        return override
    for c in candidates:
        if c in fieldnames:
            return c
    raise SystemExit(
        f"could not find a column for '{kind}' in {fieldnames}. "
        f"Tried {candidates}. Use --col-{kind} to set it explicitly."
    )


def parse_occ(occ: str) -> tuple[int, int] | None:
    occ = occ.strip()
    if not occ or "," in occ or ";" in occ:
        return None
    m = re.fullmatch(r"(\d+)-(\d+)", occ)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def load_parent_seq(parents_dir: Path, acc: str) -> str | None:
    p = parents_dir / f"{acc}.fasta"
    if not p.exists():
        # also try .fa
        p2 = parents_dir / f"{acc}.fa"
        if p2.exists():
            p = p2
        else:
            return None
    with p.open() as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">"))
    return seq.upper()


def window_for(start: int, end: int, parent_len: int, window: int,
               epi_len: int) -> tuple[int, int]:
    """1-based inclusive window of length `window` centred on [start, end].

    If parent_len <= window, returns (1, parent_len) — embed the whole
    protein.
    """
    if parent_len <= window:
        return 1, parent_len
    flank = (window - epi_len) // 2
    ws = max(1, start - flank)
    we = ws + window - 1
    if we > parent_len:
        we = parent_len
        ws = we - window + 1
    return ws, we


def discover_featurize_script(override: str | None) -> Path:
    if override:
        p = Path(override).expanduser()
        if not p.exists():
            raise SystemExit(f"--featurize-script not found: {p}")
        return p
    here = Path(__file__).resolve().parent
    # 1. Next to this script (if the skill bundles its own copy).
    bundled = here / "featurize.py"
    if bundled.exists():
        return bundled
    # 2. Standalone user-level skill (legacy location).
    standalone = Path.home() / ".claude" / "skills" / "esm-featurize" / "scripts" / "featurize.py"
    if standalone.exists():
        return standalone
    # 3. Inside any installed Claude Code plugin: ~/.claude/plugins/<plugin>/skills/esm-featurize/scripts/featurize.py
    plugins_root = Path.home() / ".claude" / "plugins"
    if plugins_root.is_dir():
        hits = sorted(plugins_root.glob("*/skills/esm-featurize/scripts/featurize.py"))
        if hits:
            return hits[0]
    raise SystemExit(
        "could not find featurize.py from the esm-featurize skill. "
        "Install the esm-featurize skill (e.g. via the structural-bioinfo plugin) "
        "or pass --featurize-script."
    )


# ── row parsing / filtering ───────────────────────────────────────────
def load_rows(csv_path: Path, parents_dir: Path,
              col_epi: str, col_uni: str, col_occ: str,
              window: int,
              strict: str = "drop") -> tuple[list[dict], dict[str, int]]:
    """Returns (kept_rows, drop_counts).

    Each kept row dict carries the parsed metadata plus the parent
    sequence (so we don't re-read the FASTA later). Drop counts are
    bucketed by reason for the report.
    """
    drop = {
        "empty_epitope": 0,
        "non_canonical_epitope": 0,
        "missing_parent_fasta": 0,
        "non_canonical_parent": 0,
        "bad_occurrence": 0,
        "length_mismatch": 0,
        "position_mismatch": 0,
        "out_of_bounds": 0,
    }
    kept = []
    parent_cache: dict[str, str | None] = {}

    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            epi = (row.get(col_epi) or "").upper().strip()
            acc = (row.get(col_uni) or "").strip()
            occ_raw = (row.get(col_occ) or "").strip()
            if not epi:
                drop["empty_epitope"] += 1; continue
            if set(epi) - CANON:
                drop["non_canonical_epitope"] += 1; continue
            occ = parse_occ(occ_raw)
            if occ is None:
                drop["bad_occurrence"] += 1; continue
            start, end = occ
            if end - start + 1 != len(epi):
                drop["length_mismatch"] += 1; continue
            if acc not in parent_cache:
                parent_cache[acc] = load_parent_seq(parents_dir, acc)
            pseq = parent_cache[acc]
            if pseq is None:
                drop["missing_parent_fasta"] += 1; continue
            if set(pseq) - CANON:
                drop["non_canonical_parent"] += 1; continue
            if not (1 <= start <= end <= len(pseq)):
                drop["out_of_bounds"] += 1; continue
            if pseq[start - 1:end] != epi:
                if strict == "warn":
                    print(f"[warn] position mismatch for {epi} in {acc} "
                          f"at {start}-{end}; kept anyway", file=sys.stderr)
                else:
                    drop["position_mismatch"] += 1; continue

            kept.append({
                "epitope": epi,
                "uniprot_acc": acc,
                "start": start,
                "end": end,
                "parent_length": len(pseq),
                "parent_seq": pseq,
            })
    return kept, drop


# ── FASTA writing ─────────────────────────────────────────────────────
def parent_record_id(acc: str, ws: int, we: int) -> str:
    return f"{acc}__w{ws}_{we}"


def write_fastas(rows: list[dict], out_dir: Path, window: int,
                 epi_len_default: int = 9) -> tuple[Path, Path, dict]:
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    epi_fa = fasta_dir / "epitopes.fasta"
    par_fa = fasta_dir / "parents.fasta"

    # Unique epitope sequences
    epis: dict[str, str] = {}
    for r in rows:
        epis.setdefault(r["epitope"], r["epitope"])
    with epi_fa.open("w") as f:
        for epi in epis:
            f.write(f">{epi}\n{epi}\n")

    # Unique (acc, window) parent records
    par_records: dict[str, str] = {}
    row_to_pid: dict[tuple[str, str, int, int], str] = {}
    for r in rows:
        ws, we = window_for(r["start"], r["end"], r["parent_length"],
                            window=window, epi_len=len(r["epitope"]))
        pid = parent_record_id(r["uniprot_acc"], ws, we)
        par_records.setdefault(pid, r["parent_seq"][ws - 1:we])
        row_to_pid[(r["epitope"], r["uniprot_acc"], r["start"], r["end"])] = pid
    with par_fa.open("w") as f:
        for pid, seq in par_records.items():
            f.write(f">{pid}\n{seq}\n")

    return epi_fa, par_fa, row_to_pid


# ── featurize.py drivers ──────────────────────────────────────────────
def run_featurize(featurize: Path, fasta: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(featurize), "--fasta", str(fasta),
           "--out", str(out_dir)]
    print(f"[run] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def run_featurize_parallel(featurize: Path,
                           fasta_a: Path, out_a: Path,
                           fasta_b: Path, out_b: Path,
                           log_a: Path, log_b: Path) -> None:
    out_a.mkdir(parents=True, exist_ok=True)
    out_b.mkdir(parents=True, exist_ok=True)
    cmd_a = [sys.executable, str(featurize), "--fasta", str(fasta_a),
             "--out", str(out_a)]
    cmd_b = [sys.executable, str(featurize), "--fasta", str(fasta_b),
             "--out", str(out_b)]
    print(f"[run-a] {' '.join(cmd_a)} > {log_a}", flush=True)
    print(f"[run-b] {' '.join(cmd_b)} > {log_b}", flush=True)
    with open(log_a, "w") as fa, open(log_b, "w") as fb:
        pa = subprocess.Popen(cmd_a, stdout=fa, stderr=subprocess.STDOUT)
        pb = subprocess.Popen(cmd_b, stdout=fb, stderr=subprocess.STDOUT)
        rc_a = pa.wait()
        rc_b = pb.wait()
    if rc_a != 0:
        raise RuntimeError(f"direct featurize failed (rc={rc_a}); see {log_a}")
    if rc_b != 0:
        raise RuntimeError(f"parent featurize failed (rc={rc_b}); see {log_b}")


# ── manifest ──────────────────────────────────────────────────────────
def write_manifest(rows: list[dict], row_to_pid: dict, window: int,
                   direct_dir: Path, parent_dir: Path, manifest_path: Path):
    """One row per input row joining (epitope, parent occurrence) to the
    two npz files and the 0-based slice into the parent window's
    per_residue array.
    """
    with manifest_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "epitope", "uniprot_acc", "occurrence_start", "occurrence_end",
            "parent_length",
            "window_start", "window_end",
            "slice_start_0based", "slice_end_0based",
            "direct_npz", "parent_npz",
        ])
        for r in rows:
            key = (r["epitope"], r["uniprot_acc"], r["start"], r["end"])
            pid = row_to_pid[key]
            ws = int(pid.split("__w")[1].split("_")[0])
            we = int(pid.rsplit("_", 1)[1])
            ss = r["start"] - ws        # 0-based into parent window
            se = ss + len(r["epitope"])
            direct_npz = direct_dir / f"{r['epitope']}.npz"
            parent_npz = parent_dir / f"{pid}.npz"
            w.writerow([
                r["epitope"], r["uniprot_acc"], r["start"], r["end"],
                r["parent_length"],
                ws, we, ss, se,
                str(direct_npz), str(parent_npz),
            ])


def write_summary(out_dir: Path, rows: list[dict], drop: dict,
                  n_unique_epis: int, n_unique_par_windows: int,
                  window: int, args) -> None:
    summary = {
        "n_input_rows_kept": len(rows),
        "n_unique_epitopes": n_unique_epis,
        "n_unique_parent_windows": n_unique_par_windows,
        "window_size": window,
        "drop_counts": drop,
        "csv": str(Path(args.csv).resolve()),
        "parents_dir": str(Path(args.parents_dir).resolve()),
        "out": str(out_dir.resolve()),
        "sequential": args.sequential,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    # Markdown report
    md = [
        "# epitope-parent-featurize summary",
        "",
        f"- input CSV: `{summary['csv']}`",
        f"- parents dir: `{summary['parents_dir']}`",
        f"- output: `{summary['out']}`",
        f"- window size: {window} residues (centred on the occurrence "
        f"for parents longer than that; full protein otherwise)",
        "",
        "## Row accounting",
        "",
        f"- kept: **{len(rows)}**",
        f"- unique epitope sequences: {n_unique_epis}",
        f"- unique parent windows:    {n_unique_par_windows}",
        "",
        "### Dropped rows by reason",
        "",
        "| reason | count |",
        "|---|---:|",
    ]
    for k, v in drop.items():
        md.append(f"| {k} | {v} |")
    md += [
        "",
        "## Files",
        "",
        "- `manifest.csv` — one row per input row, with the slice offsets",
        "- `direct/<epitope>.npz` — direct ESM C 300M embedding of each "
        "unique 9-mer (or other length)",
        "- `parent/<acc>__w<ws>_<we>.npz` — parent-context window "
        "embedding; per_residue has shape (window_len, 960)",
        "- `fasta/{epitopes,parents}.fasta` — what was fed to the "
        "underlying featurize.py",
        "- `featurize_{direct,parent}.log` — logs from each model run",
        "",
        "## Loading the per-residue context vector for an epitope",
        "",
        "```python",
        "import csv, numpy as np",
        "manifest = list(csv.DictReader(open('manifest.csv')))",
        "row = manifest[0]",
        "d = np.load(row['parent_npz'])",
        "context_per_res = d['per_residue'][int(row['slice_start_0based'])"
        ":int(row['slice_end_0based'])]   # (L, 960)",
        "direct = np.load(row['direct_npz'])['per_residue']               "
        "      # (L, 960)",
        "```",
        "",
    ]
    (out_dir / "report.md").write_text("\n".join(md))
    write_html_report(out_dir, rows, drop, n_unique_epis,
                      n_unique_par_windows, window, summary)


def write_html_report(out_dir: Path, rows: list[dict], drop: dict,
                      n_unique_epis: int, n_unique_par_windows: int,
                      window: int, summary: dict) -> None:
    """Self-contained HTML twin of report.md. No JS, no external assets."""
    total_seen = len(rows) + sum(drop.values())

    def esc(x):
        return html.escape(str(x))

    drop_rows_html = ""
    for k, v in drop.items():
        pct = (100 * v / total_seen) if total_seen else 0
        drop_rows_html += (
            f"<tr><td>{esc(k.replace('_', ' '))}</td>"
            f"<td>{v}</td>"
            f"<td>{pct:.1f}%</td></tr>"
        )
    kept_pct = (100 * len(rows) / total_seen) if total_seen else 0
    drop_rows_html += (
        f"<tr class='kept-row'><td><b>kept</b></td>"
        f"<td><b>{len(rows)}</b></td>"
        f"<td><b>{kept_pct:.1f}%</b></td></tr>"
    )

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    mode = "sequential" if summary["sequential"] else "parallel"

    load_snippet = (
        "import csv, numpy as np\n"
        "manifest = list(csv.DictReader(open('manifest.csv')))\n"
        "row = manifest[0]\n"
        "# parent-context per-residue vectors at the epitope's positions\n"
        "p = np.load(row['parent_npz'])\n"
        "ss = int(row['slice_start_0based'])\n"
        "se = int(row['slice_end_0based'])\n"
        "context_per_res = p['per_residue'][ss:se]      # (L, 960)\n"
        "context_mean    = context_per_res.mean(axis=0) # (960,)\n"
        "# direct (peptide-only) embedding\n"
        "d = np.load(row['direct_npz'])\n"
        "direct_per_res  = d['per_residue']             # (L, 960)\n"
        "direct_mean     = d['mean_pooled']             # (960,)"
    )

    html_doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>epitope-parent-featurize — run summary</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 880px; margin: 2rem auto; padding: 0 1rem; color: #222;
         line-height: 1.45; }}
  h1 {{ border-bottom: 2px solid #eee; padding-bottom: .3rem; }}
  h2 {{ margin-top: 2rem; border-bottom: 1px solid #eee; padding-bottom: .2rem; }}
  table {{ border-collapse: collapse; margin: .8rem 0; width: 100%; font-size: 0.93rem; }}
  th, td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: right; }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ background: #f6f6f6; }}
  tr:nth-child(even) td {{ background: #fafafa; }}
  tr.kept-row td {{ background: #ecfdf5; }}
  code {{ background: #f3f3f3; padding: 1px 4px; border-radius: 3px;
          font-size: 0.92em; word-break: break-all; }}
  pre {{ background: #0f172a; color: #e2e8f0; padding: .8rem 1rem;
         border-radius: 6px; overflow-x: auto; font-size: 0.85rem;
         line-height: 1.4; }}
  .kpi {{ display: flex; gap: 1rem; flex-wrap: wrap; margin: 1rem 0; }}
  .kpi div {{ flex: 1 1 160px; background: #f6f8fb; border: 1px solid #e3e8ef;
              border-radius: 6px; padding: .6rem .8rem; }}
  .kpi div b {{ display: block; font-size: 1.4rem; }}
  .kpi div small {{ color: #555; }}
  .meta {{ background: #fafafa; border: 1px solid #eee; border-radius: 6px;
           padding: .6rem .8rem; font-size: 0.9rem; }}
  .meta dt {{ font-weight: 600; color: #444; margin-top: .35rem; }}
  .meta dd {{ margin-left: 0; word-break: break-all; }}
  .note {{ background: #fff4e5; border-left: 4px solid #f59e0b;
           padding: .6rem .8rem; border-radius: 4px; margin: 1rem 0;
           font-size: 0.93rem; }}
  ul li {{ margin: .15rem 0; }}
</style></head>
<body>
<h1>epitope-parent-featurize — run summary</h1>
<p>Dual ESM C 300M embeddings for each (epitope, parent) pair: a
<b>direct</b> embedding of the peptide alone, plus a <b>parent-context</b>
embedding from a {window}-residue window centred on the occurrence
(or the whole parent if it is shorter than {window} aa). The
<code>manifest.csv</code> joins each input row to both <code>.npz</code>
files and includes the 0-based slice offsets into the parent window so
the caller can recover per-residue context vectors.</p>

<div class="kpi">
  <div><b>{len(rows)}</b><small>input rows kept</small></div>
  <div><b>{n_unique_epis}</b><small>unique epitope sequences</small></div>
  <div><b>{n_unique_par_windows}</b><small>unique parent windows</small></div>
  <div><b>{window} aa</b><small>context window size</small></div>
</div>

<h2>Run metadata</h2>
<div class="meta"><dl>
  <dt>input CSV</dt><dd><code>{esc(summary['csv'])}</code></dd>
  <dt>parents directory</dt><dd><code>{esc(summary['parents_dir'])}</code></dd>
  <dt>output directory</dt><dd><code>{esc(summary['out'])}</code></dd>
  <dt>featurize mode</dt><dd>{mode}</dd>
  <dt>report generated</dt><dd>{now}</dd>
</dl></div>

<h2>Row accounting</h2>
<p>Total CSV rows seen: <b>{total_seen}</b>. Per-reason breakdown:</p>
<table>
  <thead><tr><th>reason</th><th>count</th><th>% of seen</th></tr></thead>
  <tbody>{drop_rows_html}</tbody>
</table>

<h2>Files produced</h2>
<ul>
  <li><code>manifest.csv</code> — one row per input row, with the slice offsets (read this first).</li>
  <li><code>direct/&lt;epitope&gt;.npz</code> — direct ESM C 300M embedding of each unique epitope sequence. Contains <code>per_residue</code> ((L, 960) float32), <code>mean_pooled</code> ((960,) float32), <code>sequence</code>, <code>id</code>, <code>model</code>.</li>
  <li><code>parent/&lt;acc&gt;__w&lt;ws&gt;_&lt;we&gt;.npz</code> — parent-context window embedding; <code>per_residue</code> has shape (window_len, 960). The window's start/end (1-based, inclusive) is encoded in the filename.</li>
  <li><code>fasta/epitopes.fasta</code>, <code>fasta/parents.fasta</code> — what was fed to the underlying <code>featurize.py</code>.</li>
  <li><code>featurize_direct.log</code>, <code>featurize_parent.log</code> — logs from each model run.</li>
  <li><code>summary.json</code> — same row-accounting numbers as above, machine-readable.</li>
  <li><code>report.md</code>, <code>report.html</code> — this report.</li>
</ul>

<h2>Loading the per-residue context vector for an epitope</h2>
<pre>{esc(load_snippet)}</pre>

<div class="note">
<b>Why a window and not the full protein?</b> ESM C's stated max input is
2048 residues; embedding a multi-kilobase polyprotein wastes compute on
residues far from the epitope and risks model edge effects. A
{window}-residue window centred on the occurrence keeps every
"context" embedding the same length so downstream comparisons are
apples-to-apples; parents shorter than the window are embedded in
full.<br><br>
<b>Slice verification.</b> The script always checks that
<code>parent[start-1:end] == epitope</code> before writing the
manifest. Silent off-by-ones produce embeddings that <i>look</i>
fine but are sliced at the wrong residues — a uniquely insidious
failure mode for this kind of workflow.
</div>
</body></html>
"""
    (out_dir / "report.html").write_text(html_doc)


# ── main ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=__doc__)
    ap.add_argument("--csv", required=True, help="input CSV path")
    ap.add_argument("--parents-dir", required=True,
                    help="directory of <ACC>.fasta files for each parent")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--col-epitope", default=None,
                    help="column name for the epitope sequence "
                         f"(default: auto-detect, tried {EPITOPE_COL_CANDIDATES})")
    ap.add_argument("--col-uniprot", default=None,
                    help="column name for the parent UniProt accession")
    ap.add_argument("--col-occurrences", default=None,
                    help="column name for the occurrence range (e.g. '251-259')")
    ap.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                    help=f"window size in residues for long parents (default {DEFAULT_WINDOW})")
    ap.add_argument("--n", type=int, default=0,
                    help="if >0, sample N rows after filtering (for smoke "
                         "tests). Default 0 = all eligible rows.")
    ap.add_argument("--seed", type=int, default=20260517,
                    help="RNG seed for --n sampling")
    ap.add_argument("--position-check", choices=["drop", "warn"], default="drop",
                    help="what to do when parent[start-1:end] != epitope "
                         "(default: drop the row)")
    ap.add_argument("--skip-featurize", action="store_true",
                    help="reuse existing npz outputs; only rebuild manifest+report")
    ap.add_argument("--sequential", action="store_true",
                    help="run direct + parent featurize sequentially "
                         "(default: parallel)")
    ap.add_argument("--featurize-script", default=None,
                    help="path to esm-featurize's featurize.py "
                         "(default: auto-discover)")
    ap.add_argument("--no-report", action="store_true",
                    help="skip writing summary.json / report.md / report.html")
    args = ap.parse_args()

    csv_path = Path(args.csv).resolve()
    parents_dir = Path(args.parents_dir).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        raise SystemExit(f"CSV not found: {csv_path}")
    if not parents_dir.is_dir():
        raise SystemExit(f"parents dir not found: {parents_dir}")

    # Sniff CSV header to find columns
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
    col_epi = pick_column(fieldnames, EPITOPE_COL_CANDIDATES, args.col_epitope, "epitope")
    col_uni = pick_column(fieldnames, UNIPROT_COL_CANDIDATES, args.col_uniprot, "uniprot")
    col_occ = pick_column(fieldnames, OCC_COL_CANDIDATES, args.col_occurrences, "occurrences")
    print(f"[1/5] reading rows from {csv_path}", flush=True)
    print(f"      epitope col: {col_epi}", flush=True)
    print(f"      uniprot col: {col_uni}", flush=True)
    print(f"      occ col:     {col_occ}", flush=True)

    rows, drop = load_rows(csv_path, parents_dir,
                           col_epi, col_uni, col_occ,
                           window=args.window,
                           strict=args.position_check)
    print(f"      kept: {len(rows)} rows after filtering", flush=True)
    if drop:
        print(f"      drops: {drop}", flush=True)

    if args.n and args.n < len(rows):
        rng = random.Random(args.seed)
        rng.shuffle(rows)
        rows = rows[:args.n]
        print(f"      sampled down to {len(rows)} rows (seed={args.seed})", flush=True)

    if not rows:
        raise SystemExit("no rows left after filtering; nothing to do.")

    print(f"[2/5] writing FASTAs to {out_dir/'fasta'}", flush=True)
    epi_fa, par_fa, row_to_pid = write_fastas(rows, out_dir, args.window)
    n_unique_epis = sum(1 for _ in open(epi_fa)) // 2
    n_unique_pwin = sum(1 for _ in open(par_fa)) // 2
    print(f"      unique epitope records: {n_unique_epis}", flush=True)
    print(f"      unique parent windows:  {n_unique_pwin}", flush=True)

    direct_dir = out_dir / "direct"
    parent_dir = out_dir / "parent"
    if not args.skip_featurize:
        featurize = discover_featurize_script(args.featurize_script)
        print(f"[3-4/5] running featurize ({'sequential' if args.sequential else 'parallel'})",
              flush=True)
        if args.sequential:
            run_featurize(featurize, epi_fa, direct_dir)
            run_featurize(featurize, par_fa, parent_dir)
        else:
            run_featurize_parallel(
                featurize,
                epi_fa, direct_dir,
                par_fa, parent_dir,
                out_dir / "featurize_direct.log",
                out_dir / "featurize_parent.log",
            )
    else:
        print("[3-4/5] --skip-featurize: reusing existing npz dirs", flush=True)

    manifest_path = out_dir / "manifest.csv"
    print(f"[5/5] writing manifest to {manifest_path}", flush=True)
    write_manifest(rows, row_to_pid, args.window,
                   direct_dir, parent_dir, manifest_path)

    if not args.no_report:
        write_summary(out_dir, rows, drop, n_unique_epis, n_unique_pwin,
                      args.window, args)
        print(f"      summary: {out_dir/'summary.json'}", flush=True)
        print(f"      report:  {out_dir/'report.md'} + {out_dir/'report.html'}",
              flush=True)
    print(f"done. {out_dir}")


if __name__ == "__main__":
    main()
