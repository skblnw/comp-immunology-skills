#!/usr/bin/env python3
"""
parse_netmhcpan.py — robust, version-aware parser for NetMHCpan output.

The single biggest footgun when reading NetMHCpan output programmatically is that
the `.xls` layout differs between versions:

  * NetMHCpan-4.2  : line 1 is a `#<command>` line, line 2 is the allele-group row,
                     line 3 is the column header  -> pandas needs `skiprows=2`.
                     Column labels use underscores: EL_score / EL_rank / BA_score / BA_rank.
                     Position column starts at 1.
  * NetMHCpan-4.1  : NO command line; line 1 is the allele-group row, line 2 is the
                     header                          -> pandas needs `skiprows=1`.
                     Column labels use hyphens/mixed case: EL-score / EL_Rank / BA-score / BA_Rank.
                     Position column starts at 0.

Multi-allele `.xls` files are also "wide": each requested allele gets its own block
of columns (core, icore, EL, BA...), with the allele name sitting in the skipped
allele-group row above its block. This script auto-detects the version, reads that
allele-group row, and melts the wide table into tidy long records — one row per
(peptide, allele) — so you never have to hand-align columns again.

It also parses the space-aligned `.out` (stdout) format, where the allele is already
on each row. Format (.xls vs .out) is auto-detected from the presence of tab characters.

Usage (CLI):
    python parse_netmhcpan.py OUT.xls                     # tidy CSV to stdout
    python parse_netmhcpan.py a.xls b.out --format json   # JSON to stdout
    python parse_netmhcpan.py OUT.xls --out tidy.csv      # write CSV file
    python parse_netmhcpan.py OUT.xls --sb 0.5 --wb 2.0   # custom %Rank_EL bind thresholds

Usage (import):
    from parse_netmhcpan import parse_file
    records = parse_file("OUT.xls")        # list[dict]
    # optional: df = to_dataframe(records)  # pandas, if installed

Output columns (tidy/long):
    pos, peptide, id, allele, core, icore,
    el_score, el_rank, ba_score, ba_rank, aff_nM,
    bind_level, version, source_file

`bind_level` is taken from the .out annotation when present; otherwise it is derived
from %Rank_EL using the default NetMHCpan thresholds (SB <= 0.5, WB <= 2.0). Override
with --sb/--wb if you ran NetMHCpan with non-default -rth/-rlt.

Stdlib only.
"""
import sys, csv, json, argparse

# canonical sub-column names within an allele block, keyed by normalized header token
_SUBCOL = {
    "core": "core",
    "icore": "icore",
    "el_score": "el_score", "elscore": "el_score",
    "el_rank": "el_rank",
    "ba_score": "ba_score", "bascore": "ba_score",
    "ba_rank": "ba_rank",
}
_TRAILING = {"ave", "nb"}          # per-peptide summary columns at the far right of an .xls
_LEADING = {"pos", "peptide", "id"}


def _norm(tok):
    """Normalize a header token so 4.1 'EL-score'/'BA_Rank' and 4.2 'EL_score'/'BA_rank' match."""
    return tok.strip().lower().replace("-", "_").replace(".", "_")


def _bind_level(el_rank, sb, wb):
    if el_rank is None:
        return ""
    if el_rank <= sb:
        return "SB"
    if el_rank <= wb:
        return "WB"
    return ""


def _f(v):
    """float or None."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _detect_format(text):
    # The header row (the one naming "Peptide") is tab-delimited in .xls but
    # space-aligned in .out. NetMHCpan .out preambles can contain stray tabs, so
    # keying off the header line is far more reliable than "any tab in file".
    for ln in text.splitlines():
        if "Peptide" in ln:
            return "xls" if ln.count("\t") >= 2 else "out"
    return "xls" if any(ln.count("\t") >= 3 for ln in text.splitlines()) else "out"


def _detect_xls_version(lines):
    return "4.2" if lines and lines[0].startswith("#") else "4.1"


def parse_xls(text, sb=0.5, wb=2.0, source_file=""):
    lines = [ln for ln in text.splitlines() if ln != ""]
    version = _detect_xls_version(lines)
    start = 1 if version == "4.2" else 0          # index of the allele-group row
    allele_row = lines[start].split("\t")
    header = lines[start + 1].split("\t")
    data = lines[start + 2:]

    # alleles sit at the column index where their block begins
    alleles = [(i, tok.strip()) for i, tok in enumerate(allele_row) if tok.strip()]
    if not alleles:
        raise ValueError("No alleles found in allele-group row; is this a NetMHCpan .xls?")

    first = alleles[0][0]
    # trailing summary columns (Ave, NB) mark the end of the allele blocks
    trailing_start = len(header)
    for i in range(len(header) - 1, first - 1, -1):
        if _norm(header[i]) in _TRAILING:
            trailing_start = i
        else:
            break
    width = (trailing_start - first) // len(alleles)

    # leading shared columns (pos/peptide/id) by name
    lead_idx = {}
    for i in range(first):
        n = _norm(header[i])
        if n in _LEADING:
            lead_idx[n] = i

    records = []
    for row in data:
        f = row.split("\t")
        if len(f) < trailing_start:
            continue
        pos = f[lead_idx["pos"]].strip() if "pos" in lead_idx else ""
        pep = f[lead_idx["peptide"]].strip() if "peptide" in lead_idx else ""
        rid = f[lead_idx["id"]].strip() if "id" in lead_idx else ""
        if not pep:
            continue
        for ai, (col0, allele) in enumerate(alleles):
            block_start = first + ai * width
            rec = {"pos": pos, "peptide": pep, "id": rid, "allele": allele,
                   "core": None, "icore": None, "el_score": None, "el_rank": None,
                   "ba_score": None, "ba_rank": None, "aff_nM": None}
            for j in range(width):
                canon = _SUBCOL.get(_norm(header[block_start + j]))
                if canon is None:
                    continue
                val = f[block_start + j].strip()
                rec[canon] = val if canon in ("core", "icore") else _f(val)
            rec["bind_level"] = _bind_level(rec["el_rank"], sb, wb)
            rec["version"] = version
            rec["source_file"] = source_file
            records.append(rec)
    return records


def parse_out(text, sb=0.5, wb=2.0, source_file=""):
    # leading schema is fixed across versions; BA block is present only with -BA
    LEAD = 11   # Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity
    version = ""
    has_ba = False
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("# NetMHCpan version"):
            version = s.split("version", 1)[1].strip()
        if s.startswith("Pos") and "Peptide" in s:
            has_ba = "Score_BA" in s

    records = []
    for ln in text.splitlines():
        s = ln.rstrip("\n")
        t = s.split()
        # data rows: first token is an int position, second looks like an allele (has a '*' or '-')
        if len(t) < LEAD + 2 or not t[0].isdigit():
            continue
        if not (("*" in t[1]) or t[1].startswith("HLA") or t[1].startswith("H-2")
                or "-" in t[1] or ":" in t[1]):
            continue
        allele = t[1].replace("*", "")
        rec = {"pos": t[0], "peptide": t[2], "id": t[10], "allele": allele,
               "core": t[3], "icore": t[9],
               "el_score": _f(t[LEAD]), "el_rank": _f(t[LEAD + 1]),
               "ba_score": None, "ba_rank": None, "aff_nM": None}
        idx = LEAD + 2
        if has_ba and len(t) >= idx + 3:
            rec["ba_score"] = _f(t[idx]); rec["ba_rank"] = _f(t[idx + 1]); rec["aff_nM"] = _f(t[idx + 2])
            idx += 3
        rest = " ".join(t[idx:])
        if "SB" in rest:
            rec["bind_level"] = "SB"
        elif "WB" in rest:
            rec["bind_level"] = "WB"
        else:
            rec["bind_level"] = _bind_level(rec["el_rank"], sb, wb)
        rec["version"] = version
        rec["source_file"] = source_file
        records.append(rec)
    return records


def parse_file(path, fmt=None, sb=0.5, wb=2.0):
    with open(path, errors="ignore") as fh:
        text = fh.read()
    fmt = fmt or _detect_format(text)
    if fmt == "xls":
        return parse_xls(text, sb, wb, source_file=path)
    return parse_out(text, sb, wb, source_file=path)


def to_dataframe(records):
    """Convenience: return a pandas DataFrame (pandas must be installed)."""
    import pandas as pd
    return pd.DataFrame(records)


COLUMNS = ["pos", "peptide", "id", "allele", "core", "icore",
           "el_score", "el_rank", "ba_score", "ba_rank", "aff_nM",
           "bind_level", "version", "source_file"]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Parse NetMHCpan .xls/.out into tidy long records (auto-detects 4.1 vs 4.2).")
    ap.add_argument("files", nargs="+", help="NetMHCpan .xls or .out file(s)")
    ap.add_argument("--format", choices=["auto", "xls", "out"], default="auto")
    ap.add_argument("--out", help="write to this file instead of stdout")
    ap.add_argument("--json", dest="as_json", action="store_true", help="emit JSON instead of CSV")
    ap.add_argument("--sb", type=float, default=0.5, help="%%Rank_EL strong-binder threshold (default 0.5)")
    ap.add_argument("--wb", type=float, default=2.0, help="%%Rank_EL weak-binder threshold (default 2.0)")
    args = ap.parse_args(argv)

    fmt = None if args.format == "auto" else args.format
    records = []
    for p in args.files:
        records.extend(parse_file(p, fmt=fmt, sb=args.sb, wb=args.wb))

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    try:
        if args.as_json:
            json.dump(records, out, indent=1)
            out.write("\n")
        else:
            w = csv.DictWriter(out, fieldnames=COLUMNS)
            w.writeheader()
            for r in records:
                w.writerow(r)
    finally:
        if args.out:
            out.close()
    sys.stderr.write(f"parsed {len(records)} records from {len(args.files)} file(s)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
