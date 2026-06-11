#!/usr/bin/env python3
# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""``iedb`` CLI -- efficient client for the IEDB Query API (IQ-API).

Subcommands: ``endpoints`` (discover grain + cursor keys), ``schema`` (an endpoint's columns),
``count`` (estimated/exact), ``query`` (auto keyset-paginated, streamed CSV/JSON), and ``resolve``
(name -> IRI + ancestor chain + serotype). The client lives in :mod:`iedb_toolkit.core`; this module
is only argument parsing + dispatch, wrapped by :func:`iedb_toolkit.cli._run.run`.
"""

import argparse
import json
import sys

import core
from core import *  # noqa: F401,F403 — re-export client API for `import iedb`
__all__ = list(core.__all__)
from _run import run


def _parse_where(where_list):
    """['col=op.value', ...] -> {col: 'op.value'}. Supports 'or=(a.eq.1,b.eq.2)' (key 'or')."""
    filters = {}
    for w in where_list or []:
        if "=" not in w:
            raise ValueError(f"--where must be 'col=op.value' (e.g. structure_id=eq.123), got {w!r}")
        k, v = w.split("=", 1)
        filters[k.strip()] = v.strip()
    return filters


def _main(argv=None):
    ap = argparse.ArgumentParser(
        prog="iedb",
        description="Efficient client/CLI for the IEDB Query API (IQ-API).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("endpoints", help="list main endpoints, their grain, and cursor keys")

    p_sc = sub.add_parser("schema", help="show an endpoint's columns (one sample row)")
    p_sc.add_argument("endpoint")
    p_sc.add_argument("--timeout", type=int, default=90)

    p_ct = sub.add_parser("count", help="count matching rows (estimated by default)")
    p_ct.add_argument("endpoint")
    p_ct.add_argument("--where", action="append", default=[], help="PostgREST filter col=op.value (repeatable)")
    p_ct.add_argument("--exact", action="store_true", help="exact count (slower) instead of estimated")
    p_ct.add_argument("--timeout", type=int, default=90)

    p_q = sub.add_parser("query", help="query an endpoint, auto keyset-paginated and streamed")
    p_q.add_argument("endpoint")
    p_q.add_argument("--where", action="append", default=[], help="PostgREST filter col=op.value (repeatable)")
    p_q.add_argument("--select", help="comma-separated columns (default: all)")
    p_q.add_argument("--order", help="cursor column (default: endpoint's primary key)")
    p_q.add_argument("--max-rows", type=int, default=None, help="cap rows fetched (logs if it truncates)")
    p_q.add_argument("--format", choices=["json", "csv"], default="csv")
    p_q.add_argument("--out", help="output file (default: stdout)")
    p_q.add_argument("--timeout", type=int, default=90)

    p_r = sub.add_parser("resolve", help="resolve an entity name to its IRI + ancestor chain")
    p_r.add_argument("kind", choices=sorted(core.ONTOLOGY))
    p_r.add_argument("name")
    p_r.add_argument("--timeout", type=int, default=90)

    args = ap.parse_args(argv)

    if args.cmd == "endpoints":
        print(f"{'endpoint':16} {'cursor key':28} grain")
        print("-" * 78)
        for ep, grain in core.ENDPOINTS.items():
            print(f"{ep:16} {core.CURSOR_KEYS.get(ep, '(needs --order)'):28} {grain}")
        print("\nOntology kinds for `resolve` / subtree filters:", ", ".join(sorted(core.ONTOLOGY)))
        return 0

    if args.cmd == "schema":
        rec = core.fetch_one(args.endpoint, {}, select="*", timeout=args.timeout)
        if not rec:
            sys.stderr.write(f"no rows returned for {args.endpoint}\n")
            return 1
        print(f"{args.endpoint}: {len(rec)} columns")
        for k in sorted(rec):
            sample = repr(rec[k])
            print(f"  {k:42} {sample[:60]}")
        return 0

    if args.cmd == "count":
        filters = _parse_where(args.where)
        n = core.count(args.endpoint, filters, mode="exact" if args.exact else "estimated",
                       timeout=args.timeout)
        mode = "exact" if args.exact else "estimated"
        print(f"{n}   ({mode})")
        return 0

    if args.cmd == "query":
        filters = _parse_where(args.where)
        try:
            est = core.count(args.endpoint, filters, mode="estimated", timeout=args.timeout)
            sys.stderr.write(f"  estimated rows: {est:,}\n")
            if args.max_rows and est > args.max_rows:
                sys.stderr.write(
                    f"  NOTE: --max-rows={args.max_rows:,} will cap the pull "
                    f"(~{est - args.max_rows:,} matching rows NOT fetched)\n")
        except Exception as e:  # count is advisory; never block the pull
            sys.stderr.write(f"  (size preflight skipped: {e})\n")

        rows = core.iter_all(args.endpoint, filters, select=args.select, order_key=args.order,
                             timeout=args.timeout, max_rows=args.max_rows)
        out_fh = open(args.out, "w", newline="") if args.out else sys.stdout
        try:
            core.write_rows(rows, args.format, out_fh)
        finally:
            if args.out:
                out_fh.close()
                sys.stderr.write(f"  wrote {args.out}\n")
        return 0

    if args.cmd == "resolve":
        try:
            info = core.resolve_entity(args.kind, args.name, timeout=args.timeout)
        except ValueError as e:
            sys.stderr.write(f"error: {e}\n")
            return 1
        print(json.dumps(info, indent=2, default=str))
        if args.kind == "allele" and info.get("serotype_name"):
            sys.stderr.write(
                f"\n  subtree filter (allele + sub-alleles): "
                f"--where '{core.ONTOLOGY['allele'][0]}=cs.{{{info['iri']}}}'\n"
                f"  serotype-level filter:                  "
                f"--where 'mhc_allele_name=eq.{info['serotype_name']}'\n")
        return 0

    ap.error("unknown command")


def main(argv=None):
    return run(_main, argv)


if __name__ == "__main__":
    sys.exit(main())
