#!/usr/bin/env python3
# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""epitopes_by_source.py — vendored CLI entry. Domain logic in report.py; IQ-API client in core.py."""
import sys

from _run import run
from report import main

if __name__ == "__main__":
    sys.exit(run(main))
