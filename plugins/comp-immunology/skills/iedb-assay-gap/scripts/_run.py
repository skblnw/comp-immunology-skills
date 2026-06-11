# vendored from the iedb-toolkit package — DO NOT EDIT HERE.
# Edit src/iedb_toolkit/ in the iedb-toolkit repo and re-run tools/vendor.py.
"""Shared CLI error wrapper.

All three console scripts route their ``main`` through :func:`run`, which turns the IQ-API's
error modes into one friendly stderr line + exit code instead of a traceback:

  * PostgREST 4xx (bad column/filter/endpoint) -> its ``{message, hint, details}`` body.
  * network/DNS/timeout -> a "could not reach" line.
  * expected ValueError/RuntimeError (user/API mistakes) -> a clean ``error: ...`` line.
  * a broken pipe (e.g. ``... | head``) -> silent success; Ctrl-C -> 130.
"""

import json
import sys
import urllib.error


def run(main_fn, argv=None):
    """Invoke ``main_fn(argv)`` and translate its exceptions into stderr + an exit code."""
    try:
        return main_fn(argv)
    except urllib.error.HTTPError as exc:            # 4xx from PostgREST (bad column/filter/endpoint)
        body = ""
        try:
            body = exc.read().decode("utf-8", "replace").strip()
        except Exception:
            pass
        try:                                          # PostgREST returns {message, hint, details}
            j = json.loads(body)
            body = "; ".join(str(j[k]) for k in ("message", "hint", "details") if j.get(k))
        except Exception:
            pass
        sys.stderr.write(f"error: HTTP {exc.code} from IEDB API: {body or exc.reason}\n")
        return 1
    except urllib.error.URLError as exc:             # network/DNS/timeout
        sys.stderr.write(f"error: could not reach IEDB API: {exc.reason}\n")
        return 1
    except (ValueError, RuntimeError) as exc:        # expected user/API errors -> clean message
        sys.stderr.write(f"error: {exc}\n")
        return 1
    except BrokenPipeError:                          # e.g. `query ... | head`
        return 0
    except KeyboardInterrupt:
        return 130
