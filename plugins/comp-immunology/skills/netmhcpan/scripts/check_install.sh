#!/usr/bin/env bash
# check_install.sh — verify the local NetMHCpan installs actually run.
#
# Checks both versions end-to-end: command on PATH, prerequisites, and a live
# prediction on a couple of peptides. Portable — uses the PATH commands, not
# hardcoded install paths. Override the command names with env vars if yours differ:
#   NETMHCPAN_42_CMD (default: netMHCpan)   NETMHCPAN_41_CMD (default: netMHCpan-4.1)
#
# Exit 0 if every install that is present passes; non-zero if a present install is broken.

set -u
CMD42="${NETMHCPAN_42_CMD:-netMHCpan}"
CMD41="${NETMHCPAN_41_CMD:-netMHCpan-4.1}"
GREEN=$'\033[32m'; RED=$'\033[31m'; YEL=$'\033[33m'; DIM=$'\033[2m'; RST=$'\033[0m'
ok(){ echo "  ${GREEN}✓${RST} $*"; }
bad(){ echo "  ${RED}✗${RST} $*"; }
warn(){ echo "  ${YEL}!${RST} $*"; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
printf '%s\n' NLVPMVATV GILGFVFTL > "$TMP/p.pep"   # CMV pp65 + flu M1: both strong HLA-A*02:01 binders
rc=0

echo "Host arch: $(uname -m)"
echo "Prerequisites:"
[ -x /bin/tcsh ] && ok "/bin/tcsh present" || { bad "/bin/tcsh missing (both launchers are tcsh scripts)"; rc=1; }
if command -v gawk >/dev/null 2>&1; then ok "gawk present ($(command -v gawk)) — needed by 4.1 -xls"; else warn "gawk missing — NetMHCpan-4.1 '-xls' writes a 0-byte file without it (brew install gawk)"; fi
if [ "$(uname -m)" = "arm64" ]; then
  if /usr/bin/pgrep oahd >/dev/null 2>&1 || [ -f /Library/Apple/usr/libexec/oah/libRosettaRuntime ]; then
    ok "Rosetta 2 present — needed to run the Intel 4.1 binary on Apple Silicon"
  else
    warn "Rosetta 2 not detected — NetMHCpan-4.1 (x86_64) needs it: softwareupdate --install-rosetta"
  fi
fi

# $1 = label, $2 = command
check_version () {
  local label="$1" cmd="$2"
  echo ""
  echo "${label}  (command: ${cmd})"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    warn "not on PATH — skipping (set ${3:-CMD} or add its launcher symlink to ~/.local/bin)"
    return 0
  fi
  ok "resolves to $(command -v "$cmd")"
  # live EL prediction
  if "$cmd" -p "$TMP/p.pep" -a HLA-A02:01 > "$TMP/${label}.out" 2>"$TMP/${label}.err"; then
    local n; n=$(grep -cE "HLA-A\*?02:01.*(NLVPMVATV|GILGFVFTL)" "$TMP/${label}.out")
    if [ "$n" -ge 2 ]; then ok "live prediction OK ($n rows; SB calls expected for these peptides)"
    else bad "ran but produced $n prediction rows (expected 2)"; rc=1; fi
  else
    bad "prediction failed (exit $?)"; sed 's/^/      /' "$TMP/${label}.err" | head -4; rc=1
  fi
  # xls path (this is where the 4.1 gawk dependency bites)
  if "$cmd" -p "$TMP/p.pep" -a HLA-A02:01 -xls -xlsfile "$TMP/${label}.xls" >/dev/null 2>&1 && [ -s "$TMP/${label}.xls" ]; then
    ok "-xls output OK ($(wc -c < "$TMP/${label}.xls" | tr -d ' ') bytes)"
  else
    bad "-xls produced no/empty file (4.1: install gawk; check -xlsfile path is writable)"; rc=1
  fi
}

check_version "NetMHCpan-4.2c" "$CMD42" "NETMHCPAN_42_CMD"
check_version "NetMHCpan-4.1b" "$CMD41" "NETMHCPAN_41_CMD"

echo ""
if [ "$rc" -eq 0 ]; then echo "${GREEN}All present NetMHCpan installs are working.${RST}"
else echo "${RED}One or more present installs are broken — see ✗ above.${RST}"; fi
echo "${DIM}Reminder: 4.1 and 4.2 EL scores differ (EL net was retrained) — don't mix versions in one analysis.${RST}"
exit $rc
