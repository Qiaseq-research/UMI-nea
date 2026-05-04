#!/usr/bin/env bash
# test_reproducibility.sh
# Verify that UMI-nea produces identical clustering results across different
# -t (threads), -p (pool_size), and -d (prod_size) settings.
#
# Strategy: run a reference configuration, then compare every other config's
# sorted output against it.  Sorting by (UMI, founder) normalises the output
# ordering because threading changes the write order but must not change the
# clustering assignments.
#
# Usage:  bash tests/test_reproducibility.sh [--quick]
#   --quick  run a reduced config matrix (faster CI smoke test)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/UMI-nea/UMI-nea"
INPUT="$REPO_ROOT/datafiles/UMI-nea_testfiles/sim1.input"

PASS=0
FAIL=0
ERRORS=()

# ── helpers ──────────────────────────────────────────────────────────────────

die() { echo "FATAL: $*" >&2; exit 1; }

sorted_output() {
    # Sort by col2 (UMI), then col3 (founder) for stable comparison.
    sort -k2,2 -k3,3 "$1"
}

run_one() {
    local label="$1"; shift
    local outfile="$TMPDIR/${label}.clustered"
    "$BINARY" "$@" -i "$INPUT" -o "$outfile" -f 1 > /dev/null 2>&1
    sorted_output "$outfile"
}

compare() {
    local label="$1"
    local result="$2"      # file path of sorted result
    if diff -q "$REFERENCE_SORTED" "$result" > /dev/null 2>&1; then
        echo "  PASS  $label"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $label"
        FAIL=$(( FAIL + 1 ))
        ERRORS+=("$label")
        # Show first differing lines for diagnosis
        diff "$REFERENCE_SORTED" "$result" | head -20 || true
    fi
}

# ── setup ────────────────────────────────────────────────────────────────────

[[ -x "$BINARY" ]] || die "binary not found: $BINARY  (run 'make' in UMI-nea/)"
[[ -f "$INPUT"  ]] || die "input not found: $INPUT"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

QUICK=false
[[ "${1:-}" == "--quick" ]] && QUICK=true

# ── reference run ────────────────────────────────────────────────────────────
# t=4, pool=4000 (just under total 4165 → one multi-batch pass), prod=200
echo "=== UMI-nea reproducibility test ==="
echo "Input: $INPUT  ($(wc -l < "$INPUT") UMIs)"
echo ""
echo "Building reference  t=4  pool=4000  prod=200 ..."
REFERENCE_SORTED="$TMPDIR/reference.sorted"
run_one "reference" -l 18 -e 0.01 -t 4 -p 4000 -d 200 > "$REFERENCE_SORTED"
echo "Reference has $(wc -l < "$REFERENCE_SORTED") output lines."
echo ""

# ── thread-count sweep ───────────────────────────────────────────────────────
echo "--- Thread count sweep (pool=4000, prod=200) ---"
for t in 2 6 8 12 16 24 32 48; do
    $QUICK && [[ "$t" -gt 8 ]] && continue
    result="$TMPDIR/t${t}_p4000_d200.sorted"
    run_one "t${t}_p4000_d200" -l 18 -e 0.01 -t "$t" -p 4000 -d 200 > "$result"
    compare "threads=$t  pool=4000  prod=200" "$result"
done
echo ""

# ── pool_size sweep: single batch (pool > total UMIs) ───────────────────────
echo "--- Pool size sweep: single batch (pool >= 4165, all UMIs in one batch) ---"
for p in 4165 5000 8000 16000 48000; do
    $QUICK && [[ "$p" -gt 8000 ]] && continue
    result="$TMPDIR/t4_p${p}_d200.sorted"
    run_one "t4_p${p}_d200" -l 18 -e 0.01 -t 4 -p "$p" -d 200 > "$result"
    compare "threads=4  pool=$p  prod=200" "$result"
done
echo ""

# ── pool_size sweep: multi-batch ─────────────────────────────────────────────
echo "--- Pool size sweep: multi-batch (pool < 4165) ---"
for p in 200 500 1000 2000 3000; do
    $QUICK && [[ "$p" -lt 500 ]] && continue
    result="$TMPDIR/t4_p${p}_d200.sorted"
    run_one "t4_p${p}_d200" -l 18 -e 0.01 -t 4 -p "$p" -d 200 > "$result"
    compare "threads=4  pool=$p  prod=200" "$result"
done
echo ""

# ── prod_size sweep ───────────────────────────────────────────────────────────
echo "--- Prod size sweep (threads=4, pool=4000) ---"
for d in 50 100 250 500 1000 2000; do
    $QUICK && { [[ "$d" -lt 100 ]] || [[ "$d" -gt 1000 ]]; } && continue
    result="$TMPDIR/t4_p4000_d${d}.sorted"
    run_one "t4_p4000_d${d}" -l 18 -e 0.01 -t 4 -p 4000 -d "$d" > "$result"
    compare "threads=4  pool=4000  prod=$d" "$result"
done
echo ""

# ── combined sweep ────────────────────────────────────────────────────────────
echo "--- Combined parameter sweep ---"
# Each row: threads pool prod
COMBOS=(
    "2  500   50"
    "4  1000  100"
    "8  2000  250"
    "16 3000  500"
    "32 500   1000"
    "48 2000  2000"
    "8  200   50"
    "16 4000  2000"
)
for combo in "${COMBOS[@]}"; do
    read -r t p d <<< "$combo"
    $QUICK && [[ "$t" -gt 8 ]] && continue
    label="t${t}_p${p}_d${d}"
    result="$TMPDIR/${label}.sorted"
    run_one "$label" -l 18 -e 0.01 -t "$t" -p "$p" -d "$d" > "$result"
    compare "threads=$t  pool=$p  prod=$d" "$result"
done
echo ""

# ── same-config repeated runs (non-determinism within a single config) ────────
echo "--- Same-config repeated runs (t=8, pool=2000, prod=300) ---"
REPEAT_REF="$TMPDIR/repeat_ref.sorted"
run_one "repeat_ref" -l 18 -e 0.01 -t 8 -p 2000 -d 300 > "$REPEAT_REF"
for run in 2 3 4 5; do
    result="$TMPDIR/repeat_run${run}.sorted"
    run_one "repeat_run${run}" -l 18 -e 0.01 -t 8 -p 2000 -d 300 > "$result"
    if diff -q "$REPEAT_REF" "$result" > /dev/null 2>&1; then
        echo "  PASS  same-config repeated run $run"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  same-config repeated run $run"
        FAIL=$(( FAIL + 1 ))
        ERRORS+=("same-config run $run")
        diff "$REPEAT_REF" "$result" | head -20 || true
    fi
done
echo ""

# ── max-dist / -m flag variants ───────────────────────────────────────────────
echo "--- Cross-thread consistency with -m (explicit distance) ---"
for t in 2 8 16; do
    $QUICK && [[ "$t" -gt 8 ]] && continue
    result="$TMPDIR/mdist_t${t}.sorted"
    run_one "mdist_t${t}" -l 18 -m 2 -t "$t" -p 4000 -d 200 > "$result"
done
# compare t=2 vs t=8 and (full mode) t=16
diff -q "$TMPDIR/mdist_t2.sorted" "$TMPDIR/mdist_t8.sorted" > /dev/null 2>&1 \
    && { echo "  PASS  -m 2: threads=2 vs threads=8"; PASS=$(( PASS+1 )); } \
    || { echo "  FAIL  -m 2: threads=2 vs threads=8"; FAIL=$(( FAIL+1 )); ERRORS+=("-m2 t2 vs t8"); }
if ! $QUICK; then
    diff -q "$TMPDIR/mdist_t2.sorted" "$TMPDIR/mdist_t16.sorted" > /dev/null 2>&1 \
        && { echo "  PASS  -m 2: threads=2 vs threads=16"; PASS=$(( PASS+1 )); } \
        || { echo "  FAIL  -m 2: threads=2 vs threads=16"; FAIL=$(( FAIL+1 )); ERRORS+=("-m2 t2 vs t16"); }
fi
echo ""

# ── summary ───────────────────────────────────────────────────────────────────
TOTAL=$(( PASS + FAIL ))
echo "=== Results: $PASS/$TOTAL passed ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED configurations:"
    for e in "${ERRORS[@]}"; do
        echo "  - $e"
    done
    exit 1
fi
echo "All reproducibility checks PASSED."
