#!/usr/bin/env bash
# test_sim2_full_sweep.sh
# For each error rate, establish a reference run on sim2.input, then sweep
# threads x pool_size x prod_size and verify all outputs are identical.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/UMI-nea/UMI-nea"
INPUT="$REPO_ROOT/tests/sim2.input"
N_UMI=$(wc -l < "$INPUT")

[[ -x "$BINARY" ]] || { echo "FATAL: binary not found: $BINARY"; exit 1; }
[[ -f "$INPUT"  ]] || { echo "FATAL: input not found: $INPUT"; exit 1; }

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

PASS=0; FAIL=0; declare -a ERRORS=()

run_one() {
    local tag="$1"; shift
    local out="$TMPDIR/${tag}.clustered"
    "$BINARY" "$@" -i "$INPUT" -o "$out" > /dev/null 2>&1
    sort -k2,2 -k3,3 "$out"
}

check() {
    local label="$1" result="$2"
    if diff -q "$REF" "$result" > /dev/null 2>&1; then
        echo "  PASS  $label"; PASS=$(( PASS+1 ))
    else
        echo "  FAIL  $label"; FAIL=$(( FAIL+1 )); ERRORS+=("$label")
        diff "$REF" "$result" | head -20 || true
    fi
}

echo "=== sim2.input full sweep  ($N_UMI UMIs)  [$(date)] ==="
echo "Binary: $BINARY"
echo ""

# Error rates to test.  Each implies a different max_dist for 18-mer UMIs:
#   e=0.001 → dist=1   e=0.005 → dist=1   e=0.01 → dist=2
#   e=0.02  → dist=3   e=0.05  → dist=4
ERROR_RATES=(0.001 0.005 0.01 0.02 0.05)

# Sweeps (kept the same for every error rate so results are comparable)
THREADS=(2 4 6 8 12 16 24 32 48 64)

# Pool sizes split into multi-batch (< N_UMI) and single-batch (>= N_UMI)
MULTI_POOLS=(50 100 200 500 1000 2000 5000 10000)
SINGLE_POOLS=(15000 50000 100000)

PROD_SIZES=(25 50 100 200 500 1000 2000)

# Combined (t, pool, prod) combos — exercise the cross-parameter space
COMBOS=(
    "2   100   25"
    "2   500   50"
    "4   1000  100"
    "4   200   200"
    "8   2000  250"
    "8   500   500"
    "12  3000  300"
    "16  5000  500"
    "16  1000  1000"
    "24  10000 1000"
    "32  500   50"
    "32  12000 2000"
    "48  2000  200"
    "64  10000 500"
    "64  1000  100"
)

for e in "${ERROR_RATES[@]}"; do
    echo "════════════════════════════════════════════════════════"
    echo "error_rate = $e  (reference: t=8  pool=10000  prod=500)"
    echo "════════════════════════════════════════════════════════"

    REF="$TMPDIR/e${e}_ref.sorted"
    run_one "e${e}_ref" -l 18 -e "$e" -t 8 -p 10000 -d 500 -f 1 > "$REF"
    echo "Reference: $(wc -l < "$REF") output lines."

    # ── thread sweep ──────────────────────────────────────────────────
    echo "  [thread sweep  pool=10000  prod=500]"
    for t in "${THREADS[@]}"; do
        r="$TMPDIR/e${e}_t${t}_p10000_d500.sorted"
        run_one "e${e}_t${t}_p10000_d500" -l 18 -e "$e" -t "$t" -p 10000 -d 500 -f 1 > "$r"
        check "e=$e  t=$t  pool=10000  prod=500" "$r"
    done

    # ── multi-batch pool sweep ────────────────────────────────────────
    echo "  [multi-batch pool sweep  t=8  prod=500]"
    for p in "${MULTI_POOLS[@]}"; do
        r="$TMPDIR/e${e}_t8_p${p}_d500.sorted"
        run_one "e${e}_t8_p${p}_d500" -l 18 -e "$e" -t 8 -p "$p" -d 500 -f 1 > "$r"
        check "e=$e  t=8  pool=$p  prod=500" "$r"
    done

    # ── single-batch pool sweep ───────────────────────────────────────
    echo "  [single-batch pool sweep  t=8  prod=500]"
    for p in "${SINGLE_POOLS[@]}"; do
        r="$TMPDIR/e${e}_t8_p${p}_d500.sorted"
        run_one "e${e}_t8_p${p}_d500" -l 18 -e "$e" -t 8 -p "$p" -d 500 -f 1 > "$r"
        check "e=$e  t=8  pool=$p  prod=500" "$r"
    done

    # ── prod_size sweep ───────────────────────────────────────────────
    echo "  [prod sweep  t=8  pool=10000]"
    for d in "${PROD_SIZES[@]}"; do
        r="$TMPDIR/e${e}_t8_p10000_d${d}.sorted"
        run_one "e${e}_t8_p10000_d${d}" -l 18 -e "$e" -t 8 -p 10000 -d "$d" -f 1 > "$r"
        check "e=$e  t=8  pool=10000  prod=$d" "$r"
    done

    # ── combined sweep ────────────────────────────────────────────────
    echo "  [combined sweep]"
    for combo in "${COMBOS[@]}"; do
        read -r t p d <<< "$combo"
        r="$TMPDIR/e${e}_t${t}_p${p}_d${d}.sorted"
        run_one "e${e}_t${t}_p${p}_d${d}" -l 18 -e "$e" -t "$t" -p "$p" -d "$d" -f 1 > "$r"
        check "e=$e  t=$t  pool=$p  prod=$d" "$r"
    done

    echo ""
done

TOTAL=$(( PASS+FAIL ))
echo "=== Results: $PASS/$TOTAL passed  [$(date)] ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED:"
    for e in "${ERRORS[@]}"; do echo "  - $e"; done
    exit 1
fi
echo "All checks PASSED."
