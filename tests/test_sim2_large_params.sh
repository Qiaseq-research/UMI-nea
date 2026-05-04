#!/usr/bin/env bash
# test_sim2_large_params.sh
# Test reproducibility on sim2.input using pool_size > 48000 (default) and
# prod_size > 1000 (default), swept across threads and error rates.
#
# sim2.input has 13078 UMIs, so every pool_size here is a single batch.
# The interesting variation is how much of that batch goes to the serial
# producer (prod_size) vs parallel consumers, and at what thread count.

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

echo "=== sim2.input large-param sweep  ($N_UMI UMIs)  [$(date)] ==="
echo "Defaults: pool_size=48000  prod_size=1000"
echo "All pool_size and prod_size values below EXCEED the defaults."
echo ""

# Error rates → max_dist for 18-mer UMIs:
#   0.001→1  0.005→1  0.01→2  0.02→3  0.05→4
ERROR_RATES=(0.001 0.005 0.01 0.02 0.05)

THREADS=(2 4 6 8 12 16 24 32 48 64)

# pool_size values — all > 48000 (default)
POOLS=(50000 60000 80000 100000 150000 200000)

# prod_size values — all > 1000 (default)
PRODS=(2000 3000 5000 8000 10000 15000 20000)

# prod close to / equal to / exceeding n_umi (→ all-serial path)
# n_umi = 13078, so prod >= 13078 means no consumers fire
PROD_ALLSERIAL=(13000 13078 15000 20000 50000)

for e in "${ERROR_RATES[@]}"; do
    echo "════════════════════════════════════════════════════════"
    echo "error_rate = $e   (reference: t=8  pool=50000  prod=2000)"
    echo "════════════════════════════════════════════════════════"

    REF="$TMPDIR/e${e}_ref.sorted"
    run_one "e${e}_ref" -l 18 -e "$e" -t 8 -p 50000 -d 2000 -f 1 > "$REF"
    echo "Reference: $(wc -l < "$REF") output lines."

    # ── thread sweep  (pool=50000, prod=2000) ────────────────────────
    echo "  [thread sweep  pool=50000  prod=2000]"
    for t in "${THREADS[@]}"; do
        r="$TMPDIR/e${e}_t${t}_p50000_d2000.sorted"
        run_one "e${e}_t${t}_p50000_d2000" \
            -l 18 -e "$e" -t "$t" -p 50000 -d 2000 -f 1 > "$r"
        check "e=$e  t=$t  pool=50000  prod=2000" "$r"
    done

    # ── pool_size sweep  (t=8, prod=2000) ────────────────────────────
    echo "  [pool sweep  t=8  prod=2000]"
    for p in "${POOLS[@]}"; do
        r="$TMPDIR/e${e}_t8_p${p}_d2000.sorted"
        run_one "e${e}_t8_p${p}_d2000" \
            -l 18 -e "$e" -t 8 -p "$p" -d 2000 -f 1 > "$r"
        check "e=$e  t=8  pool=$p  prod=2000" "$r"
    done

    # ── prod_size sweep  (t=8, pool=50000) ───────────────────────────
    echo "  [prod sweep  t=8  pool=50000]"
    for d in "${PRODS[@]}"; do
        r="$TMPDIR/e${e}_t8_p50000_d${d}.sorted"
        run_one "e${e}_t8_p50000_d${d}" \
            -l 18 -e "$e" -t 8 -p 50000 -d "$d" -f 1 > "$r"
        check "e=$e  t=8  pool=50000  prod=$d" "$r"
    done

    # ── prod >= n_umi: all-serial path (no consumers) ─────────────────
    echo "  [all-serial path: prod >= n_umi, pool=100000, t=8]"
    for d in "${PROD_ALLSERIAL[@]}"; do
        r="$TMPDIR/e${e}_t8_p100000_d${d}.sorted"
        run_one "e${e}_t8_p100000_d${d}" \
            -l 18 -e "$e" -t 8 -p 100000 -d "$d" -f 1 > "$r"
        check "e=$e  t=8  pool=100000  prod=$d (all-serial)" "$r"
    done

    # ── combined large-param sweep ────────────────────────────────────
    echo "  [combined large-param sweep]"
    # Each row: threads pool prod
    COMBOS=(
        "2   50000  2000"
        "4   60000  3000"
        "8   80000  5000"
        "12  100000 8000"
        "16  150000 10000"
        "24  200000 15000"
        "32  50000  20000"
        "48  100000 2000"
        "64  200000 5000"
        "2   100000 13000"   # prod > n_umi → all serial, low thread count
        "16  80000  15000"   # large prod, mid threads
        "32  150000 10000"   # high threads, very large pool
        "64  50000  13078"   # prod exactly = n_umi
        "48  200000 20000"   # extreme pool + prod, high threads
        "8   50000  50000"   # prod > pool (prod clamped to pool → all serial)
    )
    for combo in "${COMBOS[@]}"; do
        read -r t p d <<< "$combo"
        label="e${e}_t${t}_p${p}_d${d}"
        r="$TMPDIR/${label}.sorted"
        run_one "$label" -l 18 -e "$e" -t "$t" -p "$p" -d "$d" -f 1 > "$r"
        check "e=$e  t=$t  pool=$p  prod=$d" "$r"
    done

    echo ""
done

TOTAL=$(( PASS+FAIL ))
echo "=== Results: $PASS/$TOTAL passed  [$(date)] ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED:"
    for err in "${ERRORS[@]}"; do echo "  - $err"; done
    exit 1
fi
echo "All checks PASSED."
