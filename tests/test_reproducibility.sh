#!/usr/bin/env bash
# test_reproducibility.sh
# Verify that UMI-nea produces identical clustering results across different
# -t (threads), -p (pool_size), and -d (prod_size) settings.
#
# Strategy: for each dataset, run a reference configuration, then compare
# every other config's sorted output against it.  Sorting by (UMI, founder)
# normalises write-order differences — threading must not change assignments.
#
# Datasets:
#   sim1  – real data, 4 165 UMIs
#   sim2  – synthetic, ~13 k UMIs (1 500 founders)
#   sim3  – synthetic, ~35 k UMIs (4 000 founders)
#
# Usage:
#   bash tests/test_reproducibility.sh           # full suite
#   bash tests/test_reproducibility.sh --quick   # reduced matrix, CI smoke test

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/UMI-nea/UMI-nea"
TESTS_DIR="$REPO_ROOT/tests"

PASS=0
FAIL=0
ERRORS=()

# ── helpers ──────────────────────────────────────────────────────────────────

die() { echo "FATAL: $*" >&2; exit 1; }

sorted_output() { sort -k2,2 -k3,3 "$1"; }

# run_one <label> [UMI-nea args...]  → prints sorted clustered output to stdout
run_one() {
    local label="$1"; shift
    local outfile="$TMPDIR/${label}.clustered"
    "$BINARY" "$@" -o "$outfile" > /dev/null 2>&1
    sorted_output "$outfile"
}

compare() {
    local label="$1" result="$2"
    if diff -q "$REFERENCE_SORTED" "$result" > /dev/null 2>&1; then
        echo "  PASS  $label"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $label"
        FAIL=$(( FAIL + 1 ))
        ERRORS+=("$label")
        diff "$REFERENCE_SORTED" "$result" | head -30 || true
    fi
}

section() { echo ""; echo "--- $* ---"; }

# run_sweep <input> <ref_t> <ref_p> <ref_d>
# Runs a full parameter sweep for one input file against its own reference.
run_sweep() {
    local input="$1" ref_t="$2" ref_p="$3" ref_d="$4"
    local n_umi base_label
    n_umi=$(wc -l < "$input")
    base_label="$(basename "$input" .input)"
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "Dataset: $input  ($n_umi UMIs)"
    echo "Reference: t=$ref_t  pool=$ref_p  prod=$ref_d"
    echo "════════════════════════════════════════════════════════════"

    REFERENCE_SORTED="$TMPDIR/${base_label}_ref.sorted"
    run_one "${base_label}_ref" -l 18 -e 0.01 -t "$ref_t" -p "$ref_p" -d "$ref_d" -f 1 -i "$input" \
        > "$REFERENCE_SORTED"
    echo "Reference has $(wc -l < "$REFERENCE_SORTED") output lines."

    # ── thread-count sweep ───────────────────────────────────────────────
    section "Thread sweep  pool=$ref_p  prod=$ref_d"
    for t in "${THREAD_LIST[@]}"; do
        local result="$TMPDIR/${base_label}_t${t}_p${ref_p}_d${ref_d}.sorted"
        run_one "${base_label}_t${t}_p${ref_p}_d${ref_d}" \
            -l 18 -e 0.01 -t "$t" -p "$ref_p" -d "$ref_d" -f 1 -i "$input" > "$result"
        compare "t=$t  pool=$ref_p  prod=$ref_d  ($base_label)" "$result"
    done

    # ── pool_size: single batch (pool ≥ total UMIs) ──────────────────────
    section "Pool ≥ n_umi (single batch)  t=$ref_t  prod=$ref_d"
    for p in "${SINGLE_BATCH_POOLS[@]}"; do
        [[ "$p" -lt "$n_umi" ]] && continue   # skip if not actually single-batch
        local result="$TMPDIR/${base_label}_t${ref_t}_p${p}_d${ref_d}.sorted"
        run_one "${base_label}_t${ref_t}_p${p}_d${ref_d}" \
            -l 18 -e 0.01 -t "$ref_t" -p "$p" -d "$ref_d" -f 1 -i "$input" > "$result"
        compare "t=$ref_t  pool=$p  prod=$ref_d  ($base_label)" "$result"
    done

    # ── pool_size: multi-batch (pool < total UMIs) ───────────────────────
    section "Pool < n_umi (multi-batch)  t=$ref_t  prod=$ref_d"
    for p in "${MULTI_BATCH_POOLS[@]}"; do
        [[ "$p" -ge "$n_umi" ]] && continue   # skip if not actually multi-batch
        local result="$TMPDIR/${base_label}_t${ref_t}_p${p}_d${ref_d}.sorted"
        run_one "${base_label}_t${ref_t}_p${p}_d${ref_d}" \
            -l 18 -e 0.01 -t "$ref_t" -p "$p" -d "$ref_d" -f 1 -i "$input" > "$result"
        compare "t=$ref_t  pool=$p  prod=$ref_d  ($base_label)" "$result"
    done

    # ── prod_size sweep ──────────────────────────────────────────────────
    section "Prod sweep  t=$ref_t  pool=$ref_p"
    for d in "${PROD_LIST[@]}"; do
        local result="$TMPDIR/${base_label}_t${ref_t}_p${ref_p}_d${d}.sorted"
        run_one "${base_label}_t${ref_t}_p${ref_p}_d${d}" \
            -l 18 -e 0.01 -t "$ref_t" -p "$ref_p" -d "$d" -f 1 -i "$input" > "$result"
        compare "t=$ref_t  pool=$ref_p  prod=$d  ($base_label)" "$result"
    done

    # ── combined sweep ───────────────────────────────────────────────────
    section "Combined sweep"
    for combo in "${COMBOS[@]}"; do
        read -r t p d <<< "$combo"
        [[ "$p" -gt "$n_umi" ]] && continue   # skip absurdly large pools for small datasets
        local label="${base_label}_t${t}_p${p}_d${d}"
        local result="$TMPDIR/${label}.sorted"
        run_one "$label" -l 18 -e 0.01 -t "$t" -p "$p" -d "$d" -f 1 -i "$input" > "$result"
        compare "t=$t  pool=$p  prod=$d  ($base_label)" "$result"
    done

    # ── prod ≥ pool edge case (everything serial) ────────────────────────
    section "Prod ≥ pool edge case  t=$ref_t"
    for p in 100 500 1000; do
        [[ "$p" -gt "$n_umi" ]] && continue
        local d="$p"
        local result="$TMPDIR/${base_label}_t${ref_t}_p${p}_d${d}_edge.sorted"
        run_one "${base_label}_t${ref_t}_p${p}_d${d}_edge" \
            -l 18 -e 0.01 -t "$ref_t" -p "$p" -d "$d" -f 1 -i "$input" > "$result"
        compare "t=$ref_t  pool=$p  prod=$d (prod=pool edge)  ($base_label)" "$result"
    done

    # ── same-config repeated runs (within-config non-determinism) ────────
    section "Same-config repeated runs  t=8  pool=1000  prod=200"
    local rep_ref="$TMPDIR/${base_label}_repeat_ref.sorted"
    run_one "${base_label}_repeat_ref" \
        -l 18 -e 0.01 -t 8 -p 1000 -d 200 -f 1 -i "$input" > "$rep_ref"
    for run in 2 3 4 5; do
        local result="$TMPDIR/${base_label}_repeat_run${run}.sorted"
        run_one "${base_label}_repeat_run${run}" \
            -l 18 -e 0.01 -t 8 -p 1000 -d 200 -f 1 -i "$input" > "$result"
        if diff -q "$rep_ref" "$result" > /dev/null 2>&1; then
            echo "  PASS  same-config run $run  ($base_label)"
            PASS=$(( PASS + 1 ))
        else
            echo "  FAIL  same-config run $run  ($base_label)"
            FAIL=$(( FAIL + 1 ))
            ERRORS+=("same-config run $run ($base_label)")
            diff "$rep_ref" "$result" | head -20 || true
        fi
    done
}

# ── setup ────────────────────────────────────────────────────────────────────

[[ -x "$BINARY" ]] || die "binary not found: $BINARY  (run 'make' in UMI-nea/)"

SIM1="$REPO_ROOT/datafiles/UMI-nea_testfiles/sim1.input"
SIM2="$TESTS_DIR/sim2.input"
SIM3="$TESTS_DIR/sim3.input"

# Generate synthetic datasets if missing
if [[ ! -f "$SIM2" ]] || [[ ! -f "$SIM3" ]]; then
    echo "Generating synthetic test datasets..."
    python3 "$TESTS_DIR/generate_test_input.py" 1500 18 42 "$SIM2"
    python3 "$TESTS_DIR/generate_test_input.py" 4000 18 99 "$SIM3"
fi

for f in "$SIM1" "$SIM2" "$SIM3"; do
    [[ -f "$f" ]] || die "input not found: $f"
done

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

QUICK=false
[[ "${1:-}" == "--quick" ]] && QUICK=true

# ── parameter matrices ────────────────────────────────────────────────────────

if $QUICK; then
    THREAD_LIST=(2 4 8 16)
    SINGLE_BATCH_POOLS=(50000 100000)
    MULTI_BATCH_POOLS=(100 500 2000)
    PROD_LIST=(50 200 1000)
    COMBOS=(
        "2   500   50"
        "4   1000  100"
        "8   2000  250"
        "16  5000  500"
    )
else
    THREAD_LIST=(2 4 6 8 12 16 24 32 48 64)
    SINGLE_BATCH_POOLS=(50000 100000)
    MULTI_BATCH_POOLS=(50 100 200 500 1000 2000 5000 10000)
    PROD_LIST=(25 50 100 200 500 1000 2000)
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
        "32  20000 2000"
        "48  2000  200"
        "48  50000 2000"
        "64  10000 500"
        "64  1000  100"
    )
fi

echo "=== UMI-nea reproducibility test  [$(date)] ==="
$QUICK && echo "(quick mode)"

# ── run sweeps ───────────────────────────────────────────────────────────────

# sim1: real data 4165 UMIs  — reference at t=4, pool=4000, prod=200
run_sweep "$SIM1" 4 4000 200

# sim2: synthetic ~13k UMIs  — reference at t=8, pool=10000, prod=500
run_sweep "$SIM2" 8 10000 500

# sim3: synthetic ~35k UMIs  — reference at t=8, pool=30000, prod=1000
run_sweep "$SIM3" 8 30000 1000

# ── cross-dataset: -m flag (explicit distance) ────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Cross-thread check with explicit -m 2"
echo "════════════════════════════════════════════════════════════"
for input_f in "$SIM1" "$SIM2" "$SIM3"; do
    base="$(basename "$input_f" .input)"
    n_umi=$(wc -l < "$input_f")
    ref_m="$TMPDIR/${base}_m2_ref.sorted"
    run_one "${base}_m2_t4" -l 18 -m 2 -t 4 -p 4000 -d 200 -f 1 -i "$input_f" > "$ref_m"
    for t in 8 16 32; do
        $QUICK && [[ "$t" -gt 16 ]] && continue
        result="$TMPDIR/${base}_m2_t${t}.sorted"
        run_one "${base}_m2_t${t}" -l 18 -m 2 -t "$t" -p 4000 -d 200 -f 1 -i "$input_f" > "$result"
        if diff -q "$ref_m" "$result" > /dev/null 2>&1; then
            echo "  PASS  -m 2  t=4 vs t=$t  ($base)"
            PASS=$(( PASS + 1 ))
        else
            echo "  FAIL  -m 2  t=4 vs t=$t  ($base)"
            FAIL=$(( FAIL + 1 ))
            ERRORS+=("-m2 t4 vs t$t ($base)")
            diff "$ref_m" "$result" | head -20 || true
        fi
    done
done

# ── summary ───────────────────────────────────────────────────────────────────

TOTAL=$(( PASS + FAIL ))
echo ""
echo "=== Results: $PASS/$TOTAL passed  [$(date)] ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED configurations:"
    for e in "${ERRORS[@]}"; do
        echo "  - $e"
    done
    exit 1
fi
echo "All reproducibility checks PASSED."
