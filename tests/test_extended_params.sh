#!/usr/bin/env bash
# test_extended_params.sh
# Reproducibility tests covering parameter combinations not in existing suites:
#   - UMI length variation (-l 8/12/15/20/25): tests trimming path
#   - Min-founder threshold (-f 2/3/5): shifts low-reads boundary
#   - Explicit distance (-m 1/2/3) on sim3 with above-default pool/prod
#   - Greedy mode (-g) on sim2 and sim3 with above-default pool/prod
#   - sim3 (35k UMIs) with above-default pool/prod across error rates
#   - sim1 (real data) with above-default pool/prod
#
# Each config: reference at t=8 p=50000 d=2000, then sweep:
#   threads (2/8/16/32/48), pool (50000/100000/200000),
#   prod (2000/5000/10000), all-serial (prod>=n_umi), combined combos.
#
# Usage:  bash tests/test_extended_params.sh [--quick]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/UMI-nea/UMI-nea"
TESTS_DIR="$REPO_ROOT/tests"

PASS=0; FAIL=0; declare -a ERRORS=()

die()          { echo "FATAL: $*" >&2; exit 1; }
sorted_output(){ sort -k2,2 -k3,3 "$1"; }

run_one() {
    local label="$1"; shift
    local out="$TMPDIR/${label}.clustered"
    "$BINARY" "$@" -o "$out" >/dev/null 2>&1
    sorted_output "$out"
}

compare() {
    local label="$1" result="$2" ref="$3"
    if diff -q "$ref" "$result" >/dev/null 2>&1; then
        echo "  PASS  $label"; PASS=$(( PASS+1 ))
    else
        echo "  FAIL  $label"; FAIL=$(( FAIL+1 )); ERRORS+=("$label")
        diff "$ref" "$result" | head -20 || true
    fi
}

# sweep_exec <base_label> <ref_sorted> <input> <n_umi> <algo_args...>
# Sweeps t/p/d including above-default values, comparing to pre-built ref.
sweep_exec() {
    local base_label="$1" ref="$2" input="$3" n_umi="$4"; shift 4
    local algo_args=("$@")

    # Thread sweep (pool=50000, prod=2000 — both above default)
    for t in "${THREAD_LIST[@]}"; do
        local r="$TMPDIR/${base_label}_t${t}.sorted"
        run_one "${base_label}_t${t}" "${algo_args[@]}" -t "$t" -p 50000 -d 2000 -i "$input" > "$r"
        compare "t=$t  p=50000  d=2000  [$base_label]" "$r" "$ref"
    done

    # Pool sweep (t=8, prod=2000)
    for p in "${POOL_LIST[@]}"; do
        local r="$TMPDIR/${base_label}_p${p}.sorted"
        run_one "${base_label}_p${p}" "${algo_args[@]}" -t 8 -p "$p" -d 2000 -i "$input" > "$r"
        compare "t=8  p=$p  d=2000  [$base_label]" "$r" "$ref"
    done

    # Prod sweep (t=8, pool=50000)
    for d in "${PROD_LIST[@]}"; do
        local r="$TMPDIR/${base_label}_d${d}.sorted"
        run_one "${base_label}_d${d}" "${algo_args[@]}" -t 8 -p 50000 -d "$d" -i "$input" > "$r"
        compare "t=8  p=50000  d=$d  [$base_label]" "$r" "$ref"
    done

    # All-serial path: prod >= n_umi (no consumers fire)
    for d in "$n_umi" $(( n_umi + 1000 )) $(( n_umi * 2 )); do
        local r="$TMPDIR/${base_label}_allserial_d${d}.sorted"
        run_one "${base_label}_allserial_d${d}" "${algo_args[@]}" -t 8 -p 200000 -d "$d" -i "$input" > "$r"
        compare "t=8  p=200000  d=$d (all-serial)  [$base_label]" "$r" "$ref"
    done

    # Combined: varied t + large p + large d
    for combo in "${COMBO_LIST[@]}"; do
        read -r t p d <<< "$combo"
        local r="$TMPDIR/${base_label}_t${t}_p${p}_d${d}.sorted"
        run_one "${base_label}_t${t}_p${p}_d${d}" "${algo_args[@]}" -t "$t" -p "$p" -d "$d" -i "$input" > "$r"
        compare "t=$t  p=$p  d=$d  [$base_label]" "$r" "$ref"
    done
}

# run_config <base_label> <input> <algo_args...>
run_config() {
    local base_label="$1" input="$2"; shift 2
    local algo_args=("$@")
    local n_umi
    n_umi=$(wc -l < "$input")
    echo ""
    echo "  Config: ${algo_args[*]}  |  $(basename "$input")  ($n_umi UMIs)"
    local ref="$TMPDIR/${base_label}_ref.sorted"
    run_one "${base_label}_ref" "${algo_args[@]}" -t 8 -p 50000 -d 2000 -i "$input" > "$ref"
    echo "  Reference: $(wc -l < "$ref") output lines."
    sweep_exec "$base_label" "$ref" "$input" "$n_umi" "${algo_args[@]}"
}

# ── setup ─────────────────────────────────────────────────────────────────────

[[ -x "$BINARY" ]] || die "binary not found: $BINARY"

SIM1="$REPO_ROOT/datafiles/UMI-nea_testfiles/sim1.input"
SIM2="$TESTS_DIR/sim2.input"
SIM3="$TESTS_DIR/sim3.input"

if [[ ! -f "$SIM2" ]] || [[ ! -f "$SIM3" ]]; then
    echo "Generating synthetic datasets..."
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

if $QUICK; then
    THREAD_LIST=(2 8 16 32)
    POOL_LIST=(50000 100000)
    PROD_LIST=(2000 5000)
    COMBO_LIST=("16 100000 5000" "32 200000 10000" "2 200000 50000")
else
    THREAD_LIST=(2 4 8 16 32 48)
    POOL_LIST=(50000 100000 200000)
    PROD_LIST=(2000 5000 10000)
    COMBO_LIST=(
        "4  50000  5000"
        "8  100000 10000"
        "16 200000 5000"
        "32 50000  10000"
        "48 100000 2000"
        "2  200000 50000"
        "16 50000  50000"   # prod > pool → all serial
    )
fi

echo "=== UMI-nea extended parameter reproducibility test  [$(date)] ==="
$QUICK && echo "(quick mode)"

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 1: UMI length sweep  (sim2, -l 8/12/15/20/25)"
echo "Tests trimming path: input UMIs are 18-mers; lengths <18 trim them."
echo "════════════════════════════════════════════════════════════"
for l in 8 12 15 20 25; do
    run_config "sim2_l${l}" "$SIM2" -l "$l" -e 0.01 -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 2: Min-founder threshold  (-f 2/3/5) on sim2 and sim3"
echo "Shifts which UMIs go to low-reads pool vs producer."
echo "════════════════════════════════════════════════════════════"
for f in 2 3 5; do
    run_config "sim2_f${f}"  "$SIM2" -l 18 -e 0.01 -f "$f"
    run_config "sim3_f${f}"  "$SIM3" -l 18 -e 0.01 -f "$f"
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 3: Explicit distance  (-m 1/2/3) on sim2 and sim3"
echo "Above-default pool/prod stress-test for each distance setting."
echo "════════════════════════════════════════════════════════════"
for m in 1 2 3; do
    run_config "sim2_m${m}" "$SIM2" -l 18 -m "$m" -f 1
    run_config "sim3_m${m}" "$SIM3" -l 18 -m "$m" -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 4: Greedy mode  (-g) on sim2 and sim3"
echo "First-found founder rather than nearest: tests different code path."
echo "════════════════════════════════════════════════════════════"
run_config "sim2_greedy" "$SIM2" -l 18 -e 0.01 -f 1 -g
run_config "sim3_greedy" "$SIM3" -l 18 -e 0.01 -f 1 -g

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 5: sim3 (35k UMIs) — error rate sweep with above-default pool/prod"
echo "Large dataset: tests multi-batch boundary behaviour at scale."
echo "════════════════════════════════════════════════════════════"
for e in 0.001 0.005 0.01 0.02 0.05; do
    tag="e${e//./p}"
    run_config "sim3_${tag}" "$SIM3" -l 18 -e "$e" -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 6: sim1 (real data, 4165 UMIs) — above-default pool/prod"
echo "Validates on real sequencing data with large execution params."
echo "════════════════════════════════════════════════════════════"
run_config "sim1_e0p01"  "$SIM1" -l 18 -e 0.01  -f 1
run_config "sim1_e0p001" "$SIM1" -l 18 -e 0.001 -f 1
run_config "sim1_m2"     "$SIM1" -l 18 -m 2     -f 1
run_config "sim1_f2"     "$SIM1" -l 18 -e 0.01  -f 2
run_config "sim1_greedy" "$SIM1" -l 18 -e 0.01  -f 1 -g

# ══════════════════════════════════════════════════════════════════════════════
TOTAL=$(( PASS+FAIL ))
echo ""
echo "=== Results: $PASS/$TOTAL passed  [$(date)] ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED configurations:"
    for e in "${ERRORS[@]}"; do echo "  - $e"; done
    exit 1
fi
echo "All extended parameter reproducibility checks PASSED."
