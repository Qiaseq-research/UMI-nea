#!/usr/bin/env bash
# test_extended_params_sim2.sh
# Same reproducibility test as test_extended_params.sh but run on the large
# benchmark dataset: /home/jdeng/claude/benchmark/sim2.input (519k 24-mer UMIs).
#
# NOTE: Each run takes ~6 min on this dataset. To keep runtime feasible:
#   - All-serial paths are skipped (prod=519k would serialize 519k UMIs)
#   - Pool/prod sweeps are skipped in default mode; use --full to enable them
#   - Default thread list: 2 8 16 32; --quick mode: 2 16
#
# Usage:  bash tests/test_extended_params_sim2.sh [--quick] [--full]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/UMI-nea/UMI-nea"
BIG_SIM2="/home/jdeng/claude/benchmark/sim2.input"

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

# sweep_exec <base_label> <ref_sorted> <input> <algo_args...>
sweep_exec() {
    local base_label="$1" ref="$2" input="$3"; shift 3
    local algo_args=("$@")

    # Thread sweep (pool=50000, prod=2000)
    for t in "${THREAD_LIST[@]}"; do
        local r="$TMPDIR/${base_label}_t${t}.sorted"
        run_one "${base_label}_t${t}" "${algo_args[@]}" -t "$t" -p 50000 -d 2000 -i "$input" > "$r"
        compare "t=$t  p=50000  d=2000  [$base_label]" "$r" "$ref"
    done

    if $FULL; then
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

        # Combined
        for combo in "${COMBO_LIST[@]}"; do
            read -r t p d <<< "$combo"
            local r="$TMPDIR/${base_label}_t${t}_p${p}_d${d}.sorted"
            run_one "${base_label}_t${t}_p${p}_d${d}" "${algo_args[@]}" -t "$t" -p "$p" -d "$d" -i "$input" > "$r"
            compare "t=$t  p=$p  d=$d  [$base_label]" "$r" "$ref"
        done
    fi
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
    local t0=$SECONDS
    run_one "${base_label}_ref" "${algo_args[@]}" -t 8 -p 50000 -d 2000 -i "$input" > "$ref"
    echo "  Reference built in $(( SECONDS - t0 ))s: $(wc -l < "$ref") output lines."
    sweep_exec "$base_label" "$ref" "$input" "${algo_args[@]}"
}

# ── setup ─────────────────────────────────────────────────────────────────────

[[ -x "$BINARY" ]]   || die "binary not found: $BINARY"
[[ -f "$BIG_SIM2" ]] || die "input not found: $BIG_SIM2"

QUICK=false; FULL=false
for arg in "${@}"; do
    [[ "$arg" == "--quick" ]] && QUICK=true
    [[ "$arg" == "--full"  ]] && FULL=true
done

if $QUICK; then
    THREAD_LIST=(2 16)
    POOL_LIST=(50000 100000)
    PROD_LIST=(2000 5000)
    COMBO_LIST=("16 100000 5000" "2 200000 50000")
else
    THREAD_LIST=(2 8 16 32)
    POOL_LIST=(50000 100000 200000)
    PROD_LIST=(2000 5000 10000)
    COMBO_LIST=(
        "4  50000  5000"
        "16 200000 5000"
        "32 50000  10000"
        "2  200000 50000"
    )
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== UMI-nea extended param test — benchmark/sim2.input  [$(date)] ==="
echo "    Dataset: $BIG_SIM2 ($(wc -l < "$BIG_SIM2") UMIs, 24-mer)"
$QUICK && echo "(quick mode)"
$FULL  && echo "(full mode: pool/prod sweeps enabled)"
echo "    Thread list: ${THREAD_LIST[*]}"

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 1: UMI length sweep  (-l 8/12/18/24/30)"
echo "Tests trimming path; input UMIs are 24-mers."
echo "════════════════════════════════════════════════════════════"
for l in 8 12 18 24 30; do
    run_config "big_sim2_l${l}" "$BIG_SIM2" -l "$l" -e 0.01 -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 2: Min-founder threshold  (-f 2/3/5)"
echo "════════════════════════════════════════════════════════════"
for f in 2 3 5; do
    run_config "big_sim2_f${f}" "$BIG_SIM2" -l 24 -e 0.01 -f "$f"
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 3: Explicit distance  (-m 1/2/3)"
echo "════════════════════════════════════════════════════════════"
for m in 1 2 3; do
    run_config "big_sim2_m${m}" "$BIG_SIM2" -l 24 -m "$m" -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 4: Greedy mode  (-g)"
echo "════════════════════════════════════════════════════════════"
run_config "big_sim2_greedy" "$BIG_SIM2" -l 24 -e 0.01 -f 1 -g

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════"
echo "Section 5: Error rate sweep  (-e 0.001/0.005/0.01/0.02/0.05)"
echo "════════════════════════════════════════════════════════════"
for e in 0.001 0.005 0.01 0.02 0.05; do
    tag="e${e//./p}"
    run_config "big_sim2_${tag}" "$BIG_SIM2" -l 24 -e "$e" -f 1
done

# ══════════════════════════════════════════════════════════════════════════════
TOTAL=$(( PASS+FAIL ))
echo ""
echo "=== Results: $PASS/$TOTAL passed  [$(date)] ==="
if [[ $FAIL -gt 0 ]]; then
    echo "FAILED configurations:"
    for e in "${ERRORS[@]}"; do echo "  - $e"; done
    exit 1
fi
echo "All reproducibility checks PASSED."
