# =============================================================================
# hali/env/smoke.sh -- tiny design (2 folds, p <= 30, few trees): checks the
# whole campaign end to end in minutes. Use it once after building the image.
# =============================================================================
export HALI_ROOT="${HALI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export RUN_LABEL=smoke
export MEXGB_SMOKE=1
export MEXGB_N_REP=1
export MEXGB_N_REP_REAL=1
export MEXGB_K=2
export OUT_DIR="${HALI_ROOT}/output_smoke"
export FIT_MEM=4G
export FIT_TIME=01:00:00
echo "HALI_ROOT = $HALI_ROOT"
echo "OUT_DIR   = $OUT_DIR"
echo "K = $MEXGB_K folds | N_REP = $MEXGB_N_REP | FIT_MEM = $FIT_MEM | FIT_TIME = $FIT_TIME"
echo "(to change FIT_MEM / FIT_TIME, export them AFTER sourcing this file)"
