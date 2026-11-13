# =============================================================================
# hali/env/full.sh -- the paper: 10 repetitions x 10 folds = 100 splits per study
# =============================================================================
#   source env/full.sh
# Set FIT_MEM / FIT_TIME from the REHEARSAL's measured maxima
# (output_rehearsal/results/tables/tab_task_timing.csv): per-task cost does not
# grow with the number of repetitions -- only the NUMBER of tasks does (10x).
# =============================================================================
export HALI_ROOT="${HALI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export RUN_LABEL=full
export MEXGB_N_REP=10
export MEXGB_N_REP_REAL=10
export MEXGB_SMOKE=
export MEXGB_K=10                 # set explicitly: a value left over from env/smoke.sh must not leak in
export OUT_DIR="${HALI_ROOT}/output_full"
export FIT_MEM=8G
export FIT_TIME=12:00:00
echo "HALI_ROOT = $HALI_ROOT"
echo "OUT_DIR   = $OUT_DIR"
echo "K = $MEXGB_K folds | N_REP = $MEXGB_N_REP | FIT_MEM = $FIT_MEM | FIT_TIME = $FIT_TIME"
echo "(to change FIT_MEM / FIT_TIME, export them AFTER sourcing this file)"
