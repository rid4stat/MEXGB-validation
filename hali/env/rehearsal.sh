# =============================================================================
# hali/env/rehearsal.sh -- 1 repetition x 10 folds per study (10 splits)
# =============================================================================
#   source env/rehearsal.sh
# Runs every study at one tenth of the paper's size to MEASURE the per-task time
# and memory (tables/tab_task_timing.csv) before the full run. Every variable is
# EXPORTED, because submit.sh passes the environment to sbatch with --export=ALL.
# =============================================================================
export HALI_ROOT="${HALI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export RUN_LABEL=rehearsal
export MEXGB_N_REP=1
export MEXGB_N_REP_REAL=1
export MEXGB_SMOKE=
export MEXGB_K=10                 # set explicitly: a value left over from env/smoke.sh must not leak in
export OUT_DIR="${HALI_ROOT}/output_rehearsal"
export FIT_MEM=16G
export FIT_TIME=24:00:00
echo "HALI_ROOT = $HALI_ROOT"
echo "OUT_DIR   = $OUT_DIR"
echo "K = $MEXGB_K folds | N_REP = $MEXGB_N_REP | FIT_MEM = $FIT_MEM | FIT_TIME = $FIT_TIME"
echo "(to change FIT_MEM / FIT_TIME, export them AFTER sourcing this file)"
