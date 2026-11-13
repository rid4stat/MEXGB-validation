# Sourced by every SLURM script. Loads modules, sets threading, echoes context.
set -euo pipefail
source "${HALI_ROOT}/hali/config.sh"
source "${HALI_ROOT}/hali/run_r.sh"
hali_load_modules
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "=== $(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"
echo "=== run=${RUN_LABEL} out=${OUT_DIR}"
