#!/bin/bash
#SBATCH --job-name=mexgb_fit
#SBATCH --output=logs/slurm/%x_%A_%a.out
#SBATCH --error=logs/slurm/%x_%A_%a.err
source "${HALI_ROOT}/hali/slurm/_common.sh"
cd "${HALI_ROOT}/hali"
# The array index points at a line of the task list written by submit.sh
# (index + offset, so indices stay below SLURM's MaxArraySize).
LINE=$(( SLURM_ARRAY_TASK_ID + ${TASK_OFFSET:-0} ))
TASK_ID=$(sed -n "${LINE}p" "${TASK_LIST:?TASK_LIST not set}")
[[ -n "$TASK_ID" ]] || { echo "no task on line $LINE of $TASK_LIST"; exit 0; }
run_r R/02_fit_task.R "${TASK_ID}"
