#!/bin/bash
#SBATCH --job-name=mexgb_combine
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --error=logs/slurm/%x_%j.err
source "${HALI_ROOT}/hali/slurm/_common.sh"
cd "${HALI_ROOT}/hali"
run_r R/03_combine.R
