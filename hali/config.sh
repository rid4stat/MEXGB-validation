#!/usr/bin/env bash
# =============================================================================
# hali/config.sh -- central configuration for the MEXGB campaign on Hali
# =============================================================================
# Sourced by submit.sh and by every SLURM job. Edit here, nowhere else.
# The analysis itself (methods, hyper-parameters, scenarios) is configured in
# analysis/config.R; the variables below only choose the run SIZE, the output
# location and the cluster resources, and reach R through MEXGB_* environment
# overrides read at the end of analysis/config.R -- so the shell and R views
# cannot disagree.
# =============================================================================

# ---- run identity -----------------------------------------------------------
# REHEARSAL : 1 repetition  x 10 folds =  10 splits per study (measures cost)
# FULL      : 10 repetitions x 10 folds = 100 splits per study (the paper)
export RUN_LABEL="${RUN_LABEL:-rehearsal}"
export MEXGB_N_REP="${MEXGB_N_REP:-1}"          # simulation repetitions
export MEXGB_N_REP_REAL="${MEXGB_N_REP_REAL:-1}" # real-data repetitions
export MEXGB_K="${MEXGB_K:-10}"                  # folds
export MEXGB_SMOKE="${MEXGB_SMOKE:-}"            # 1 = tiny smoke design (minutes)

# Which blocks to run (comma-separated): csim, bsim, ssim, real, sreal
#   csim = continuous simulation      bsim = binary simulation
#   ssim = survival simulation        real = BodyWeight/Soybean/Sprockett (cont + binary)
#   sreal = Retinopathy/Readmission/BreastMulti/OvarianMulti (survival)
export BLOCKS="${BLOCKS:-csim,bsim,ssim,real,sreal}"

# A study whose estimated cost per split exceeds this many seconds is split into
# one task per METHOD (e.g. Sprockett, where MERF/REEMforest on p = 1879 dominate);
# otherwise one task runs every method of one (study, repetition, fold).
export SPLIT_THRESHOLD="${SPLIT_THRESHOLD:-1800}"

# ---- paths ------------------------------------------------------------------
# HALI_ROOT is the directory containing MEXGB/, analysis/ and hali/.
export HALI_ROOT="${HALI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export OUT_DIR="${OUT_DIR:-$HALI_ROOT/output_${RUN_LABEL}}"
export MEXGB_OUT="$OUT_DIR/results"              # raw/, tables/, figures/ (as locally)
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/hali/4.4.3}"

# ---- SLURM ------------------------------------------------------------------
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"
export SLURM_QOS="${SLURM_QOS:-}"
export MAIL_USER="${MAIL_USER:-}"

# prep: calibration of the binary/survival simulations (seconds) + loading and
# saving the real data sets once (Sprockett CSV, Bioconductor expression sets).
export PREP_MEM="${PREP_MEM:-16G}"
export PREP_CPUS="${PREP_CPUS:-1}"
export PREP_TIME="${PREP_TIME:-02:00:00}"

# Per array task: ONE core (every method is single-threaded, and one core keeps
# the timing table honest). Most tasks need < 2 GB; the high-dimensional real
# data (Sprockett, multi-study expression) need more.
export FIT_MEM="${FIT_MEM:-8G}"
export FIT_CPUS="${FIT_CPUS:-1}"
export FIT_TIME="${FIT_TIME:-12:00:00}"

# combine: assembles the per-task files, re-runs the small MEXGB convergence
# study (05_convergence.R) and writes every table and figure.
export COMBINE_MEM="${COMBINE_MEM:-32G}"
export COMBINE_CPUS="${COMBINE_CPUS:-4}"
export COMBINE_TIME="${COMBINE_TIME:-06:00:00}"
export RUN_CONVERGENCE="${RUN_CONVERGENCE:-1}"

# SLURM caps both the array size (MaxArraySize, often 1001) and concurrency.
# ARRAY_CHUNK splits the task table into contiguous submitted ranges;
# ARRAY_CONCURRENCY is the %N cap on simultaneously running tasks.
export ARRAY_CHUNK="${ARRAY_CHUNK:-1000}"
export ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-200}"

# ---- container --------------------------------------------------------------
# Hali's R modules are unusable (linked against libraries that are not installed),
# so -- as for the mmCTBN campaign -- everything runs inside an Apptainer image
# that carries its own R stack, MEXGB and every competitor package.
export USE_CONTAINER="${USE_CONTAINER:-1}"
export CONTAINER_IMAGE="${CONTAINER_IMAGE:-$HALI_ROOT/hali/container/mexgb.sif}"
# Apptainer maps $HOME and the CWD but NOT /gpfs.
export CONTAINER_BINDS="${CONTAINER_BINDS:-/gpfs}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRATCH:-/tmp}/${USER:-$(id -un)}/apptainer_cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${SCRATCH:-/tmp}/${USER:-$(id -un)}/apptainer_tmp}"

# ---- module environment -----------------------------------------------------
export HALI_MODULES="${HALI_MODULES:-apptainer/1.4.0-rc.1}"

hali_load_modules() {
  if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null || true
    for m in $HALI_MODULES; do module load "$m" || echo "WARN: module $m failed"; done
  fi
}
