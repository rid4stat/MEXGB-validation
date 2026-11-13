# =============================================================================
# hali/run_r.sh -- single definition of "how do I run R here"
# =============================================================================
# Sourced by submit.sh and by every SLURM job. run_r() executes an Rscript
# inside the Apptainer image (USE_CONTAINER=1, required on Hali) or natively.
# One definition, so the login node, the array tasks and the debug path cannot
# disagree about which R they use.
# =============================================================================

run_r() {
  if [[ "${USE_CONTAINER:-1}" != "1" ]]; then
    Rscript "$@"
    return
  fi

  if [[ ! -f "${CONTAINER_IMAGE}" ]]; then
    echo "ERROR: container image not found: ${CONTAINER_IMAGE}" >&2
    echo "Build it first:  cd hali/container && ./build_image.sh" >&2
    return 1
  fi

  local binds=() b
  for b in "${HALI_ROOT}" "${OUT_DIR}" ${CONTAINER_BINDS:-}; do
    [[ -n "$b" && -e "$b" ]] && binds+=(--bind "$b")
  done

  # Host-library mode (fakeroot unavailable at build time): packages live outside
  # the image, so bind that directory and put it on the library path.
  local kind_file="${HALI_ROOT}/hali/container/image_kind"
  local rlibs_env=()
  if [[ -f "$kind_file" && "$(cat "$kind_file")" == "hostlib" ]]; then
    [[ -d "${R_LIBS_USER}" ]] && binds+=(--bind "${R_LIBS_USER}")
    rlibs_env+=(--env "R_LIBS_USER=${R_LIBS_USER}")
  fi

  # R needs a writable temporary directory, otherwise it dies at start-up with
  # "Fatal error: cannot create 'R_TempDir'" (a compute node whose /tmp is full or not
  # writable). --cleanenv drops the host TMPDIR, so choose one explicitly: the scheduler's
  # TMPDIR, then /tmp, then the run's output directory on GPFS -- the first one in which a
  # file can actually be written -- and bind it into the container.
  local rtmp="" c
  for c in "${TMPDIR:-}" /tmp "${OUT_DIR}/tmp"; do
    [[ -n "$c" ]] || continue
    local d="$c/mexgb_${SLURM_JOB_ID:-login}_${SLURM_ARRAY_TASK_ID:-0}_$$"
    if mkdir -p "$d" 2>/dev/null && ( : > "$d/.w" ) 2>/dev/null; then rm -f "$d/.w"; rtmp="$d"; break; fi
    rmdir "$d" 2>/dev/null || true
  done
  if [[ -z "$rtmp" ]]; then
    echo "ERROR: no writable temporary directory (tried TMPDIR, /tmp, ${OUT_DIR}/tmp) on $(hostname)" >&2
    return 1
  fi
  binds+=(--bind "$rtmp")

  local rc=0
  apptainer exec --cleanenv \
    "${binds[@]}" \
    --env TMPDIR="$rtmp" \
    --env HALI_ROOT="${HALI_ROOT}" \
    --env OUT_DIR="${OUT_DIR}" \
    --env RUN_LABEL="${RUN_LABEL}" \
    --env MEXGB_OUT="${MEXGB_OUT}" \
    --env MEXGB_N_REP="${MEXGB_N_REP:-}" \
    --env MEXGB_N_REP_REAL="${MEXGB_N_REP_REAL:-}" \
    --env MEXGB_K="${MEXGB_K:-}" \
    --env MEXGB_SMOKE="${MEXGB_SMOKE:-}" \
    --env MEXGB_N_CORES="${MEXGB_N_CORES:-1}" \
    --env BLOCKS="${BLOCKS:-}" \
    --env SPLIT_THRESHOLD="${SPLIT_THRESHOLD:-1800}" \
    --env RUN_CONVERGENCE="${RUN_CONVERGENCE:-1}" \
    --env ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-}" \
    --env FIT_TIME="${FIT_TIME:-}" \
    --env FORCE="${FORCE:-}" \
    --env STRICT_TASK="${STRICT_TASK:-}" \
    --env SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-}" \
    --env SLURM_JOB_ID="${SLURM_JOB_ID:-}" \
    --env SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-}" \
    --env OMP_NUM_THREADS=1 \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    "${rlibs_env[@]}" \
    "${CONTAINER_IMAGE}" Rscript "$@" || rc=$?
  rm -rf "$rtmp" 2>/dev/null || true
  return $rc
}
