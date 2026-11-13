#!/usr/bin/env bash
# =============================================================================
# hali/submit.sh -- single entry point for the MEXGB campaign on Hali
# =============================================================================
#
#   ./submit.sh plan               preflight + task table + resource projection (foreground)
#   ./submit.sh prep               calibrations + real data, computed once (job)
#   ./submit.sh fit                submit every task without a .done marker
#   ./submit.sh fit --only 12,88   resubmit specific task ids
#   ./submit.sh combine            assemble -> convergence study -> tables -> figures (job)
#   ./submit.sh all                plan, then prep -> fit -> combine chained with
#                                  SLURM dependencies (one command for the whole campaign)
#   ./submit.sh status             done / failed / outstanding
#   ./submit.sh local <id>         run one task in the foreground (debugging)
#
# Without SLURM (no sbatch on PATH) every stage runs in the foreground, in order.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export HALI_ROOT="${HALI_ROOT:-$(cd "$HERE/.." && pwd)}"
source "$HERE/config.sh"
source "$HERE/run_r.sh"

mkdir -p "$OUT_DIR"/{prep,fits,logs/tasks} "$MEXGB_OUT"/{raw,tables,figures}
mkdir -p "$HERE/logs/slurm"
cd "$HERE"

have_slurm() { command -v sbatch >/dev/null 2>&1; }

sb_common() {
  local args=()
  [[ -n "${SLURM_ACCOUNT}"   ]] && args+=(--account="${SLURM_ACCOUNT}")
  [[ -n "${SLURM_PARTITION}" ]] && args+=(--partition="${SLURM_PARTITION}")
  [[ -n "${SLURM_QOS}"       ]] && args+=(--qos="${SLURM_QOS}")
  [[ -n "${MAIL_USER}"       ]] && args+=(--mail-user="${MAIL_USER}" --mail-type=FAIL,END)
  printf '%s\n' "${args[@]+"${args[@]}"}"
}

submit() {
  local script="$1"; shift
  mapfile -t common < <(sb_common)
  sbatch --parsable "${common[@]+"${common[@]}"}" "$@" \
    --export=ALL,HALI_ROOT="$HALI_ROOT" "$script"
}

n_tasks() {
  local f="$OUT_DIR/prep/tab_plan_tasks.csv"
  [[ -f "$f" ]] || { echo 0; return; }
  echo $(( $(wc -l < "$f") - 1 ))
}

stage_plan() { echo ">>> PLAN"; run_r R/01_plan.R; }

stage_prep() {
  echo ">>> PREP" >&2
  if have_slurm; then
    submit slurm/prep.slurm --mem="$PREP_MEM" --cpus-per-task="$PREP_CPUS" \
      --time="$PREP_TIME" --chdir="$HERE" "$@"
  else
    run_r R/00_prep.R
  fi
}

# Submit the task ids listed in a file as array(s). Array indices are always
# 1..chunk (below SLURM's MaxArraySize); the job maps index + offset to a line of
# the list. Prints the submitted job ids on stdout.
submit_list() {
  local list="$1"; shift
  local n; n=$(wc -l < "$list")
  local res=(--mem="$FIT_MEM" --cpus-per-task="$FIT_CPUS" --time="$FIT_TIME" --chdir="$HERE")
  local from=1
  while [[ $from -le $n ]]; do
    local to=$(( from + ARRAY_CHUNK - 1 )); [[ $to -gt $n ]] && to=$n
    local id
    # Passed through the environment (--export=ALL in submit()). A second --export
    # option would REPLACE the first one -- sbatch keeps only the last.
    export TASK_LIST="$list" TASK_OFFSET=$(( from - 1 ))
    id=$(submit slurm/fit_array.slurm --array="1-$(( to - from + 1 ))%${ARRAY_CONCURRENCY}" \
          "${res[@]}" "$@")
    echo "    list lines ${from}-${to} -> job $id" >&2
    echo "$id"
    from=$(( to + 1 ))
  done
}

# fit: submit every task that has no .done marker (all of them on the first call).
# --only 3,17,88 submits exactly those task ids. Remaining args go to sbatch
# (e.g. --dependency=afterok:<prep job>).
stage_fit() {
  local only="" dep=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --only) only="$2"; shift 2 ;;
      *) dep+=("$1"); shift ;;
    esac
  done
  local N; N=$(n_tasks)
  if [[ "$N" -eq 0 ]]; then
    echo "!! no task table; run './submit.sh plan' (and prep) first." >&2
    exit 1
  fi
  local list="$OUT_DIR/prep/submit_$(date +%Y%m%d_%H%M%S)_$$.txt"
  if [[ -n "$only" ]]; then
    tr ',' '\n' <<< "$only" | sed '/^$/d' > "$list"
  elif find "$OUT_DIR/fits" -name '*.done' -print -quit 2>/dev/null | grep -q .; then
    run_r R/04_outstanding.R "$list" > /dev/null      # resubmission: outstanding tasks only
  else
    seq 1 "$N" > "$list"                              # first submission: everything
  fi
  local n; n=$(wc -l < "$list")
  echo ">>> FIT ($n of $N tasks)" >&2
  [[ "$n" -eq 0 ]] && { echo "    nothing outstanding" >&2; return; }

  if ! have_slurm; then
    local i
    while read -r i; do echo "--- task $i" >&2; run_r R/02_fit_task.R "$i"; done < "$list"
    return
  fi
  submit_list "$list" "${dep[@]+"${dep[@]}"}"
}

stage_combine() {
  echo ">>> COMBINE" >&2
  if have_slurm; then
    submit slurm/combine.slurm --mem="$COMBINE_MEM" --cpus-per-task="$COMBINE_CPUS" \
      --time="$COMBINE_TIME" --chdir="$HERE" "$@"
  else
    run_r R/03_combine.R
  fi
}

stage_status() {
  local slurm_dir="$HERE/logs/slurm"
  local N; N=$(n_tasks)
  local done_n fail_n
  done_n=$(find "$OUT_DIR/fits" -name '*.done' 2>/dev/null | wc -l)
  fail_n=$(find "$OUT_DIR/fits" -name '*.fail' 2>/dev/null | wc -l)
  echo "run           : $RUN_LABEL   (N_REP=$MEXGB_N_REP, N_REP_REAL=$MEXGB_N_REP_REAL, K=$MEXGB_K)"
  echo "out           : $OUT_DIR"
  echo "prep          : $(ls "$OUT_DIR"/prep/registry.rds 2>/dev/null || echo MISSING)"
  echo "tasks         : $N"
  echo "  done        : $done_n"
  echo "  failed      : $fail_n"
  echo "  outstanding : $(( N - done_n - fail_n ))"
  if [[ "$fail_n" -gt 0 ]]; then
    echo "failed tasks (first 30):"
    find "$OUT_DIR/fits" -name '*.fail' -printf '  %P\n' | head -30
  fi
  have_slurm && squeue -u "${USER:-$(id -un)}" -o '%.12i %.18j %.9T %.10M %R' 2>/dev/null | head -20 || true
  # Nothing done, nothing failed. Either the fit array has not started yet (normal while
  # prep is running), or its tasks died before R started -- only then show an error log,
  # and only one written by THIS run (newer than prep/registry.rds). SLURM's own
  # energy-accounting messages (slurm_get_node_energy / _get_joules_task) are harmless
  # node-plugin noise and are filtered out.
  local me="${USER:-$(id -un)}" noise='slurm_get_node_energy|_get_joules_task|energy'
  # (1) prep finished without writing prep/registry.rds -> it FAILED; show its log
  if [[ ! -f "$OUT_DIR/prep/registry.rds" ]] && have_slurm && \
     ! squeue -u "$me" -h -n mexgb_prep -o '%T' 2>/dev/null | grep -q .; then
    local plog; plog=$(ls -t "$slurm_dir"/mexgb_prep_*.err 2>/dev/null | head -1 || true)
    echo "!! prep is not queued or running and prep/registry.rds is MISSING: the prep job failed."
    if [[ -n "$plog" ]]; then
      local pid; pid=$(basename "$plog" .err); pid=${pid##*_}
      echo "   prep job $pid -- last lines of $plog:"
      grep -Ev "$noise" "$plog" | tail -15 | sed 's/^/  | /'
      [[ -f "${plog%.err}.out" ]] && { echo "   ... and of ${plog%.err}.out:"; tail -8 "${plog%.err}.out" | sed 's/^/  | /'; }
      command -v sacct >/dev/null 2>&1 && sacct -j "$pid" -X --format=JobID,State,ExitCode,Elapsed,Timelimit,MaxRSS,ReqMem 2>/dev/null | sed 's/^/   /'
    fi
  fi
  # (2) jobs whose dependency failed never start
  if have_slurm && squeue -u "$me" -h -o '%R' 2>/dev/null | grep -q DependencyNeverSatisfied; then
    echo "!! jobs are waiting on a dependency that FAILED and will never start. Cancel them,"
    echo "   fix the cause shown above, then resubmit:"
    echo "     scancel -n mexgb_fit; scancel -n mexgb_combine; ./submit.sh all"
    return 0
  fi
  if [[ "$done_n" -eq 0 && "$fail_n" -eq 0 && "$N" -gt 0 ]]; then
    local fit_active=0
    if have_slurm && squeue -u "${USER:-$(id -un)}" -h -n mexgb_fit,mexgb_prep -o '%T' 2>/dev/null \
         | grep -q .; then fit_active=1; fi
    local reg="$OUT_DIR/prep/registry.rds" last=""
    if [[ -f "$reg" ]]; then
      last=$(find "$slurm_dir" -name 'mexgb_fit_*.err' -newer "$reg" -printf '%T@ %p\n' 2>/dev/null \
               | sort -rn | head -1 | cut -d' ' -f2- || true)
    fi
    if [[ -n "$last" ]] && grep -Evq "$noise" "$last" 2>/dev/null; then
      echo "no task has written a result yet; latest error log of this run ($last):"
      grep -Ev "$noise" "$last" | tail -5 | sed 's/^/  | /'
    elif [[ "$fit_active" -eq 1 ]]; then
      echo "no task has finished yet (prep / fit jobs queued or running -- normal)."
    fi
  fi
}

case "${1:-}" in
  plan)    stage_plan ;;
  prep)    stage_prep ;;
  fit)     shift; stage_fit "$@" ;;
  combine) stage_combine ;;
  status)  stage_status ;;
  local)   shift; run_r R/02_fit_task.R "${1:?usage: submit.sh local <task_id>}" ;;
  all)
    stage_plan
    if have_slurm; then
      p=$(stage_prep | tail -1)
      echo "    prep    -> job $p"
      mapfile -t fids < <(stage_fit --dependency=afterok:"$p")
      deps=$(IFS=:; echo "${fids[*]}")
      c=$(stage_combine --dependency=afterany:"$deps" | tail -1)
      echo "    combine -> job $c (runs when every array task has finished)"
      echo "Monitor with:  ./submit.sh status"
    else
      stage_prep; stage_fit; stage_combine
    fi
    ;;
  *) sed -n '2,18p' "$0"; exit 1 ;;
esac
