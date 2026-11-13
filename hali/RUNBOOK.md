# MEXGB on Hali — runbook

Follow the steps in order. Each says what to run, what you should see, and what
to do if it goes wrong. The structure is the same as the mmCTBN campaign: one
config file, one `submit.sh`, everything inside an Apptainer image, an array of
small idempotent tasks, and a combine stage.

**What runs where.** The analysis code is the same as the local version
(`analysis/`). `analysis/R/registry.R` defines every study (simulation scenario or
real data set + outcome), its seeds and its methods, and both the local scripts
and the Hali tasks call it. So a split computed on Hali is identical to the same
split computed on your PC: this was checked on the smoke design, where all
studies matched exactly (every metric and every prediction).

**What a task is.** One task = one (study, repetition, fold), running every method
of that study on one core. Studies that are expensive per split (Sprockett: MERF and
REEMforest on p = 1879) are split further into one task per method. The full paper
design is 36 studies × 100 splits ≈ **5,200 tasks** (~430 core-hours, mostly
Sprockett). At 200 tasks running at once, that is a few hours of wall-clock time
instead of days on a PC.

Do the **smoke test (Step 6)** once, the **rehearsal (Steps 7–10)** next, and then
the **full run (Step 11)**.

---

## Step 1 — Put the code on Hali

```bash
cd ~
unzip MEXGB_paper1_code.zip
cd MEXGB_paper1
ls
# expect: MEXGB/  MEXGB_0.5.0.tar.gz  analysis/  hali/  README.md
```

`HALI_ROOT` is this directory (the one containing `MEXGB/`, `analysis/` and `hali/`).
`submit.sh` works it out from its own location; you never set it by hand.

---

## Step 2 — Load Apptainer

Hali's R modules do not work (their binaries are linked against libraries that are not
installed), so everything runs inside a container with its own R stack.

```bash
module purge
module load apptainer/1.4.0-rc.1
apptainer --version
```

---

## Step 3 — Build the image (once)

```bash
cd ~/MEXGB_paper1/hali/container
./build_image.sh            # uses ../../MEXGB_0.5.0.tar.gz
```

Run it on the **login node** (compute nodes usually have no network). It downloads
`rocker/r-ver:4.4.3`, installs MEXGB, every competitor package (xgboost, MEGB, gbm,
randomForest, LongituRF, gpboost, mboost, glmmLasso, lme4, nlme, survival,
randomForestSRC, glmnet, coxme) and the Bioconductor data packages for the
multi-study survival sets (breastCancerVDX/TRANSBIG/MAINZ/UNT, curatedOvarianData).
gpboost and xgboost compile large C++ sources: allow 30–60 minutes and ≥ 8 GB RAM.

On success:

```
BUILT: /.../hali/container/mexgb.sif
R version 4.4.3 (2025-02-28)
MEXGB 0.5.0 | MEGB 0.4 (crossfit: yes )
```

The build stops with `FAILED TO INSTALL: ...` if any package is missing. That is
deliberate: a missing package is cheaper to find now than at task 3,000.

**If it says "fakeroot build failed; Route B"**: fakeroot is disabled on many sites.
The script falls back to pulling the base image and installing the packages into a
host library (`$R_LIBS_USER`), which needs no privilege. Keep `R_LIBS_USER` as set in
`config.sh`. If the Bioconductor data packages fail in Route B, only `BreastMulti` and
`OvarianMulti` are skipped; everything else runs.

**Disk:** if `$HOME` is quota'd, point the build cache at scratch first:

```bash
export APPTAINER_CACHEDIR=/gpfs/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/gpfs/scratch/$USER/apptainer_tmp
```

---

## Step 4 — Edit `hali/config.sh`

This is the only file you edit. Set:

```bash
export SLURM_ACCOUNT="your-account"   # "" if Hali does not use accounts
export SLURM_PARTITION="compute"      # the real partition name
export SLURM_QOS=""                   # "" if unused
export MAIL_USER="you@uea.ac.uk"      # "" for no mail
export ARRAY_CONCURRENCY=200          # max tasks running at once (your fair share)
```

Everything about the *analysis* (methods, hyper-parameters, scenarios, target AUC,
survival design) stays in `analysis/config.R`, exactly as for the local run.
`config.sh` only chooses the run size, the output location and the resources.

---

## Step 5 — Preflight

```bash
cd ~/MEXGB_paper1/hali
source env/smoke.sh
./submit.sh plan
```

Runs in the foreground in seconds and submits nothing. In the `ENVIRONMENT` block:

- `execution environment` must say **Apptainer container**. If it says *native*, the
  image was not found or `USE_CONTAINER` is off.
- every required package must say `OK`. The Bioconductor sets may say *missing*:
  then only `BreastMulti` / `OvarianMulti` are skipped.
- `SprockettData.csv` and `readmission.rda` must say `OK` (they ship in `analysis/data/`).

---

## Step 6 — Smoke test (minutes)

```bash
source env/smoke.sh
./submit.sh all
```

`all` runs `plan`, then submits `prep` → the task array (starts when prep succeeds) →
`combine` (starts when every array task has finished), and prints the job ids.
The smoke design is tiny (2 folds, p ≤ 30, 30 trees; 46 tasks), so it finishes in a few
minutes once scheduled. Check:

```bash
./submit.sh status                         # done 46, failed 0, outstanding 0
ls ../output_smoke/results/tables/         # *.tex, *.csv
ls ../output_smoke/results/figures/        # *.pdf
```

It checks the image, the bindings, SLURM and the whole pipeline before you spend
real allocation. *If tasks fail*, read one task log in `output_smoke/logs/tasks/`
(each task writes its own) and the SLURM output in `hali/logs/slurm/`.

---

## Step 7 — Rehearsal: select it

```bash
source env/rehearsal.sh
```

It echoes:

```
HALI_ROOT = /gpfs/home/YOU/MEXGB_paper1
OUT_DIR   = /gpfs/home/YOU/MEXGB_paper1/output_rehearsal
```

The rehearsal runs **every** study with 1 repetition × 10 folds (10 splits per study,
~520 tasks) to **measure** the time and memory of every task. **Re-run this `source`
in every new shell**: otherwise the next stage writes to the wrong `OUT_DIR`.

---

## Step 8 — Rehearsal: plan and submit

```bash
./submit.sh plan     # read TASK ALGEBRA and RESOURCE PROJECTION
./submit.sh all      # prep -> array -> combine, chained
```

`plan` prints, per block, the number of studies and tasks, the estimated core-hours
and the longest task. These are estimates from single-core timings; the rehearsal
replaces them with measurements.

You can also run the stages one at a time: `./submit.sh prep`, then (after it
finishes) `./submit.sh fit`, then `./submit.sh combine`.

---

## Step 9 — Monitor and resubmit

```bash
./submit.sh status
```

```
tasks         : 520
  done        : 498
  failed      : 0
  outstanding : 22
```

One task's detail: `tail -f $OUT_DIR/logs/tasks/task_00042_*.log`.

Tasks are **idempotent**: `./submit.sh fit` submits only the tasks without a `.done`
marker, so a cut-short array (walltime, memory, node failure) is simply resubmitted.

```bash
./submit.sh fit                      # everything outstanding
./submit.sh fit --only 88,142,377    # specific task ids
FORCE=1 ./submit.sh fit --only 88    # recompute a completed task
./submit.sh local 88                 # run task 88 in the foreground (debugging)
```

If tasks fail on **memory**, raise `FIT_MEM`; on **walltime**, raise `FIT_TIME`. Then
`./submit.sh fit` and `./submit.sh combine` again.

A failure of one *method* inside a task (e.g. GLMM-PQL not converging on a separable
split) is **not** a task failure: it is recorded in the results (column `error`) and
counted in the tables, exactly as in the local run.

---

## Step 10 — Combine, and read the measurements

`all` already submitted combine. To rerun it after resubmissions:

```bash
./submit.sh combine
```

It assembles the task files into `results/raw/<study>.rds` (the same files a local run
writes), runs the MEXGB convergence study, then `analysis/04_tables_figures.R`.
Everything lands in `$OUT_DIR/results/`:

| file | contents |
|---|---|
| `results_summary.txt` | every table in plain text: **send this back** |
| `tables/*.tex`, `*.csv` | the paper's tables (same names and labels as before) |
| `figures/*.pdf` | the paper's figures |
| `tables/tab_task_timing.csv` | measured seconds and peak RSS **per task** |
| `tables/tab_task_timing_by_study.csv` | the same, per study |
| `tables/tab_task_status.csv` | done / failed per task |

If combine reports missing tasks, it prints a ready-made `./submit.sh fit --only ...`.

**Record the measured maxima.** In `tab_task_timing_by_study.csv`, take the largest
`max_sec` and `max_rss_gb`. The per-task cost does **not** grow in the full run (the
splits are the same size); only the number of tasks grows (×10). So:

```bash
export FIT_TIME=...   # measured max seconds ×2, as hh:mm:ss
export FIT_MEM=...    # measured max RSS ×1.5, rounded up (e.g. 6G)
```

Every env file sets K, the repetitions, `FIT_MEM` and `FIT_TIME` itself, so nothing
leaks from one run to the next: export your own values **after** `source env/...`.

---

## Step 11 — The full run

```bash
cd ~/MEXGB_paper1/hali
source env/full.sh           # 10 repetitions x 10 folds = 100 splits per study
export FIT_TIME=04:00:00     # from Step 10 -- AFTER sourcing (the env file sets its own)
export FIT_MEM=8G            # from Step 10
./submit.sh plan
./submit.sh all
./submit.sh status           # repeat until outstanding = 0
./submit.sh combine          # if anything was resubmitted after the chained combine
```

The full run writes to `output_full/` and does not touch the rehearsal's output.
Send back `output_full/results/results_summary.txt`, `tables/` and `figures/`.

### Running only part of the campaign

`BLOCKS` selects blocks, e.g. only the survival work:

```bash
export BLOCKS=ssim,sreal
source env/full.sh && ./submit.sh plan && ./submit.sh all
```

Blocks: `csim` continuous simulation · `bsim` binary simulation · `ssim` survival
simulation · `real` BodyWeight/Soybean/Sprockett (continuous + binary) · `sreal`
Retinopathy/Readmission/BreastMulti/OvarianMulti (survival).
Run different blocks under different `RUN_LABEL`s (and so different `OUT_DIR`s),
or re-run `plan` + `combine` in the same `OUT_DIR` with the union of blocks.

---

## Quick reference

| Command | Does |
|---|---|
| `./submit.sh plan` | preflight, task table, resource projection. Submits nothing |
| `./submit.sh prep` | calibrations + real data, computed once |
| `./submit.sh fit` | submit every task without a `.done` marker |
| `./submit.sh fit --only 1,2,3` | submit specific task ids |
| `./submit.sh combine` | assemble → convergence → tables → figures |
| `./submit.sh all` | plan, then prep → fit → combine chained by SLURM dependencies |
| `./submit.sh status` | done / failed / outstanding |
| `./submit.sh local 88` | run task 88 in the foreground |

| Variable (`config.sh` / env files) | Meaning |
|---|---|
| `RUN_LABEL`, `OUT_DIR` | name of the run and its output directory |
| `MEXGB_N_REP`, `MEXGB_N_REP_REAL`, `MEXGB_K` | repetitions (simulation / real data) and folds |
| `BLOCKS` | which blocks to run |
| `SPLIT_THRESHOLD` | studies costing more than this (s per split) get one task per method |
| `FIT_MEM`, `FIT_TIME` | per-task resources (one core each) |
| `ARRAY_CONCURRENCY`, `ARRAY_CHUNK` | max running tasks; tasks per submitted array |
| `FORCE=1`, `STRICT_TASK=1` | recompute completed tasks; make failed tasks exit non-zero |
| `USE_CONTAINER`, `CONTAINER_IMAGE`, `CONTAINER_BINDS` | container use, image path, extra binds |
| `RUN_CONVERGENCE=0` | skip the convergence study in combine |

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `no task table` | run `./submit.sh plan` first |
| `task table and registry disagree` | `analysis/config.R` or `BLOCKS` changed after `plan`: run `plan` again, then `prep` |
| Stage writes to the wrong directory | new shell without `source env/<run>.sh` |
| `container image not found` | build it (Step 3) |
| `libopenblas.so.0: cannot open shared object file` | a stage ran outside the container: `USE_CONTAINER=1`, and `plan` must say *Apptainer container* |
| Job cannot see files on /gpfs | add the filesystem root to `CONTAINER_BINDS` in `config.sh` |
| `apptainer: command not found` in a job | `HALI_MODULES` must list the apptainer module |
| Many tasks fail quickly | read one task log in `$OUT_DIR/logs/tasks/`; usually memory (`FIT_MEM`) |
| Tasks fail with an empty task log | killed by SLURM: see `hali/logs/slurm/` for the OOM or timeout notice |
| `BreastMulti`/`OvarianMulti` tasks fail with *install Bioconductor package* | the image lacks the Bioconductor data packages: rebuild (Step 3), or run without them |
| sbatch rejects the array size | lower `ARRAY_CHUNK` below the site's `MaxArraySize` |
