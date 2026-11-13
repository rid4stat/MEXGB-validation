# MEXGB — package and reproducible analysis (continuous, binary and survival outcomes)

```
MEXGB/                 R package source (xgboost back-end); version 1.0.0
MEGB/                  MEGB 0.3 source (same cross-fitted E-step); MEGB_0.3.tar.gz built package
diagnostics/           before/after benchmark
analysis/
  config.R             ALL settings (design, hyper-parameters, cores, smoke mode)
  00_install.R         installs every package used (CRAN, LongituRF from CRAN archive, Bioconductor, MEXGB)
  01_sim_continuous.R  continuous simulation: main study + sensitivity sub-study
  02_sim_binary.R      binary simulation (calibrated target AUC + imposed prevalence)
  03_real_data.R       BodyWeight, Soybean, Sprockett — continuous and binary
  06_sim_survival.R    survival simulation (frailty Cox model, calibrated target C-index)
  07_real_survival.R   Retinopathy, Readmission, BreastMulti, OvarianMulti — survival
  08_data_description.R descriptive table of the real data (tab_real_data.csv/.tex)
  04_tables_figures.R  LaTeX tables + PDF figures + results_summary.txt (all outcomes)
  05_convergence.R     MEXGB log-likelihood evolution (analogue of MEGB paper Fig. 1): one figure per
                       outcome, fig_mexgb_loglik_{continuous,binary,survival}.pdf
  run_all.R            runs 01 → 02 → 03 → 06 → 07 → 05 → 04
  R/methods.R          one wrapper per method (each calls its own package)
  R/cv.R               metrics (MSE, AUC, Brier, C-index, IBS) + repeated K-fold CV driver
  R/real_data.R        BodyWeight / Soybean / Sprockett loaders
  R/real_survival.R    survival data loaders
  R/registry.R         THE list of studies (scenarios / data sets), seeds and methods,
                       shared by the local scripts and the Hali campaign
hali/                  SLURM + Apptainer campaign for the Hali cluster (see hali/RUNBOOK.md)
  data/SprockettData.csv   your Sprockett export (used directly)
  data/readmission.rda     frailtypack::readmission (so frailtypack need not be installed)
```

## Running on Hali (recommended for the full design)

`hali/` runs the same analysis as a SLURM campaign inside an Apptainer image,
`hali/submit.sh plan | prep | fit | combine | all | status | local <id>`, env files for
`smoke`, `rehearsal` (10 splits per study) and `full` (100 splits per study), and a
step-by-step **`hali/RUNBOOK.md`**. One array task = one (study, repetition, fold);
the full design is ~5,200 one-core tasks (~430 core-hours), a few hours at 200
concurrent tasks. Every study, seed and setting comes from `analysis/R/registry.R`,
shared with the local scripts. On the smoke design the Hali pipeline (run without
SLURM) reproduced the local results exactly for all studies.

```bash
cd hali/container && ./build_image.sh          # once, on the login node
cd .. && source env/smoke.sh && ./submit.sh all  # minutes: end-to-end check
source env/rehearsal.sh && ./submit.sh all      # measures time/memory per task
source env/full.sh && ./submit.sh all           # the paper
```

## How to run locally

```r
setwd("analysis")
source("00_install.R")              # once
Sys.setenv(MEXGB_SMOKE = "1"); source("run_all.R")   # ~1 min end-to-end check (writes results_smoke/)
Sys.setenv(MEXGB_SMOKE = "");  source("run_all.R")   # full run (writes results/)
```

From a terminal: `Rscript run_all.R` (Linux/macOS smoke: `MEXGB_SMOKE=1 Rscript run_all.R`).
Every scenario is cached in `results/raw/*.rds`; an interrupted run resumes where it stopped.
Workers: `CFG$N_CORES` (default = cores − 1).

## Design
* n = 50 subjects × 20 occasions (v0.4; was 20 × 10), t rescaled to [0,1]; x ~ N_p(0, Σ_W), (Σ_W)_gg' = 0.6^|g−g'|;
  p ∈ {6, 170, 2000}, 6 active predictors; linear f = 2x1−2x2+1.5x3−1.5x4+x5−x6, nonlinear
  f = 2 sin(πx1x2)+4(x3−0.5)²+2x4+exp(x5)+x6x1; correlated random intercept/slope
  (τ0² = 0.5, τ1² = 3, ρ_Z = 0.6); AR(1) within-subject errors σ = 0.5 (lag in visit number).
  Binary: f centred/scaled to SD 1, intercept solved so mean P(y=1) = 0.5. → `MEXGB::simLong()`.
* Sensitivity (continuous, nonlinear, p = 170): ρ_W ∈ {0, 0.2, 0.6, 0.8} × σ ∈ {0.5, 1.5}.
* Evaluation: every data set, simulated or real, is randomly split K times into a
  training and a test set. Each test set holds 30% of the measurements of EVERY subject
  (round(0.3 × n_i), at least 1, at most n_i − 1; simulation: 6 of 20 for testing, 14 for
  training). 10 simulated data sets × K = 10 splits = 100 splits per scenario; real data 10 × 10.
  Test subjects are seen in training, so every mixed model predicts with f(x) + z'b̂_i.
  Settings: `CFG$test_frac = 0.3`, `CFG$sim_cv_unit`, `CFG$real_cv_unit`, `CFG$real_surv$cv_unit`.
* Real data, all with ≥ 10 measurements per subject or cluster: continuous BodyWeight,
  Soybean (log outcome), Sprockett; binary PBC (hepatomegaly, low-dim) and OvarianDebulking
  (high-dim); survival Insemination (dairy herds, low-dim), BreastMulti and OvarianMulti
  (high-dim).
* Metrics: MSE; AUC, Brier, classification TPR/FPR (0.5); variable-selection TPR/FPR of the
  top-6 ranked predictors (non-zero importance only); training time.
* Real data: BodyWeight, Soybean (`nlme`), Sprockett Tsimane cohort
  (`microbiomeDataSets::SprockettTHData()`: feces, ≥3 visits, CLR of counts + 0.5, lowest-variance
  decile removed, Shannon diversity response, visit index as time). Binary = median split.
  Real-data binary AUC/Brier are computed per split. Every test set now contains all subjects;
  with subject-level CV the predictions are pooled over the folds of a repetition instead.

## Methods and packages
| Method | Package / function | Notes |
|---|---|---|
| **MEXGB** | `MEXGB::MEXGB` (`xgboost::xgb.train`) | same EM/BLUP loop as MEGB; cross-fitted E-step (V = 5), LMM start |
| XGBoost | `xgboost::xgb.train`, identical booster settings, no random effects | continuous + binary (XGB-Cox for survival); the D → 0 limit of MEXGB |
| MEGB | `MEGB::MEGB` 0.3 (gbm) | continuous only; cross-fitted E-step |
| MERF, REEMforest | `LongituRF::MERF`, `LongituRF::REEMforest` (`sto = "none"`) | continuous only |
| RF | `randomForest` | |
| GBM | `gbm` | |
| GPBoost | `gpboost::gpb.train` + `GPModel` (random intercept + slope) | |
| mboost | `mboost::glmboost` (centred) | |
| glmmLasso | `glmmLasso::glmmLasso` (λ = cv.lambda, random intercept + slope) | p < N_train |
| LMM | `nlme::lme` (REML, random intercept + slope) | p ≤ 50 |
| GLMM-PQL | `MASS::glmmPQL` | binary, p ≤ 50 |
| GMEXGBoost | from authors' `gmexgb()` using `xgboost` + `lme4::glmer` | binary |

Common boosting budget for MEXGB/MEGB/XGBoost/GBM/GMEXGBoost/GPBoost, (tuned for all
boosting methods): 500 trees, learning rate 0.05, depth 2, min 2 obs per node, exact xgboost trees; EM: max 100 iterations, tolerance 0.001 (package defaults); forests
500 trees. All settings are in `config.R`.
