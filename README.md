# MEXGB — package and reproducible analysis (continuous, binary and survival outcomes)

```
MEXGB/                 R package source (xgboost back-end); version 0.3.1
MEGB/                  MEGB 0.3 source (same cross-fitted E-step); MEGB_0.4.tar.gz built package
paper_text/            real_data_description.tex (data-description subsection + table)
diagnostics/           before/after benchmark of the 0.3 fixes + oracle decomposition
MEXGB_0.4.0.tar.gz     built source package (install.packages("MEXGB_0.4.0.tar.gz", repos = NULL, type = "source"))
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

`hali/` runs the same analysis as a SLURM campaign inside an Apptainer image, in the
structure of the mmCTBN campaign: `hali/config.sh` (the only file to edit),
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
Workers: `CFG$N_CORES` (default = cores − 1). **Send back `results/results_summary.txt`,
`results/tables/` and `results/figures/`** and I will update the paper and thesis.

### Version 0.4 (see CHANGES_v0.4.md)
* Design n = 50 × t = 20 with a 30% within-subject test set, for the simulation and all real data.
* MEXGB/MEGB: the E-step is iterated to convergence after each boosting refit (`inner`). Cox:
  frailties centred (identifiability) and the boosted fixed effect recalibrated by its out-of-fold
  calibration slope (`calibrate`).
* Common boosting budget re-tuned for all boosters: depth 2, ≥ 2 observations per leaf, exact trees.
* New real data sets: PBC, OvarianDebulking (binary), Insemination (survival).

### Update 0.3.1
* **Survival fix.** With the cross-fitted E-step, the Cox EM could diverge at high p. The Cox
  PQL working model now estimates a dispersion phi, as `MASS::glmmPQL` does
  (`dispersion = "auto"`, which is multiplicative for Cox and phi = 1 for binary; see
  CHANGES_v0.3.md section 6). All families stop and return the best iteration if the
  likelihood becomes non-finite.
* Real data (continuous, binary, survival) use the simulation's within-subject partition: the
  test set is max(15%, 3) measurements of every subject, and seen subjects are predicted with
  their random effects.
* `05_convergence.R` writes one log-likelihood figure per outcome
  (`fig_mexgb_loglik_continuous.pdf`, `_binary.pdf`, `_survival.pdf`). It uses the 0.3
  estimation settings, caches every fit in `results/raw/loglik/`, skips a failing fit instead
  of stopping, and runs in parallel on the Hali combine job's CPUs.

### What is new in version 0.3 (see CHANGES_v0.3.md)
* **Cross-fitted E-step in MEXGB and MEGB.** The 0.2 E-step computed the BLUPs, sigma^2 and D
  from *in-sample* boosting residuals. Those are almost zero because 500 trees interpolate the
  training data, so the random effects collapsed to 0 and MEXGB/MEGB reduced to XGBoost/GBM.
  The E-step now uses out-of-fold predictions of f (`crossfit = 5`), and the EM starts from
  the LMM (`init = "lmm"`).
* **`colsample = 1`.** The old 0.5 / 0.15 column subsampling hid the 6 active predictors from
  most trees.
* **MEGB 0.3** (folder `MEGB/`, `MEGB_0.4.tar.gz`) carries the same fix and replaces CRAN's
  0.2. It also returns the best iteration instead of `NULL` and matches random effects by id
  in `predict()`.
* `diagnostics/` holds the before/after benchmark and the oracle decomposition.

### What was new in version 0.2
* **Survival outcomes.** `MEXGB(..., family = "cox")` fits a mixed-effects (log-normal frailty)
  Cox model: the M-step is `xgboost` with `objective = "survival:cox"` and the random effects
  as offset; the E-step is the PQL BLUP of the Poisson representation of the Cox likelihood
  (Breslow baseline). `predict(..., type = "link" | "response" | "survival", times =)`.
* **Binary simulation with a target AUC.** `calibrateSim()` finds the signal multiplier κ and
  intercept c0 so that the population (fixed-effect) AUC equals `target_auc` and the event rate
  equals the imposed `prevalence`; `simLong(..., calibration =)` then keeps the DGP fixed across
  replicates. Default: AUC 0.85, prevalence 0.5, plus a prevalence study (0.2, 0.3).
* **Survival simulation** with the same covariates, fixed effects and random intercept+slope,
  Weibull baseline, κ calibrated to a target C-index (0.75), exponential censoring solved to a
  target fraction (30%; 50% in a censoring study).
* **EM safeguard** (`patience = 5`): with column subsampling the likelihood jitters around its
  limit; the EM also stops after 5 iterations without improvement and returns the
  best-likelihood state (Gaussian fits that converge normally are unchanged).
* **Sprockett** is read from `data/SprockettData.csv` (feces, ≥3 visits → 193 subjects, 1815
  samples; ASVs centred within sample, lowest-variance decile removed → p = 1879; outcomes =
  the file's `shannon` and `shannon_binary`). The random-slope time is the exact follow-up
  (months since first sample from `Age_Months`); `time_var = "visit"` restores the paper's
  visit index.

### Survival data and methods
| Data set | Source | Clusters | N | p | CV |
|---|---|---|---|---|---|
| Insemination | `parfm::insem` (shipped in data/) | 174 herds | 10476 | 4 | within-cluster 30% |
| BreastMulti | Bioconductor breastCancerVDX/TRANSBIG/MAINZ/UNT (DMFS) | 4 cohorts | 867 | 2000 of 22283 probes | stratified by cohort |
| OvarianMulti | Bioconductor curatedOvarianData (OS) | all eligible studies | ~2000+ | 2000 common genes | stratified by study |

Survival competitors: XGB-Cox (`xgboost`, no frailty), GBM (`gbm`, coxph), RSF
(`randomForestSRC`), mboost (`glmboost`, CoxPH), Cox-Lasso (`glmnet`), CoxFrailty (`coxme`,
p ≤ 50). Metrics: Harrell's C-index, integrated Brier score (IPCW; Breslow-calibrated survival
curves for every method), variable-selection TPR/FPR, time.

### Expected run time (full design, 100 splits per scenario)
Per split at p = 2000 (one core), version 0.3 (cross-fitted E-step, colsample = 1, measured):
MEXGB ≈ 280–300 s, MEGB ≈ 215–235 s, MERF ≈ 160 s, REEMforest ≈ 130 s, RF ≈ 13 s,
XGBoost ≈ 8 s, GBM ≈ 6 s, GPBoost ≈ 5 s, mboost ≈ 3 s. In 0.2, MEXGB took ≈ 45–75 s and MEGB
≈ 75–160 s. The cross-fitting multiplies the boosting work by about (V + 1) = 6. The
continuous/binary study is now roughly 100–120 core-hours, dominated by the p = 2000
continuous cells and Sprockett. Use the Hali campaign, and the rehearsal to measure
Sprockett. `CFG$em$crossfit <- 3` cuts the MEXGB/MEGB time by about a third. Survival adds about 15–20 core-hours (MEXGB-Cox
≈ 40–60 s per split at p = 2000, the other survival methods < 10 s). To cut time,
lower `CFG$N_REP` (e.g. 5 → 50 splits) or `CFG$real$forest_ntree_hd`.

## Design (as in paper/thesis and the published MEGB paper)
* n = 50 subjects × 20 occasions (v0.4; was 20 × 10), t rescaled to [0,1]; x ~ N_p(0, Σ_W), (Σ_W)_gg' = 0.6^|g−g'|;
  p ∈ {6, 170, 2000}, 6 active predictors; linear f = 2x1−2x2+1.5x3−1.5x4+x5−x6, nonlinear
  f = 2 sin(πx1x2)+4(x3−0.5)²+2x4+exp(x5)+x6x1; correlated random intercept/slope
  (τ0² = 0.5, τ1² = 3, ρ_Z = 0.6); AR(1) within-subject errors σ = 0.5 (lag in visit number).
  Binary: f centred/scaled to SD 1, intercept solved so mean P(y=1) = 0.5. → `MEXGB::simLong()`.
* Sensitivity (continuous, nonlinear, p = 170): ρ_W ∈ {0, 0.2, 0.6, 0.8} × σ ∈ {0.5, 1.5}.
* Evaluation (v0.4): every data set, simulated or real, is randomly split K times into a
  training and a test set. Each test set holds 30% of the measurements of EVERY subject
  (round(0.3 × n_i), at least 1, at most n_i − 1; simulation: 6 of 20 for testing, 14 for
  training). 10 simulated data sets × K = 10 splits = 100 splits per scenario; real data 10 × 10.
  Test subjects are seen in training, so every mixed model predicts with f(x) + z'b̂_i.
  Settings: `CFG$test_frac = 0.3`, `CFG$sim_cv_unit`, `CFG$real_cv_unit`, `CFG$real_surv$cv_unit`.
* Real data (v0.4), all with ≥ 10 measurements per subject or cluster: continuous BodyWeight,
  Soybean (log outcome), Sprockett; binary PBC (hepatomegaly, low-dim) and OvarianDebulking
  (high-dim); survival Insemination (dairy herds, low-dim), BreastMulti and OvarianMulti
  (high-dim). See CHANGES_v0.4.md and paper_text/.
* Figures: one panel per scenario / data set with its own y-axis scale; within each panel the
  methods are ordered by mean performance (best on the left).
  Metrics: MSE; AUC, Brier, classification TPR/FPR (0.5); variable-selection TPR/FPR of the
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
| glmmLasso | `glmmLasso::glmmLasso` (λ = 10, random intercept) | p < N_train |
| LMM | `nlme::lme` (REML, random intercept) | p ≤ 50 |
| GLMM-PQL | `MASS::glmmPQL` | binary, p ≤ 50 |
| GMEXGBoost | port of authors' `gmexgb()` using `xgboost` + `lme4::glmer` | binary |

Common boosting budget for MEXGB/MEGB/XGBoost/GBM/GMEXGBoost/GPBoost (v0.4, tuned for all
boosters alike): 500 trees, learning rate 0.05, depth 2, min 2 obs per node, exact xgboost trees; EM: max 100 iterations, tolerance 0.001 (package defaults); forests
500 trees. All settings are in `config.R`.
