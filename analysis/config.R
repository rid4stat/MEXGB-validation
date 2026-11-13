## =============================================================================
## MEXGB paper 1 (continuous + binary) -- global configuration
## Every script sources this file. Change settings HERE only.
## =============================================================================

CFG <- list()

## ---- run size --------------------------------------------------------------
## Design of the published MEGB paper: blocked (subject-level) K-fold CV, repeated
## N_REP times -> K x N_REP = 100 train/test splits per scenario. In the simulation a
## fresh data set is generated for every repetition.
CFG$K          <- 10      # train/test splits per data set (simulation) / folds (real data)
CFG$N_REP      <- 10      # repetitions (simulation): 10 x 10 = 100 splits
CFG$N_REP_REAL <- 10      # repetitions (real data):  10 x 10 = 100 splits
CFG$SEED       <- 2026    # master seed (all other seeds are derived from it)

## Parallel workers (folds of one repetition are run in parallel).
## Uses the 'future' framework (works on Windows, macOS and Linux).
CFG$N_CORES    <- max(1, parallel::detectCores() - 1)

## SMOKE TEST: set to TRUE (or run with env var MEXGB_SMOKE=1) to run a tiny version of
## every script end-to-end in a few minutes (2 folds, 1 repetition, few trees, small p).
CFG$SMOKE <- identical(Sys.getenv("MEXGB_SMOKE"), "1")

## ---- simulation train/test design ------------------------------------------------
## "To evaluate prediction performance, data associated to one simulation are randomly split
##  K times into a training set and a test set. Each test set is obtained by randomly picking
##  three measurements of each individual." (within-subject design, Capitaine et al. 2021)
## Test subjects are therefore seen in training and every mixed model predicts them with its
## predicted random effects. "cluster" restores the previous subject-level hold-out.
CFG$sim_cv_unit        <- "within"
CFG$n_test_per_subject <- 1     # minimum test measurements of a subject (with n_i >= 2)
## Test measurements of subject i = round(test_frac * n_i), at least n_test_per_subject and at
## most n_i - 1 (every subject keeps >= 1 training measurement): 30% of every subject's
## measurements, i.e. 6 of the 20 occasions in the simulation (14 for training).
CFG$test_frac          <- 0.30

## ---- real-data train/test design -----------------------------------------------------------
## Same within-subject design as the simulation: K random splits per repetition, each test set =
## max(15%, 3) measurements of every subject (capped at n_i - 1), and every mixed model predicts the
## (seen) test subjects with f(x) + z'b_i. "cluster" restores subject-level CV (unseen subjects).
CFG$real_cv_unit       <- "within"

## ---- simulation design (paper Section 6.1) ----------------------------------
CFG$sim <- list(
  n = 50, time_points = 20,          # 50 subjects x 20 occasions (v0.4; was 20 x 10)
  p_grid = c(6, 170, 2000),
  dgp = c("linear", "nonlinear"),
  rho_W = 0.6, rho_Z = 0.6,
  tau0 = sqrt(0.5), tau1 = sqrt(3), sigma = 0.5,
  ## sensitivity sub-study (continuous, nonlinear, p = 170): all (rho_W, sigma) pairs
  ## except the main-study pair (0.6, 0.5); rho_W = 0 gives independent predictors and
  ## independent within-subject errors (correlation only through the random effects)
  sens_p = 170, sens_rho = c(0, 0.2, 0.6, 0.8), sens_sigma = c(0.5, 1.5),
  run_sensitivity_binary = FALSE,
  n_active = 6,          # active predictors = Var1..Var6
  ## BINARY outcome: the signal strength kappa and intercept c0 of
  ##   eta = c0 + kappa * f_std + b0 + b1 t
  ## are calibrated (MEXGB::calibrateSim) so that the population (fixed-effect) AUC equals
  ## target_auc and the event rate equals prevalence. Set target_auc = NULL for the
  ## original design (kappa = 1, prevalence 0.5).
  target_auc = 0.85, prevalence = 0.5,
  ## optional extra binary scenarios varying the prevalence (nonlinear, p = 170)
  bin_prev_grid = c(0.2, 0.3),
  pilot_n = 300          # pilot individuals used by the calibration
)

## ---- survival simulation design (frailty Cox model) ----------------------------
## eta = kappa * f_std + b0 + b1 t, Weibull baseline (shape 1.5), independent exponential
## censoring; kappa calibrated so that the population concordance equals target_c.
CFG$surv_sim <- list(
  n = 50, time_points = 20,          # 50 clusters x 20 correlated failure times (v0.4)
  p_grid = c(6, 170, 2000),
  dgp = c("linear", "nonlinear"),
  rho_W = 0.6, rho_Z = 0.6,
  tau0 = sqrt(0.5), tau1 = sqrt(0.5),  # frailty SDs (random intercept, random slope)
  target_c = 0.75, censoring = 0.3,
  cens_grid = c(0.5), sens_p = 170,  # extra censoring levels (nonlinear, p = sens_p)
  weibull_shape = 1.5, weibull_scale = 1,
  pilot_n = 300
)

## ---- methods ------------------------------------------------------------------
CFG$methods_cont <- c("MEXGB", "XGBoost", "MEGB", "MERF", "REEMforest", "RF", "GBM",
                      "GPBoost", "mboost", "glmmLasso", "LMM")
CFG$methods_bin  <- c("MEXGB", "XGBoost", "GMEXGBoost", "GLMM-PQL", "glmmLasso", "GBM", "RF",
                      "mboost", "GPBoost")
CFG$methods_surv <- c("MEXGB", "XGBCox", "GBM", "RSF", "mboost", "CoxLasso", "CoxFrailty")

## Parametric mixed models are fitted only when they admit a stable fit:
CFG$p_max_lmm       <- 50    # LMM / GLMM-PQL / CoxFrailty only if p <= 50
CFG$p_max_glmmlasso <- 500   # glmmLasso only if p <= 500 and p < N_train

## ---- hyper-parameters (held fixed across folds and scenarios) ---------------
## One common boosting budget for every boosting method (MEXGB, MEGB, GBM,
## GMEXGBoost, GPBoost) so that differences reflect the estimation scheme.
## v0.4 common budget (tuned for ALL boosting methods alike, see CHANGES_v0.4.md): depth 2 (at
## depth 2 gbm and xgboost grow the same 4-leaf trees; xgboost max_depth = 3 allows 8 leaves
## against gbm's 4) and >= 2 observations per leaf (5 was too restrictive for every booster).
CFG$boost <- list(ntree = 500, shrinkage = 0.05, depth = 2, minobs = 2,
                  bag.fraction = 1)          # gbm: bag.fraction = 1 as in MEGB
## XGBoost-specific regularisation (MEXGB, GMEXGBoost)
CFG$xgb <- list(lambda = 1, alpha = 0, gamma = 0, subsample = 1,
                tree_method = "auto")   # "auto": exact greedy (as gbm) up to 10000 rows (histogram binning
                                        # cost xgboost ~1 MSE unit at n = 50 x 20), hist above
## Column subsampling rate for MEXGB / GMEXGBoost (xgboost) as a function of p
## (was 0.5 for p >= 100 and 0.15 for p >= 1000: with 6 active predictors most trees then never
##  see the signal, which made XGBoost/MEXGB clearly worse than gbm at p = 170 and 2000)
CFG$colsample <- function(p) 1
## EM settings shared by MEXGB, MEGB, MERF and REEMforest (package defaults)
CFG$em <- list(iter = 100, delta = 0.001)
## MEXGB / MEGB 0.3 estimation fixes (see CHANGES_v0.3.md). crossfit = 0 and init = "zero"
## reproduce the 0.2 algorithms (in-sample E-step, which lets the boosting absorb the random
## effects so that MEXGB / MEGB collapse to plain XGBoost / GBM).
CFG$em$crossfit <- 5          # folds of the cross-fitted (out-of-fold) E-step
CFG$em$init     <- "lmm"      # EM starting values from the LMM with a constant fixed effect
CFG$em$inner    <- 50         # E-step (BLUP + variance components [+ PQL linearisation]) iterated
                              # to convergence for each boosting refit (v0.4; 1 = MEXGB 0.3)
## MEXGB 0.5 (see CHANGES_v0.5.md): fixed effect = penalised linear part (lasso, lambda chosen
## on the cross-fitting folds) + trees whose number is chosen on the same folds (out-of-fold
## loss, 0 trees allowed, growth stopped 200 rounds after the best); comparable Laplace
## marginal-likelihood criterion for the binary and Cox EM.
CFG$mexgb <- list(linear = "auto", init = "auto", honest = TRUE, honest_d = TRUE, lin_alpha = 1, ntree_select = "cv", ntree_max = 5000, shrinkage = 0.1)
## Random forests (MERF, REEMforest, RF)
CFG$forest <- list(ntree = 500,
                   mtry_reg = function(p) max(1, floor(p / 3)),
                   mtry_cls = function(p) max(1, floor(sqrt(p))))
## mboost (glmboost, linear base-learners, centred covariates)
CFG$mboost <- list(mstop = 300, nu = 0.1)
## glmmLasso
CFG$glmmlasso <- list(lambda = 10)
## GMEXGBoost (Asadi et al. 2026; port of the authors' gmexgb() function)
CFG$gmex <- list(itmax = 30, toll = 0.02, random_slope = TRUE,
                 glm_init_max_p = 15,   # p > 15: ridge-penalised GLM initialisation (glmnet)
                 eta_bound = 10)        # bound on the initial GLM linear predictor
## GPBoost
## (feature_fraction < 1 destabilises GPBoost's joint covariance estimation -> package default 1)
CFG$gpb <- list(random_slope = TRUE, feature_fraction = 1)
## Survival competitors
CFG$rsf <- list(nodesize = 15, importance = "random")   # randomForestSRC (VIMP: random-split)
CFG$coxlasso <- list(nfolds = 5, s = "lambda.min",     # glmnet::cv.glmnet
                     nlambda = 50, lambda.min.ratio = 0.05)
CFG$coxfrailty <- list(random_slope = TRUE)             # coxme: (1 + time | id) when Z has a slope

## ---- real data ----------------------------------------------------------------
CFG$real <- list(
  ## v0.4: continuous and binary outcomes use different data sets; every data set has >= 10
  ## measurements per subject (PBC: patients with >= 10 visits)
  datasets_cont = c("BodyWeight", "Soybean", "Sprockett"),
  datasets_bin  = c("PBC", "OvarianDebulking"),
  ## outcome transformations (continuous): Soybean leaf weight is strongly right-skewed
  ## (mean 6.4, SD 6.9, range 0.03-30) -> log scale
  transform = c(BodyWeight = "identity", Soybean = "log", Sprockett = "identity"),
  ## PBC (JM::pbc2, Mayo Clinic primary biliary cirrhosis trial): patients with >= min_visits
  ## visits; binary outcome = hepatomegaly at the visit; low-dimensional clinical predictors
  pbc = list(file = "data/pbc2.rda", min_visits = 10, outcome = "hepatomegaly"),
  ## Sprockett microbiome cohort (paper Section 7), read from data/SprockettData.csv:
  ## feces samples of subjects with >= min_visits visits (193 subjects, 1815 samples);
  ## features = the ASV columns (log(1 + relative abundance)), centred within sample, the
  ## lowest-variance decile of ASVs removed; outcomes = the 'shannon' and 'shannon_binary'
  ## columns of the file; longitudinal time:
  ## time_var = "followup" (default): months since the first sample (Age_Months), rescaled
  ## to [0, 1]; time_var = "visit": visit index ordered by age, rescaled to [0, 1] (paper v1).
  sprockett = list(file = "data/SprockettData.csv", min_visits = 3, sample_type = "Feces",
                   drop_low_var_frac = 0.10, centre = TRUE, time_var = "followup"),
  ## forests on the p = 1880 cohort are expensive; reduce if needed
  forest_ntree_hd = 300
)

CFG$real$datasets <- unique(c(CFG$real$datasets_cont, CFG$real$datasets_bin))

## ---- real SURVIVAL data -------------------------------------------------------------
## Low-dimensional, cluster = patient (subject-level CV, unseen patients):
##   Retinopathy  survival::retinopathy  -- 197 patients x 2 eyes, time to blindness
##   Readmission  frailtypack::readmission -- 403 colorectal-cancer patients, 861 recurrent
##                hospital readmission gap times
## High-dimensional, cluster = study (observation-level CV stratified by study, seen studies):
##   BreastMulti  Bioconductor breastCancerVDX/TRANSBIG/MAINZ/UNT (Affymetrix hgu133a),
##                distant-metastasis-free survival, 4 cohorts, ~860 patients
##   OvarianMulti Bioconductor curatedOvarianData, overall survival, one cluster per study
## v0.4: every survival data set has >= 10 observations per cluster. Retinopathy (2 eyes per
## patient) and Readmission (median 2 gap times) are kept as loaders but not run by default.
##   Insemination parfm::insem -- 10513 dairy cows in 181 herds (cluster = herd, herds with
##                >= min_cows cows), time from calving to first insemination; low-dimensional
CFG$real_surv <- list(
  datasets = c("Insemination", "BreastMulti", "OvarianMulti"),
  insem = list(file = "data/insem.RData", min_cows = 10),
  ## same within-subject design as the simulation (max(15%, 3) measurements of every cluster,
  ## capped at n_i - 1; seen clusters). Retinopathy: 1 of the 2 eyes of every patient;
  ## Readmission: patients with a single gap time stay in training; multi-study sets: 15% of
  ## every study. "cluster" / "observation" restore the previous designs.
  cv_unit = c(Insemination = "within", Retinopathy = "within", Readmission = "within",
              BreastMulti = "within", OvarianMulti = "within"),
  n_genes = 2000,                 # most variable genes kept in the multi-study sets
  max_years = 10,                 # administrative censoring of follow-up (multi-study sets)
  breast_cohorts = c("VDX", "TRANSBIG", "MAINZ", "UNT"),
  breast_local_dir = NULL,        # optional folder with vdx.rda, transbig.rda, ... (no Bioconductor)
  ovarian_local_dir = NULL,       # optional folder with curatedOvarianData *_eset.rda files
  ovarian_min_n = 40, ovarian_min_events = 15,
  readmission_random_slope = FALSE  # TRUE: random slope on the (rescaled) episode number
)

## ---- output -------------------------------------------------------------------
CFG$dir_raw     <- "results/raw"
CFG$dir_tables  <- "results/tables"
CFG$dir_figures <- "results/figures"
CFG$dir_data    <- "data"

## ---- smoke-test overrides --------------------------------------------------------
if (isTRUE(CFG$SMOKE)) {
  CFG$K <- 2; CFG$N_REP <- 1; CFG$N_REP_REAL <- 1
  CFG$sim$p_grid <- c(6, 30)
  CFG$sim$sens_p <- 30; CFG$sim$sens_rho <- c(0, 0.6); CFG$sim$sens_sigma <- c(0.5)
  CFG$boost$ntree <- 30; CFG$forest$ntree <- 30; CFG$em$iter <- 5; CFG$em$crossfit <- 3
  CFG$mexgb$ntree_max <- 60
  CFG$mboost$mstop <- 30; CFG$gmex$itmax <- 3
  CFG$sim$bin_prev_grid <- c(0.3); CFG$sim$pilot_n <- 100
  CFG$surv_sim$p_grid <- c(6, 30); CFG$surv_sim$cens_grid <- c(0.5); CFG$surv_sim$sens_p <- 30; CFG$surv_sim$pilot_n <- 100
  CFG$real_surv$n_genes <- 200
  CFG$real$forest_ntree_hd <- 20
  CFG$dir_raw <- "results_smoke/raw"; CFG$dir_tables <- "results_smoke/tables"
  CFG$dir_figures <- "results_smoke/figures"
}
## ---- environment overrides (used by the Hali campaign; see hali/config.sh) -----------
## MEXGB_OUT       results root: <root>/raw, <root>/tables, <root>/figures
## MEXGB_N_REP     repetitions per simulation scenario
## MEXGB_N_REP_REAL repetitions per real data set
## MEXGB_K         folds
## MEXGB_N_CORES   parallel workers for the LOCAL driver (array tasks always use 1)
env_int <- function(v) { x <- Sys.getenv(v, ""); if (nzchar(x)) as.integer(x) else NULL }
if (!is.null(env_int("MEXGB_N_REP")))      CFG$N_REP      <- env_int("MEXGB_N_REP")
if (!is.null(env_int("MEXGB_N_REP_REAL"))) CFG$N_REP_REAL <- env_int("MEXGB_N_REP_REAL")
if (!is.null(env_int("MEXGB_K")))          CFG$K          <- env_int("MEXGB_K")
if (!is.null(env_int("MEXGB_N_CORES")))    CFG$N_CORES    <- env_int("MEXGB_N_CORES")
if (nzchar(Sys.getenv("MEXGB_OUT"))) {
  CFG$dir_raw     <- file.path(Sys.getenv("MEXGB_OUT"), "raw")
  CFG$dir_tables  <- file.path(Sys.getenv("MEXGB_OUT"), "tables")
  CFG$dir_figures <- file.path(Sys.getenv("MEXGB_OUT"), "figures")
}

for (d in c(CFG$dir_raw, CFG$dir_tables, CFG$dir_figures, CFG$dir_data))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
