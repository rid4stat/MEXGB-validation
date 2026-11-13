## =============================================================================
## MEXGB paper 1 -- Simulation study, BINARY outcome (paper Sections 6.1 and 6.4)
##
## eta_ij = c0 + kappa * f_std(x_ij) + b0i + b1i t_ij,  y_ij ~ Bernoulli(expit(eta_ij)).
## kappa (signal strength) and c0 (intercept) are calibrated ONCE per scenario on a large
## pilot sample (MEXGB::calibrateSim) so that
##   * the population (fixed-effect) AUC equals CFG$sim$target_auc, and
##   * the event rate equals the imposed prevalence,
## and are then held fixed across the N_REP replicate data sets.
## Main study      : {linear, nonlinear} x p in {6, 170, 2000}, prevalence CFG$sim$prevalence
## Prevalence study: nonlinear, p = sens_p, prevalence in CFG$sim$bin_prev_grid
## Each scenario: N_REP simulated data sets, each randomly split K times into a training
## and a test set; each test set = 3 randomly picked measurements of every individual.
## Raw fold-level results -> results/raw/bsim_<scenario>.rds
## Calibration table      -> results/raw/bsim_calibration.csv
##
## Run from the 'analysis' folder:   Rscript 02_sim_binary.R
## =============================================================================

source("config.R")
source("R/methods.R")
source("R/cv.R")
source("R/real_data.R")
source("R/real_survival.R")
source("R/registry.R")
suppressPackageStartupMessages(library(MEXGB))
setup_parallel(CFG)

## all scenarios / data sets of this block are defined once in R/registry.R
res <- run_block_local("bsim", CFG)
all <- rbind_fill(unname(res))
saveRDS(all, file.path(CFG$dir_raw, "bsim_all.rds"))
utils::write.csv(all, file.path(CFG$dir_raw, "bsim_all.csv"), row.names = FALSE)
message("bsim finished: ", nrow(all), " fold x method rows.")
utils::sessionInfo()
