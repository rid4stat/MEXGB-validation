## =============================================================================
## MEXGB -- Simulation study, SURVIVAL outcome (clustered time-to-event, frailty Cox model)
##
## lambda_ij(t) = lambda0(t) exp(eta_ij),  eta_ij = kappa * f_std(x_ij) + b0i + b1i t_ij,
## Weibull baseline (shape CFG$surv_sim$weibull_shape), n = 20 clusters x 10 failure times,
## same covariates, fixed effects (linear / nonlinear) and correlated random intercept + slope
## as the continuous/binary designs (frailty SDs in CFG$surv_sim). kappa is calibrated once
## per scenario (MEXGB::calibrateSim) so that the population concordance of the fixed effect
## equals CFG$surv_sim$target_c; independent exponential censoring with rate solved in every
## data set so that the censored fraction equals the target.
## Main study     : {linear, nonlinear} x p in {6, 170, 2000}, censoring 30%
## Censoring study: nonlinear, p = sens_p, censoring in CFG$surv_sim$cens_grid
## Metrics: Harrell's C-index, integrated Brier score, variable-selection TPR/FPR, time.
## Raw fold-level results -> results/raw/ssim_<scenario>.rds
##
## Run from the 'analysis' folder:   Rscript 06_sim_survival.R
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
res <- run_block_local("ssim", CFG)
all <- rbind_fill(unname(res))
saveRDS(all, file.path(CFG$dir_raw, "ssim_all.rds"))
utils::write.csv(all, file.path(CFG$dir_raw, "ssim_all.csv"), row.names = FALSE)
message("ssim finished: ", nrow(all), " fold x method rows.")
utils::sessionInfo()
