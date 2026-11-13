## =============================================================================
## MEXGB paper 1 -- Simulation study, CONTINUOUS outcome (paper Sections 6.1-6.3)
##
## Main study : {linear, nonlinear} x p in {6, 170, 2000}
## Sensitivity: nonlinear, p = 170, rho_W in {0, 0.2, 0.6, 0.8} x sigma in {0.5, 1.5}
##              (the main-study pair (0.6, 0.5) is not repeated)
## Each scenario: N_REP simulated data sets, each randomly split K times into a training
## and a test set; each test set = 3 randomly picked measurements of every individual
## (default 10 x 10 = 100 splits). Raw fold-level results -> results/raw/csim_<scenario>.rds
##
## Run from the 'analysis' folder:   Rscript 01_sim_continuous.R
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
res <- run_block_local("csim", CFG)
all <- rbind_fill(unname(res))
saveRDS(all, file.path(CFG$dir_raw, "csim_all.rds"))
utils::write.csv(all, file.path(CFG$dir_raw, "csim_all.csv"), row.names = FALSE)
message("csim finished: ", nrow(all), " fold x method rows.")
utils::sessionInfo()
