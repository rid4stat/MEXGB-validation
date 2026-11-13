## =============================================================================
## MEXGB -- Real-data applications, SURVIVAL outcome (clustered time-to-event data)
##
## Data sets (see R/real_survival.R): Retinopathy, Readmission (low-dimensional; cluster =
## patient) and BreastMulti, OvarianMulti (high-dimensional gene expression; cluster = study).
## Protocol: same within-subject design as the simulation (CFG$real_surv$cv_unit = "within"):
## N_REP_REAL x K random splits; every test set holds round(30% x n_i) observations (v0.4)
## of every cluster, capped at n_i - 1. Data sets (v0.4): Insemination (herds), BreastMulti,
## OvarianMulti (studies); Retinopathy and Readmission (< 10 per cluster) are not run. Clusters are
## seen in training, so the frailty models predict with f(x) + z'b_i. Metrics: Harrell's C-index, integrated Brier score.
## Output: results/raw/sreal_<data>.rds (fold metrics + attr 'pred')
##
## Run from the 'analysis' folder:   Rscript 07_real_survival.R
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
res <- run_block_local("sreal", CFG)
all <- rbind_fill(unname(res))
saveRDS(all, file.path(CFG$dir_raw, "sreal_all.rds"))
utils::write.csv(all, file.path(CFG$dir_raw, "sreal_all.csv"), row.names = FALSE)
message("sreal finished: ", nrow(all), " fold x method rows.")
print(write_real_summary(CFG$real_surv$datasets, CFG, file.path(CFG$dir_raw, "sreal_data_summary.csv")))
utils::sessionInfo()
