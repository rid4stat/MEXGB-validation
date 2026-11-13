## =============================================================================
## MEXGB paper 1 -- Real-data applications (paper Section 7)
##
## Data sets (v0.4): continuous BodyWeight, Soybean (log), Sprockett; binary PBC, OvarianDebulking.
## Outcomes : continuous for the first three, binary (hepatomegaly, suboptimal debulking)
##            for the last two.
## Protocol : same within-subject design as the simulation (CFG$real_cv_unit = "within"):
##            N_REP_REAL x K random train/test splits (default 10 x 10 = 100); every test set
##            holds round(30% x n_i) randomly chosen measurements of EVERY subject (v0.4)
##            (capped at n_i - 1), so test subjects are seen in training and every mixed
##            model predicts them with f(x) + z'b_i. CFG$real_cv_unit = "cluster" restores
##            subject-level CV (unseen subjects, random effects = 0).
## Output   : results/raw/real_<data>_<outcome>.rds (fold metrics + attr 'pred')
##
## Run from the 'analysis' folder:   Rscript 03_real_data.R
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
res <- run_block_local("real", CFG)
all <- rbind_fill(unname(res))
saveRDS(all, file.path(CFG$dir_raw, "real_all.rds"))
utils::write.csv(all, file.path(CFG$dir_raw, "real_all.csv"), row.names = FALSE)
message("real finished: ", nrow(all), " fold x method rows.")
print(write_real_summary(real_datasets(CFG), CFG, file.path(CFG$dir_raw, "real_data_summary.csv")))
utils::sessionInfo()
