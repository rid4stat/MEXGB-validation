## =============================================================================
## MEXGB -- run the complete analysis (continuous, binary and survival)
## From the 'analysis' folder:  Rscript run_all.R
## Smoke test (few minutes):    MEXGB_SMOKE=1 Rscript run_all.R
##   (on Windows PowerShell:    $env:MEXGB_SMOKE="1"; Rscript run_all.R)
## Each step caches its results in results/raw, so an interrupted run resumes.
## =============================================================================
t0 <- Sys.time()
for (s in c("01_sim_continuous.R", "02_sim_binary.R", "03_real_data.R",
            "06_sim_survival.R", "07_real_survival.R",
            "05_convergence.R", "08_data_description.R", "04_tables_figures.R")) {
  message("\n########## ", s, " ##########")
  local(source(s, local = TRUE))
}
message("All done in ", format(round(difftime(Sys.time(), t0, units = "hours"), 2)))
