#!/usr/bin/env Rscript
# =============================================================================
# hali/R/01_plan.R -- preflight + task table + resource projection
# =============================================================================
# Runs in the foreground (seconds). Checks the R environment, builds the task
# table from analysis/config.R + R/registry.R (no data needed), and projects the
# campaign's cost. Writes prep/tasks.rds and prep/tab_plan_tasks.csv, which size
# the SLURM array. Re-run it whenever BLOCKS, RUN_LABEL or the run size change.
# =============================================================================
source(file.path(Sys.getenv("HALI_ROOT", "../.."), "hali", "R", "hali_env.R"))

cat("\n== ENVIRONMENT ================================================================\n")
cat(sprintf("  %-28s %s\n", "R", R.version.string))
cat(sprintf("  %-28s %s\n", "execution environment",
            if (in_container()) "Apptainer container" else "native (host R)"))
req <- c("MEXGB", "xgboost", "MEGB", "gbm", "randomForest", "LongituRF", "gpboost",
         "mboost", "glmmLasso", "nlme", "MASS", "lme4", "survival", "randomForestSRC",
         "glmnet", "coxme", "ggplot2")
opt <- c("Biobase", "breastCancerVDX", "breastCancerTRANSBIG", "breastCancerMAINZ",
         "breastCancerUNT", "curatedOvarianData")
for (p in c(req, opt)) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA)
  cat(sprintf("  %-28s %s\n", p, if (is.na(v)) {
    if (p %in% req) "MISSING  <-- required" else "missing (multi-study survival data skipped)"
  } else paste("OK", v)))
}
for (f in c(CFG$real$sprockett$file, file.path(CFG$dir_data, "readmission.rda")))
  cat(sprintf("  %-28s %s\n", basename(f), if (file.exists(f)) "OK" else "MISSING"))
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]

cat("\n== RUN ========================================================================\n")
cat(sprintf("  run label %s | output %s\n", Sys.getenv("RUN_LABEL"), OUT_DIR))
cat(sprintf("  blocks %s | K = %d folds | N_REP = %d (simulation), %d (real) | smoke = %s\n",
            paste(BLOCKS, collapse = ","), CFG$K, CFG$N_REP, CFG$N_REP_REAL, isTRUE(CFG$SMOKE)))

reg <- build_registry(CFG, BLOCKS, calibrate = FALSE)
tt <- build_task_table(reg, CFG, SPLIT_THRESHOLD)
saveRDS(tt, file.path(PREP_DIR, "tasks.rds"))
utils::write.csv(tt, file.path(PREP_DIR, "tab_plan_tasks.csv"), row.names = FALSE)

cat("\n== TASK ALGEBRA ===============================================================\n")
by_block <- do.call(rbind, lapply(split(tt, tt$block), function(x) data.frame(
  block = x$block[1], studies = length(unique(x$study_id)), tasks = nrow(x),
  method_split_studies = length(unique(x$study_id[x$unit != "all"])),
  est_core_hours = round(sum(x$est_sec) / 3600, 1), longest_task_min = round(max(x$est_sec) / 60, 1))))
print(by_block, row.names = FALSE)
cat(sprintf("\n  TOTAL: %d studies, %d tasks, ~%.0f core-hours (estimate)\n",
            length(reg), nrow(tt), sum(tt$est_sec) / 3600))

cat("\n== RESOURCE PROJECTION ========================================================\n")
conc <- as.integer(Sys.getenv("ARRAY_CONCURRENCY", "200"))
wall <- max(sum(tt$est_sec) / conc, max(tt$est_sec))
cat(sprintf("  with %d concurrent tasks: ~%.1f h wall-clock for the array (+ queueing)\n",
            conc, wall / 3600))
cat(sprintf("  longest single task: ~%.0f min (%s) -- FIT_TIME is %s\n",
            max(tt$est_sec) / 60, tt$study_id[1], Sys.getenv("FIT_TIME", "?")))
cat("  These are ESTIMATES from single-core timings at N = 180. The rehearsal\n",
    " measures every task: results/tables/tab_task_timing.csv after combine.\n")
if (length(miss)) {
  cat("\n!! required packages missing: ", paste(miss, collapse = ", "),
      "\n!! fit tasks will fail for those methods. Rebuild the image (hali/container).\n")
}
cat("\nplan written: ", file.path(PREP_DIR, "tab_plan_tasks.csv"), "\n")
