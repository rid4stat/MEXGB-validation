#!/usr/bin/env Rscript
# =============================================================================
# hali/R/00_prep.R -- registry (with calibrations) + real data, computed ONCE
# =============================================================================
# * builds the study registry with the binary / survival calibrations (target AUC
#   and prevalence; target C-index) -- done once here instead of in every task;
# * loads every real data set once and saves it to prep/real/<name>.rds, so the
#   array tasks never re-read the Sprockett CSV or the Bioconductor sets;
# * freezes the configuration (prep/cfg.rds) that every task then uses.
# =============================================================================
source(file.path(Sys.getenv("HALI_ROOT", "../.."), "hali", "R", "hali_env.R"))
t0 <- Sys.time()
log_msg("PREP | run %s | blocks %s", Sys.getenv("RUN_LABEL"), paste(BLOCKS, collapse = ","))

tt_file <- file.path(PREP_DIR, "tasks.rds")
if (!file.exists(tt_file)) stop("no task table: run './submit.sh plan' first", call. = FALSE)
tt <- readRDS(tt_file)

reg <- build_registry(CFG, BLOCKS, calibrate = TRUE)
if (!all(unique(tt$study_id) %in% names(reg)))
  stop("task table and registry disagree (config changed since plan?): re-run './submit.sh plan'",
       call. = FALSE)
for (st in reg) if (!is.null(st$data$cal))
  log_msg("calibration %-40s kappa = %.3f  achieved = %.3f", st$id, st$data$cal$signal,
          st$data$cal$achieved)
saveRDS(reg, file.path(PREP_DIR, "registry.rds"))
saveRDS(CFG, file.path(PREP_DIR, "cfg.rds"))
write_calibration_tables(reg, CFG$dir_raw)

real_needed <- unique(unlist(lapply(reg, function(st) if (st$data$type == "real") st$data$dataset)))
status <- list()
for (ds in real_needed) {
  f <- file.path(REAL_DIR, paste0(ds, ".rds"))
  ok <- tryCatch({
    d <- load_any_real(ds, CFG)
    atomic_save(d, f)
    log_msg("real data %-14s %s", ds, attr(d, "description") %||% "")
    TRUE
  }, error = function(e) { log_msg("real data %-14s UNAVAILABLE: %s", ds, conditionMessage(e)); FALSE })
  status[[ds]] <- data.frame(dataset = ds, available = ok)
}
if (length(status)) utils::write.csv(do.call(rbind, status), file.path(PREP_DIR, "real_status.csv"),
                                     row.names = FALSE)
if ("real" %in% BLOCKS)
  write_real_summary(real_datasets(CFG), CFG, file.path(CFG$dir_raw, "real_data_summary.csv"), REAL_DIR)
if ("sreal" %in% BLOCKS)
  write_real_summary(CFG$real_surv$datasets, CFG, file.path(CFG$dir_raw, "sreal_data_summary.csv"), REAL_DIR)

log_msg("PREP done in %.1f min | %d studies | %d tasks | peak RSS %.2f GB",
        as.numeric(difftime(Sys.time(), t0, units = "mins")), length(reg), nrow(tt), peak_rss_gb())
