#!/usr/bin/env Rscript
# =============================================================================
# hali/R/03_combine.R -- assemble the tasks, then tables and figures
# =============================================================================
# 1. task status + measured timing (tab_task_status.csv, tab_task_timing.csv);
# 2. for every study, binds its task files into results/raw/<study>.rds -- the
#    SAME file a local run writes (fold rows + attr 'pred'), so
# 3. analysis/04_tables_figures.R (and 05_convergence.R) run unchanged on it.
# Missing tasks are reported with a ready-made resubmission command.
# =============================================================================
source(file.path(Sys.getenv("HALI_ROOT", "../.."), "hali", "R", "hali_env.R"))
t0 <- Sys.time()
tt  <- readRDS(file.path(PREP_DIR, "tasks.rds"))
reg <- readRDS(file.path(PREP_DIR, "registry.rds"))
log_msg("COMBINE | %d tasks | %d studies", nrow(tt), length(reg))

## ---- 1. status and timing -------------------------------------------------------------
stems <- task_stem(FITS_DIR, tt$study_id, tt$rep, tt$fold, tt$unit)
tt$done <- file.exists(paste0(stems, ".done")) & file.exists(paste0(stems, ".rds"))
tt$failed <- file.exists(paste0(stems, ".fail"))
TAB <- CFG$dir_tables
utils::write.csv(tt, file.path(TAB, "tab_task_status.csv"), row.names = FALSE)
timing <- do.call(rbind, lapply(which(tt$done), function(i) {
  x <- readRDS(paste0(stems[i], ".rds"))
  data.frame(task_id = tt$task_id[i], study_id = tt$study_id[i], unit = tt$unit[i],
             rep = tt$rep[i], fold = tt$fold[i], est_sec = tt$est_sec[i],
             seconds = round(x$seconds, 1), peak_rss_gb = round(x$peak_rss_gb, 3), host = x$host)
}))
if (!is.null(timing)) {
  utils::write.csv(timing, file.path(TAB, "tab_task_timing.csv"), row.names = FALSE)
  by_study <- do.call(rbind, lapply(split(timing, timing$study_id), function(x) data.frame(
    study_id = x$study_id[1], tasks = nrow(x), median_sec = stats::median(x$seconds),
    max_sec = max(x$seconds), max_rss_gb = max(x$peak_rss_gb, na.rm = TRUE))))
  utils::write.csv(by_study, file.path(TAB, "tab_task_timing_by_study.csv"), row.names = FALSE)
  log_msg("measured: longest task %.1f min | peak RSS %.2f GB | total %.1f core-hours",
          max(timing$seconds) / 60, max(timing$peak_rss_gb, na.rm = TRUE), sum(timing$seconds) / 3600)
}
missing <- tt$task_id[!tt$done]
if (length(missing)) {
  log_msg("!! %d of %d tasks are not complete (%d failed). Resubmit with:",
          length(missing), nrow(tt), sum(tt$failed))
  cat("     ./submit.sh fit --only ", paste(head(missing, 400), collapse = ","), "\n", sep = "")
  if (length(missing) > 400) cat("     (first 400 shown; ./submit.sh fit resubmits all outstanding)\n")
}

## ---- 2. assemble every study into its raw results file ------------------------------------
for (st in reg) {
  ix <- which(tt$study_id == st$id & tt$done)
  if (!length(ix)) { log_msg("study %-45s no completed tasks", st$id); next }
  parts <- lapply(ix, function(i) readRDS(paste0(stems[i], ".rds")))
  res  <- rbind_fill(lapply(parts, `[[`, "res"))
  pred <- rbind_fill(lapply(parts, `[[`, "pred"))
  if (is.null(res)) next
  ## same row order as a local run: by repetition, fold, then the study's method order
  o <- order(res$rep, res$fold, match(res$method, st$methods))
  res <- res[o, , drop = FALSE]; rownames(res) <- NULL
  if (!is.null(pred)) {
    pred <- pred[order(pred$rep, pred$fold, match(pred$method, st$methods)), , drop = FALSE]
    rownames(pred) <- NULL
  }
  attr(res, "pred") <- pred
  saveRDS(res, file.path(CFG$dir_raw, st$raw_file))
  n_exp <- sum(tt$study_id == st$id)
  log_msg("study %-45s %4d/%d tasks  %5d rows%s", st$id, length(ix), n_exp, nrow(res),
          if (length(ix) < n_exp) "  (INCOMPLETE)" else "")
}
for (b in BLOCKS) {
  fs <- file.path(CFG$dir_raw, vapply(Filter(function(s) s$block == b, reg), `[[`, "", "raw_file"))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) next
  all <- rbind_fill(lapply(fs, readRDS))
  saveRDS(all, file.path(CFG$dir_raw, paste0(b, "_all.rds")))
  utils::write.csv(all, file.path(CFG$dir_raw, paste0(b, "_all.csv")), row.names = FALSE)
}

## ---- 3. convergence study, tables and figures (the local scripts, unchanged) --------------
if (identical(Sys.getenv("RUN_CONVERGENCE", "1"), "1")) {
  log_msg("05_convergence.R ...")
  tryCatch(local(source("05_convergence.R", local = TRUE)),
           error = function(e) log_msg("convergence study failed: %s", conditionMessage(e)))
}
log_msg("08_data_description.R ...")
tryCatch(local(source("08_data_description.R", local = TRUE)),
         error = function(e) log_msg("data description failed: %s", conditionMessage(e)))
log_msg("04_tables_figures.R ...")
local(source("04_tables_figures.R", local = TRUE))

log_msg("COMBINE done in %.1f min | results in %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
        Sys.getenv("MEXGB_OUT"))
if (length(missing)) log_msg("!! results are INCOMPLETE: %d tasks outstanding (see above)", length(missing))
