#!/usr/bin/env Rscript
# =============================================================================
# hali/R/02_fit_task.R -- ONE task = one (study, repetition, fold[, method])
# =============================================================================
#   Rscript 02_fit_task.R [task_id]        # or via SLURM_ARRAY_TASK_ID
#
# Writes exactly one result file at a deterministic path
#   fits/<study>/rRR_fFF_<unit>.rds   (+ .done, or .fail with the error)
# so a rerun overwrites rather than duplicates, and combine can see what is missing.
#
# IDEMPOTENT: a completed task is skipped unless FORCE=1, so a cut-short array is
# simply resubmitted. FAILURE POLICY: a failure writes .fail and exits 0 (one
# method failing on one split is a result, not a reason to lose the array); set
# STRICT_TASK=1 to exit non-zero. Failures of individual METHODS inside a task are
# already recorded per row (column 'error') and do not fail the task.
# =============================================================================
source(file.path(Sys.getenv("HALI_ROOT", "../.."), "hali", "R", "hali_env.R"))

task_id <- current_task_id()
tt  <- readRDS(file.path(PREP_DIR, "tasks.rds"))
if (task_id < 1 || task_id > nrow(tt)) stop(sprintf("task id %d outside 1..%d", task_id, nrow(tt)))
tk  <- tt[task_id, ]
reg <- readRDS(file.path(PREP_DIR, "registry.rds"))
CFG <- readRDS(file.path(PREP_DIR, "cfg.rds"))       # the configuration frozen by prep
CFG$N_CORES <- 1L
st  <- reg[[tk$study_id]]

stem  <- task_stem(FITS_DIR, tk$study_id, tk$rep, tk$fold, tk$unit)
dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
donef <- paste0(stem, ".done"); failf <- paste0(stem, ".fail")

logf <- file.path(LOG_DIR, "tasks", sprintf("task_%05d_%s_r%02d_f%02d_%s.log",
                                            task_id, tk$study_id, tk$rep, tk$fold, tk$unit))
zz <- file(logf, open = "wt"); sink(zz, split = TRUE); sink(zz, type = "message")
log_msg("task %d/%d | %s | rep %d fold %d | methods %s | est %.0f s",
        task_id, nrow(tt), tk$study_id, tk$rep, tk$fold, tk$methods, tk$est_sec)

if (file.exists(donef) && !identical(Sys.getenv("FORCE"), "1")) {
  log_msg("already complete; skipping (FORCE=1 to recompute)")
  quit(save = "no", status = 0)
}
unlink(c(donef, failf))

t0 <- proc.time()[["elapsed"]]
res <- tryCatch(
  run_split(st, CFG, tk$rep, tk$fold, methods = strsplit(tk$methods, ",")[[1]],
            real_dir = REAL_DIR),
  error = function(e) e)
secs <- proc.time()[["elapsed"]] - t0
rss  <- peak_rss_gb()

if (inherits(res, "error")) {
  msg <- conditionMessage(res)
  writeLines(c(sprintf("task=%d", task_id), sprintf("seconds=%.1f", secs),
               sprintf("peak_rss_gb=%.3f", rss), paste("error:", msg)), failf)
  log_msg("FAILED after %.1f s: %s", secs, msg)
  if (identical(Sys.getenv("STRICT_TASK"), "1")) quit(save = "no", status = 1)
  quit(save = "no", status = 0)
}

atomic_save(list(res = res, pred = attr(res, "pred"), task = tk, seconds = secs,
                 peak_rss_gb = rss, host = Sys.info()[["nodename"]],
                 finished = Sys.time()), paste0(stem, ".rds"))
writeLines(sprintf("seconds=%.1f peak_rss_gb=%.3f", secs, rss), donef)
if (!is.null(res) && nrow(res))
  for (i in seq_len(nrow(res)))
    log_msg("  %-12s %s  time %.1f s%s", res$method[i],
            paste(sprintf("%s=%.4g", intersect(c("MSE", "AUC", "Cindex", "IBS"), names(res)),
                          unlist(res[i, intersect(c("MSE", "AUC", "Cindex", "IBS"), names(res))])),
                  collapse = " "),
            res$time[i], if (!is.na(res$error[i])) paste("  ERROR:", res$error[i]) else "")
log_msg("DONE in %.1f s | peak RSS %.2f GB", secs, rss)
