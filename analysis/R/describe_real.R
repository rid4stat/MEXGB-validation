## =============================================================================
## Descriptive summary of the real data sets, computed from the SAME loaders and
## train/test design that the analysis uses (so every number matches the fitted data).
## describe_real_data(cfg) -> data.frame, one row per data set; used by
## 08_data_description.R (tab_real_data.csv / .tex, real_data_description.txt).
## =============================================================================

fmt_range <- function(x, d = 0) sprintf(paste0("%.", d, "f--%.", d, "f"), min(x), max(x))
fmt_msd   <- function(x, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x), stats::sd(x))
fmt_mrange <- function(x) sprintf("%g [%g--%g]", stats::median(x), min(x), max(x))

## expected test-set size of one split under the within-subject rule (deterministic)
test_size <- function(id, cfg) {
  ni <- as.numeric(table(id))
  m <- pmin(pmax(round((cfg$test_frac %||% 0) * ni), cfg$n_test_per_subject %||% 1), ni - 1)
  sum(pmax(m, 0))
}

describe_one <- function(name, cfg) {
  d <- load_any_real(name, cfg)
  ni <- as.numeric(table(d$id))
  surv <- !is.null(d$status)
  unit <- if (surv) cfg$real_surv$cv_unit[[name]] else cfg$real_cv_unit
  nt <- if (identical(unit, "within")) test_size(d$id, cfg) else NA
  out <- data.frame(
    dataset = name,
    outcome_type = if (surv) "survival" else if (name %in% (cfg$real$datasets_bin %||% character(0))) "binary" else "continuous",
    N = length(d$y),
    clusters = length(unique(d$id)),
    per_cluster = fmt_mrange(ni),
    singletons = sum(ni == 1),
    p = ncol(d$X),
    random_effects = if (ncol(d$Z) == 2) "intercept + slope (time)" else "intercept",
    test_per_split = if (is.na(nt)) NA_character_ else sprintf("%d (%.0f%%)", nt, 100 * nt / length(d$y)),
    stringsAsFactors = FALSE)
  if (surv) {
    out$outcome <- sprintf("follow-up median %.2f [%.2f--%.2f] %s", stats::median(d$y), min(d$y), max(d$y), if (name == "Insemination") "months" else "years")
    out$events <- sprintf("%d (%.1f%%)", sum(d$status), 100 * mean(d$status))
    out$binary_prevalence <- NA_character_
    if (!is.null(d$cluster_names))
      out$cluster_sizes <- paste(sprintf("%s=%d", d$cluster_names, as.integer(table(d$id))), collapse = ", ")
  } else {
    is_bin <- name %in% (cfg$real$datasets_bin %||% character(0))
    out$outcome <- if (is_bin) sprintf("binary; %.1f%% events", 100 * mean(d$y)) else
      sprintf("mean (SD) %s; range %s", fmt_msd(d$y), fmt_range(d$y, 2))
    out$events <- NA_character_
    out$binary_prevalence <- if (is_bin) sprintf("%.1f%%", 100 * mean(d$y)) else NA_character_
  }
  out$description <- attr(d, "description") %||% ""
  out
}

describe_real_data <- function(cfg, datasets = c(cfg$real$datasets, cfg$real_surv$datasets)) {
  rows <- lapply(datasets, function(ds) tryCatch(describe_one(ds, cfg), error = function(e) {
    message("  (", ds, " not described: ", conditionMessage(e), ")"); NULL }))
  rbind_fill(rows)
}
