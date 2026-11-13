## =============================================================================
## MEXGB paper 1 -- Tables (LaTeX + CSV) and figures (PDF) from the raw results
##
## Reads results/raw/*.rds written by 01_, 02_ and 03_ and writes
##   results/tables/*.tex   (same file names and \label{}s as in paper1.tex / thesis)
##   results/tables/*.csv   (machine-readable summaries)
##   results/figures/*.pdf  (same file names as in paper1.tex / thesis)
##   results/results_summary.txt  (all tables in plain text -- send this file back)
##
## Run from the 'analysis' folder:   Rscript 04_tables_figures.R
## =============================================================================
source("config.R")
source("R/cv.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages(library(ggplot2))

TAB <- CFG$dir_tables; FIG <- CFG$dir_figures
out_txt <- file.path(dirname(TAB), "results_summary.txt")
cat("MEXGB paper 1 -- results summary\nGenerated:", format(Sys.time()), "\n",
    file = out_txt)
log_tab <- function(title, df) {
  cat("\n\n==== ", title, " ====\n", file = out_txt, append = TRUE)
  utils::capture.output(print(df, row.names = FALSE), file = out_txt, append = TRUE)
}

read_raw <- function(pattern) {
  f <- list.files(CFG$dir_raw, pattern = pattern, full.names = TRUE)
  f <- f[!grepl("_all\\.rds$", f)]
  if (length(f) == 0) return(NULL)
  lst <- lapply(f, readRDS)
  list(res = rbind_fill(lst),
       pred = rbind_fill(lapply(lst, function(x) attr(x, "pred"))))
}

se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else stats::sd(x) / sqrt(length(x)) }
mn <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else mean(x) }

summarise_metric <- function(res, metric, by) {
  if (!metric %in% names(res)) return(NULL)
  g <- interaction(res[, by, drop = FALSE], drop = TRUE, sep = "|")
  out <- do.call(rbind, lapply(split(seq_len(nrow(res)), g), function(ix) {
    r <- res[ix, , drop = FALSE]
    data.frame(r[1, by, drop = FALSE], mean = mn(r[[metric]]), se = se(r[[metric]]),
               n_ok = sum(is.finite(r[[metric]])), n = nrow(r),
               n_fail = sum(!is.na(r$error)))
  }))
  out$metric <- metric
  rownames(out) <- NULL
  out
}

fmt <- function(m, s, digits, show_se) {
  ifelse(is.na(m), "--",
         if (show_se) ifelse(is.na(s), formatC(m, format = "f", digits = digits),
                             sprintf("%s\\,(%s)", formatC(m, format = "f", digits = digits),
                                     formatC(s, format = "f", digits = digits)))
         else formatC(m, format = "f", digits = digits))
}

## generic LaTeX table: rows = methods, columns = 'col' keys
tex_table <- function(S, methods, cols, col_heads, caption, label, best = c("min", "max"),
                      digits = 2, show_se = TRUE, file, small = "footnotesize",
                      resize = TRUE) {
  best <- match.arg(best)
  M <- matrix("--", length(methods), length(cols), dimnames = list(methods, cols))
  V <- matrix(NA_real_, length(methods), length(cols), dimnames = list(methods, cols))
  for (i in seq_len(nrow(S))) {
    if (S$method[i] %in% methods && S$col[i] %in% cols) {
      M[S$method[i], S$col[i]] <- fmt(S$mean[i], S$se[i], digits, show_se)
      V[S$method[i], S$col[i]] <- round(S$mean[i], digits)
    }
  }
  for (j in seq_along(cols)) {
    v <- V[, j]
    if (all(is.na(v))) next
    b <- if (best == "min") min(v, na.rm = TRUE) else max(v, na.rm = TRUE)
    w <- which(!is.na(v) & v == b)
    M[w, j] <- paste0("\\textbf{", M[w, j], "}")
  }
  keep <- rowSums(M != "--") > 0
  M <- M[keep, , drop = FALSE]
  lines <- c(sprintf("\\begin{table}[htbp]\\centering\\%s", small),
             sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
             "\\setlength{\\tabcolsep}{4pt}",
             if (resize) "\\resizebox{\\textwidth}{!}{%" else NULL,
             sprintf("\\begin{tabular}{l%s}\\toprule", strrep("r", length(cols))),
             paste0("Method & ", paste(col_heads, collapse = " & "), " \\\\\\midrule"),
             paste0(gsub("GLMM-PQL", "GLMM--PQL", rownames(M)), " & ",
                    apply(M, 1, paste, collapse = " & "), " \\\\"),
             paste0("\\bottomrule\\end{tabular}", if (resize) "}" else ""),
             "\\end{table}")
  writeLines(lines, file.path(TAB, file))
  invisible(M)
}

method_cols <- function(methods) {
  cols <- c("#B2182B", grDevices::hcl.colors(max(1, length(methods) - 1), "Dark 3"))
  stats::setNames(cols[seq_along(methods)], methods)
}

## Boxplots, one facet_wrap panel per scenario / data set, each with its OWN y-axis scale
## (scales = "free"); within every panel the methods are sorted by performance (best on the
## left), ranked by the mean of the metric over splits -- the quantity reported in the tables
## (white diamond). Colours stay attached to the method (MEXGB in red) across panels.
##   facet_vars    : column(s) of df defining the panels (their order = panel order)
##   higher_better : TRUE for AUC / C-index, FALSE for MSE / Brier / IBS / time
box_fig <- function(df, metric, methods, facet_vars, ylab, file, higher_better = FALSE,
                    width = 11, height = 6, ncol = NULL) {
  df <- df[is.finite(df[[metric]]) & df$method %in% methods, , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  lab <- lapply(facet_vars, function(v) as.character(df[[v]]))
  df$panel <- do.call(paste, c(lab, sep = ", "))
  ## panel order: the order of the facet variables' levels (factors) or of first appearance
  ord <- do.call(order, lapply(facet_vars, function(v)
    if (is.factor(df[[v]])) as.integer(df[[v]]) else match(df[[v]], unique(df[[v]]))))
  df$panel <- factor(df$panel, levels = unique(df$panel[ord]))
  ## rank the methods within each panel by the mean of the metric
  mm <- stats::aggregate(df[[metric]], list(panel = df$panel, method = df$method), mean)
  mm <- mm[order(mm$panel, if (higher_better) -mm$x else mm$x), ]
  key <- paste(mm$method, mm$panel, sep = "___")
  df$x <- factor(paste(df$method, df$panel, sep = "___"), levels = key)
  df$method <- factor(df$method, levels = methods)
  if (is.null(ncol)) ncol <- min(3, nlevels(df$panel))
  g <- ggplot(df, aes(x = x, y = .data[[metric]], fill = method)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, fill = "white") +
    facet_wrap(~ panel, scales = "free", ncol = ncol) +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = method_cols(methods), guide = "none") +
    labs(x = NULL, y = ylab,
         caption = sprintf("Methods ordered within each panel by mean %s (%s is better); diamond = mean.",
                           ylab, if (higher_better) "higher" else "lower")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey92"),
          plot.caption = element_text(size = 8, colour = "grey30"))
  ggsave(file.path(FIG, file), g, width = width, height = height)
  invisible(g)
}

## description of the simulation train/test design, for the table captions
sim_design_txt <- function(K, R) {
  if (identical(CFG$sim_cv_unit, "within"))
    sprintf("%d train/test splits (%d simulated data sets, each randomly split %d times; every test set holds %d\\%% of the measurements of each individual, i.e. %d of %d)",
            K * R, R, K, round(100 * CFG$test_frac), round(CFG$test_frac * CFG$sim$time_points), CFG$sim$time_points)
  else sprintf("%d subject-level train/test splits (%d repetitions of %d-fold CV)", K * R, R, K)
}

## description of the real-data train/test design, for the table captions
real_within <- function(unit = CFG$real_cv_unit) identical(unit, "within")
real_design_txt <- function() {
  if (real_within())
    sprintf("%d random within-subject train/test splits (%d repetitions $\\times$ %d splits; every test set holds %d\\%% of the measurements of each subject, at least one measurement of every subject kept in training; seen subjects predicted with $f(x)+z'\\hat b_i$)",
            CFG$K * CFG$N_REP_REAL, CFG$N_REP_REAL, CFG$K, round(100 * CFG$test_frac))
  else "subject-level train/test splits (unseen subjects)"
}

lab_dgp <- function(x) ifelse(x == "linear", "lin", "non")

## =============================================================================
## 1. Continuous simulation
## =============================================================================
cs <- read_raw("^csim_.*\\.rds$")
if (!is.null(cs)) {
  res <- cs$res
  mc <- CFG$methods_cont
  main <- res[res$study == "main", ]
  main$col <- sprintf("%s_p%d", lab_dgp(main$dgp), main$p)
  pgrid <- sort(unique(main$p))
  cols <- c(sprintf("lin_p%d", pgrid), sprintf("non_p%d", pgrid))
  heads <- c(sprintf("lin $p{=}%d$", pgrid), sprintf("non $p{=}%d$", pgrid))
  K <- CFG$K; R <- CFG$N_REP

  S <- summarise_metric(main, "MSE", c("method", "col"))
  tex_table(S, mc, cols, heads,
            sprintf("Continuous simulation: test MSE ($\\pm$ SE) over %s. Lower is better; best per column in bold.", sim_design_txt(K, R)),
            "tab:csim", "min", 2, TRUE, "csim_mse_paper.tex")
  utils::write.csv(S, file.path(TAB, "csim_mse.csv"), row.names = FALSE); log_tab("Continuous simulation: MSE", S)

  hd <- pgrid[pgrid > CFG$sim$n_active]
  cols_v <- c(sprintf("lin_p%d", hd), sprintf("non_p%d", hd))
  heads_v <- c(sprintf("lin $p{=}%d$", hd), sprintf("non $p{=}%d$", hd))
  S <- summarise_metric(main, "TPR", c("method", "col"))
  tex_table(S, mc, cols_v, heads_v,
            "Continuous simulation: variable-selection true positive rate (\\%, $\\pm$ SE), top-6 ranked predictors. Higher is better.",
            "tab:csim_tpr", "max", 1, TRUE, "csim_tpr.tex")
  utils::write.csv(S, file.path(TAB, "csim_tpr.csv"), row.names = FALSE); log_tab("Continuous simulation: TPR", S)
  S <- summarise_metric(main, "FPR", c("method", "col"))
  tex_table(S, mc, cols_v, heads_v,
            "Continuous simulation: variable-selection false positive rate (\\%, $\\pm$ SE), top-6 ranked predictors. Lower is better.",
            "tab:csim_fpr", "min", 2, TRUE, "csim_fpr.tex")
  utils::write.csv(S, file.path(TAB, "csim_fpr.csv"), row.names = FALSE); log_tab("Continuous simulation: FPR", S)

  S <- summarise_metric(main, "time", c("method", "col"))
  tex_table(S, mc, cols, heads, "Continuous simulation: mean training time (seconds) per split.",
            "tab:csim_ct", "min", 2, FALSE, "csim_ct.tex")
  utils::write.csv(S, file.path(TAB, "csim_ct.csv"), row.names = FALSE); log_tab("Continuous simulation: time (s)", S)

  ## sensitivity (include the main-study configuration at the same p as reference)
  sens <- res[(res$study == "sensitivity") |
                (res$study == "main" & res$dgp == "nonlinear" & res$p == CFG$sim$sens_p), ]
  if (any(sens$study == "sensitivity")) {
    sens$col <- sprintf("rho%.1f_sig%.1f", sens$rho_W, sens$sigma)
    cfgs <- unique(sens[order(sens$sigma, sens$rho_W), c("rho_W", "sigma", "col")])
    S <- summarise_metric(sens, "MSE", c("method", "col"))
    tex_table(S, mc, cfgs$col,
              sprintf("$\\rho_W{=}%.1f,\\sigma{=}%.1f$", cfgs$rho_W, cfgs$sigma),
              sprintf("Continuous simulation sensitivity sub-study (nonlinear, $p{=}%d$): test MSE ($\\pm$ SE) under varying within-subject correlation $\\rho_W$ and residual scale $\\sigma$. Best per column in bold.", CFG$sim$sens_p),
              "tab:csim_sens", "min", 2, TRUE, "csim_sens.tex")
    utils::write.csv(S, file.path(TAB, "csim_sens.csv"), row.names = FALSE); log_tab("Continuous sensitivity: MSE", S)
    sens <- sens[order(sens$sigma, sens$rho_W), ]
    sens$facet_lab <- factor(sprintf("rho_W = %.1f, sigma = %.1f", sens$rho_W, sens$sigma),
                             levels = unique(sprintf("rho_W = %.1f, sigma = %.1f", sens$rho_W, sens$sigma)))
    box_fig(sens, "MSE", mc, "facet_lab", "Test MSE", "fig_csim_sens.pdf", FALSE, 12, 7)
  }

  main$Scenario <- ifelse(main$dgp == "linear", "Linear", "Nonlinear")
  main$Dimension <- factor(sprintf("p = %d", main$p), levels = sprintf("p = %d", pgrid))
  box_fig(main, "MSE", mc, c("Scenario", "Dimension"), "Test MSE", "fig_csim_mse.pdf", FALSE, 12, 7, ncol = length(pgrid))
  box_fig(main, "time", mc, c("Scenario", "Dimension"), "Training time (s)", "fig_csim_time.pdf", FALSE, 12, 7, ncol = length(pgrid))
  fails <- stats::aggregate(!is.na(error) ~ method + col, data = main, FUN = sum)
  log_tab("Continuous simulation: number of failed fits", fails)
}

## =============================================================================
## 2. Binary simulation
## =============================================================================
bs <- read_raw("^bsim_.*\\.rds$")
if (!is.null(bs)) {
  res <- bs$res
  mb <- CFG$methods_bin
  main <- res[res$study == "main", ]
  main$col <- sprintf("%s_p%d", lab_dgp(main$dgp), main$p)
  pgrid <- sort(unique(main$p))
  cols <- c(sprintf("lin_p%d", pgrid), sprintf("non_p%d", pgrid))
  heads <- c(sprintf("lin $p{=}%d$", pgrid), sprintf("non $p{=}%d$", pgrid))
  K <- CFG$K; R <- CFG$N_REP

  S <- summarise_metric(main, "AUC", c("method", "col"))
  tex_table(S, mb, cols, heads,
            sprintf("Binary simulation: test AUC ($\\pm$ SE) over %s. Higher is better; best per column in bold.", sim_design_txt(K, R)),
            "tab:bsim", "max", 3, TRUE, "bsim_auc_paper.tex")
  utils::write.csv(S, file.path(TAB, "bsim_auc.csv"), row.names = FALSE); log_tab("Binary simulation: AUC", S)
  S <- summarise_metric(main, "Brier", c("method", "col"))
  tex_table(S, mb, cols, heads, "Binary simulation: Brier score ($\\pm$ SE; lower is better).",
            "tab:bsim_brier", "min", 3, TRUE, "bsim_brier.tex")
  utils::write.csv(S, file.path(TAB, "bsim_brier.csv"), row.names = FALSE); log_tab("Binary simulation: Brier", S)
  S <- summarise_metric(main, "Sens", c("method", "col"))
  tex_table(S, mb, cols, heads, "Binary simulation: classification true positive rate (\\%, sensitivity) at threshold $0.5$.",
            "tab:bsim_tpr", "max", 1, TRUE, "bsim_tpr.tex")
  utils::write.csv(S, file.path(TAB, "bsim_sens.csv"), row.names = FALSE); log_tab("Binary simulation: classification TPR", S)
  S <- summarise_metric(main, "FPRc", c("method", "col"))
  utils::write.csv(S, file.path(TAB, "bsim_class_fpr.csv"), row.names = FALSE); log_tab("Binary simulation: classification FPR", S)

  hd <- pgrid[pgrid > CFG$sim$n_active]
  cols_v <- c(sprintf("lin_p%d", hd), sprintf("non_p%d", hd))
  heads_v <- c(sprintf("lin $p{=}%d$", hd), sprintf("non $p{=}%d$", hd))
  S <- summarise_metric(main, "TPR", c("method", "col"))
  tex_table(S, mb, cols_v, heads_v,
            "Binary simulation: variable-selection true positive rate (\\%, $\\pm$ SE), top-6 ranked predictors.",
            "tab:bsim_varsel", "max", 1, TRUE, "bsim_varsel.tex")
  utils::write.csv(S, file.path(TAB, "bsim_varsel_tpr.csv"), row.names = FALSE); log_tab("Binary simulation: variable-selection TPR", S)
  S <- summarise_metric(main, "time", c("method", "col"))
  tex_table(S, mb, cols, heads, "Binary simulation: mean training time (seconds) per split.",
            "tab:bsim_ct", "min", 2, FALSE, "bsim_ct.tex")
  utils::write.csv(S, file.path(TAB, "bsim_ct.csv"), row.names = FALSE); log_tab("Binary simulation: time (s)", S)

  main$Scenario <- ifelse(main$dgp == "linear", "Linear", "Nonlinear")
  main$Dimension <- factor(sprintf("p = %d", main$p), levels = sprintf("p = %d", pgrid))
  box_fig(main, "AUC", mb, c("Scenario", "Dimension"), "Test AUC", "fig_bsim_auc.pdf", TRUE, 12, 7, ncol = length(pgrid))
  box_fig(main, "Brier", mb, c("Scenario", "Dimension"), "Brier score", "fig_bsim_brier.pdf", FALSE, 12, 7, ncol = length(pgrid))
  fails <- stats::aggregate(!is.na(error) ~ method + col, data = main, FUN = sum)
  log_tab("Binary simulation: number of failed fits", fails)

  calf <- file.path(CFG$dir_raw, "bsim_calibration.csv")
  if (file.exists(calf)) log_tab("Binary simulation: calibration (signal kappa, intercept c0)", utils::read.csv(calf))

  ## prevalence study (nonlinear, p = sens_p): main-study prevalence + extra prevalences
  if (any(res$study == "prevalence")) {
    pv <- res[(res$study == "prevalence") |
                (res$study == "main" & res$dgp == "nonlinear" & res$p == CFG$sim$sens_p), ]
    pv$col <- sprintf("prev%.2f", pv$prevalence)
    prevs <- sort(unique(pv$prevalence))
    for (met in c("AUC", "Brier")) {
      S <- summarise_metric(pv, met, c("method", "col"))
      tex_table(S, mb, sprintf("prev%.2f", prevs), sprintf("$\\pi{=}%.2f$", prevs),
                sprintf("Binary simulation (nonlinear, $p{=}%d$): test %s ($\\pm$ SE) by imposed prevalence $\\pi$ (population AUC calibrated to %.2f).",
                        CFG$sim$sens_p, met, CFG$sim$target_auc %||% NA),
                paste0("tab:bsim_prev_", tolower(met)), if (met == "AUC") "max" else "min", 3, TRUE,
                paste0("bsim_prev_", tolower(met), ".tex"), resize = FALSE)
      utils::write.csv(S, file.path(TAB, paste0("bsim_prev_", tolower(met), ".csv")), row.names = FALSE)
      log_tab(paste("Binary prevalence study:", met), S)
    }
    pv$facet_lab <- sprintf("prevalence = %.2f", pv$prevalence)
    box_fig(pv, "AUC", mb, "facet_lab", "Test AUC", "fig_bsim_prev_auc.pdf", TRUE, 12, 4.5)
  }
}

## =============================================================================
## 2b. Survival simulation
## =============================================================================
ss <- read_raw("^ssim_.*\\.rds$")
if (!is.null(ss)) {
  res <- ss$res
  msv <- CFG$methods_surv
  main <- res[res$study == "main", ]
  main$col <- sprintf("%s_p%d", lab_dgp(main$dgp), main$p)
  pgrid <- sort(unique(main$p))
  cols <- c(sprintf("lin_p%d", pgrid), sprintf("non_p%d", pgrid))
  heads <- c(sprintf("lin $p{=}%d$", pgrid), sprintf("non $p{=}%d$", pgrid))
  K <- CFG$K; R <- CFG$N_REP
  S <- summarise_metric(main, "Cindex", c("method", "col"))
  tex_table(S, msv, cols, heads,
            sprintf("Survival simulation: test concordance index ($\\pm$ SE) over %s. Higher is better; best per column in bold.", sim_design_txt(K, R)),
            "tab:ssim", "max", 3, TRUE, "ssim_cindex.tex")
  utils::write.csv(S, file.path(TAB, "ssim_cindex.csv"), row.names = FALSE); log_tab("Survival simulation: C-index", S)
  S <- summarise_metric(main, "IBS", c("method", "col"))
  tex_table(S, msv, cols, heads, "Survival simulation: integrated Brier score ($\\pm$ SE; lower is better).",
            "tab:ssim_ibs", "min", 3, TRUE, "ssim_ibs.tex")
  utils::write.csv(S, file.path(TAB, "ssim_ibs.csv"), row.names = FALSE); log_tab("Survival simulation: IBS", S)
  hd <- pgrid[pgrid > 6]
  cols_v <- c(sprintf("lin_p%d", hd), sprintf("non_p%d", hd))
  heads_v <- c(sprintf("lin $p{=}%d$", hd), sprintf("non $p{=}%d$", hd))
  for (met in c("TPR", "FPR")) {
    S <- summarise_metric(main, met, c("method", "col"))
    tex_table(S, msv, cols_v, heads_v,
              sprintf("Survival simulation: variable-selection %s (\\%%, $\\pm$ SE), top-6 ranked predictors.",
                      if (met == "TPR") "true positive rate" else "false positive rate"),
              paste0("tab:ssim_", tolower(met)), if (met == "TPR") "max" else "min",
              if (met == "TPR") 1 else 2, TRUE, paste0("ssim_", tolower(met), ".tex"))
    utils::write.csv(S, file.path(TAB, paste0("ssim_", tolower(met), ".csv")), row.names = FALSE)
    log_tab(paste("Survival simulation:", met), S)
  }
  S <- summarise_metric(main, "time", c("method", "col"))
  tex_table(S, msv, cols, heads, "Survival simulation: mean training time (seconds) per split.",
            "tab:ssim_ct", "min", 2, FALSE, "ssim_ct.tex")
  utils::write.csv(S, file.path(TAB, "ssim_ct.csv"), row.names = FALSE); log_tab("Survival simulation: time (s)", S)
  calf <- file.path(CFG$dir_raw, "ssim_calibration.csv")
  if (file.exists(calf)) log_tab("Survival simulation: calibration (signal kappa)", utils::read.csv(calf))

  if (any(res$study == "censoring")) {
    cz <- res[(res$study == "censoring") |
                (res$study == "main" & res$dgp == "nonlinear" & res$p == CFG$surv_sim$sens_p), ]
    cz$col <- sprintf("cens%.2f", cz$censoring)
    cl <- sort(unique(cz$censoring))
    S <- summarise_metric(cz, "Cindex", c("method", "col"))
    tex_table(S, msv, sprintf("cens%.2f", cl), sprintf("%d\\%% censored", round(100 * cl)),
              sprintf("Survival simulation (nonlinear, $p{=}%d$): test concordance index by censoring fraction.", CFG$surv_sim$sens_p),
              "tab:ssim_cens", "max", 3, TRUE, "ssim_cens.tex", resize = FALSE)
    utils::write.csv(S, file.path(TAB, "ssim_cens.csv"), row.names = FALSE); log_tab("Survival censoring study: C-index", S)
  }
  main$Scenario <- ifelse(main$dgp == "linear", "Linear", "Nonlinear")
  main$Dimension <- factor(sprintf("p = %d", main$p), levels = sprintf("p = %d", pgrid))
  box_fig(main, "Cindex", msv, c("Scenario", "Dimension"), "Test C-index", "fig_ssim_cindex.pdf", TRUE, 12, 7, ncol = length(pgrid))
  box_fig(main, "IBS", msv, c("Scenario", "Dimension"), "Integrated Brier score", "fig_ssim_ibs.pdf", FALSE, 12, 7, ncol = length(pgrid))
  fails <- stats::aggregate(!is.na(error) ~ method + col, data = main, FUN = sum)
  log_tab("Survival simulation: number of failed fits", fails)
}

## =============================================================================
## 3. Real data
## =============================================================================
rr <- read_raw("^real_.*\\.rds$")
if (!is.null(rr)) {
  res <- rr$res; pred <- rr$pred
  dsets_c <- intersect(CFG$real$datasets_cont %||% CFG$real$datasets, unique(res$dataset[res$outcome == "continuous"]))
  dsets_b <- intersect(CFG$real$datasets_bin %||% CFG$real$datasets, unique(res$dataset[res$outcome == "binary"]))

  ## continuous: fold-level MSE
  rc <- res[res$outcome == "continuous", ]
  Sc <- summarise_metric(rc, "MSE", c("method", "dataset")); Sc$col <- Sc$dataset
  tex_table(Sc, CFG$methods_cont, dsets_c, dsets_c,
            paste0("Real longitudinal data, continuous outcome: test MSE ($\\pm$ SE) over ", real_design_txt(), ". Lower is better; best per data set in bold. ``--'' denotes a method not applicable (or without a stable fit) in that regime."),
            "tab:creal", "min", 2, TRUE, "creal.tex", small = "small", resize = FALSE)
  utils::write.csv(Sc, file.path(TAB, "creal_mse.csv"), row.names = FALSE); log_tab("Real data: continuous MSE", Sc)
  St <- summarise_metric(rc, "time", c("method", "dataset"))
  utils::write.csv(St, file.path(TAB, "creal_time.csv"), row.names = FALSE); log_tab("Real data: continuous time (s)", St)

  ## binary. Within-subject design: every test set holds measurements of every subject, so the
  ## AUC/Brier are computed per split (the random splits overlap and must not be pooled).
  ## Subject-level CV: AUC of the predictions pooled over the K folds of each repetition
  ## (single-subject test folds may contain one class only).
  rb <- res[res$outcome == "binary", ]
  pb <- pred[pred$outcome == "binary", ]
  pooled <- if (real_within()) rb[, c("dataset", "method", "rep", "AUC", "Brier", "error")] else do.call(rbind, lapply(split(pb, list(pb$dataset, pb$method, pb$rep), drop = TRUE),
                                  function(x) data.frame(dataset = x$dataset[1], method = x$method[1],
                                                         rep = x$rep[1], AUC = auc(x$y, x$pred),
                                                         Brier = brier(x$y, x$pred), error = NA)))
  Sb <- summarise_metric(pooled, "AUC", c("method", "dataset")); Sb$col <- Sb$dataset
  tex_table(Sb, CFG$methods_bin, dsets_b, dsets_b,
            if (real_within()) paste0("Real longitudinal data, binary outcome (median split): test AUC ($\\pm$ SE) over ", real_design_txt(), ". Higher is better; best per data set in bold.") else
            "Real longitudinal data, binary outcome (median split): test AUC ($\\pm$ SE over repetitions) of the subject-level cross-validated predictions. Higher is better; best per data set in bold.",
            "tab:breal", "max", 3, TRUE, "breal.tex", small = "small", resize = FALSE)
  utils::write.csv(Sb, file.path(TAB, if (real_within()) "breal_auc.csv" else "breal_auc_pooled.csv"), row.names = FALSE)
  log_tab(if (real_within()) "Real data: binary AUC (per split)" else "Real data: binary AUC (pooled per repetition)", Sb)
  Sbf <- summarise_metric(rb, "AUC", c("method", "dataset"))
  utils::write.csv(Sbf, file.path(TAB, "breal_auc_foldlevel.csv"), row.names = FALSE); log_tab("Real data: binary AUC (fold level; single-class folds excluded)", Sbf)
  Sbb <- summarise_metric(pooled, "Brier", c("method", "dataset"))
  utils::write.csv(Sbb, file.path(TAB, "breal_brier.csv"), row.names = FALSE); log_tab(if (real_within()) "Real data: binary Brier (per split)" else "Real data: binary Brier (pooled)", Sbb)

  ## combined table as in the paper (tab:real)
  allm <- unique(c(CFG$methods_cont, CFG$methods_bin))
  cc <- rbind(transform(Sc, col = paste0("c_", dataset)), transform(Sb, col = paste0("b_", dataset)))
  M <- matrix("--", length(allm), length(dsets_c) + length(dsets_b),
              dimnames = list(allm, c(paste0("c_", dsets_c), paste0("b_", dsets_b))))
  for (i in seq_len(nrow(cc))) if (is.finite(cc$mean[i]))
    M[cc$method[i], cc$col[i]] <- formatC(cc$mean[i], format = "f",
                                          digits = ifelse(grepl("^b_", cc$col[i]), 3, 2))
  for (j in colnames(M)) {
    v <- suppressWarnings(as.numeric(M[, j])); if (all(is.na(v))) next
    b <- if (grepl("^c_", j)) min(v, na.rm = TRUE) else max(v, na.rm = TRUE)
    w <- which(!is.na(v) & v == b); M[w, j] <- paste0("\\textbf{", M[w, j], "}")
  }
  M <- M[rowSums(M != "--") > 0, , drop = FALSE]
  writeLines(c("\\begin{table}[htbp]\\centering\\small",
               paste0("\\caption{Real-data results over ", real_design_txt(), ". Continuous: test MSE (lower better). Binary: test AUC (higher better). Best per column in bold. ``--'' denotes a method not applicable in that dimensional regime or to that outcome.}"),
               "\\label{tab:real}",
               sprintf("\\begin{tabular}{l%s}", strrep("r", ncol(M))), "\\toprule",
               sprintf("& \\multicolumn{%d}{c}{Continuous (test MSE)} & \\multicolumn{%d}{c}{Binary (test AUC)}\\\\",
                       length(dsets_c), length(dsets_b)),
               sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}", length(dsets_c) + 1,
                       length(dsets_c) + 2, length(dsets_c) + length(dsets_b) + 1),
               paste0("Method & ", paste(c(dsets_c, dsets_b), collapse = " & "), "\\\\"), "\\midrule",
               paste0(gsub("GLMM-PQL", "GLMM--PQL", rownames(M)), " & ", apply(M, 1, paste, collapse = " & "), "\\\\"),
               "\\bottomrule", "\\end{tabular}", "\\end{table}"),
             file.path(TAB, "real_combined.tex"))

  rc$dataset <- factor(rc$dataset, levels = dsets_c)
  box_fig(rc, "MSE", CFG$methods_cont, "dataset", "Test MSE", "fig_creal_mse.pdf", FALSE, 12, 4.5)
  pooled$dataset <- factor(pooled$dataset, levels = dsets_b)
  box_fig(pooled, "AUC", CFG$methods_bin, "dataset", if (real_within()) "Test AUC" else "Test AUC (pooled per repetition)",
          "fig_breal_auc.pdf", TRUE, 12, 4.5)

  ## selection stability (high-dimensional cohort): share of splits in which a taxon is in
  ## a method's top-10 relative-influence list
  hd <- res[res$dataset == "Sprockett" & !is.na(res$top_features), ]
  if (nrow(hd) > 0) {
    stab <- do.call(rbind, lapply(split(hd, list(hd$method, hd$outcome), drop = TRUE), function(x) {
      tops <- lapply(strsplit(x$top_features, ";"), utils::head, 10)
      tb <- sort(table(unlist(tops)), decreasing = TRUE)
      data.frame(method = x$method[1], outcome = x$outcome[1], feature = names(tb)[1:min(15, length(tb))],
                 freq = as.numeric(tb)[1:min(15, length(tb))] / nrow(x))
    }))
    utils::write.csv(stab, file.path(TAB, "real_sprockett_selection_stability.csv"), row.names = FALSE)
    log_tab("Sprockett: selection stability (top-10 frequency)", stab[stab$method == "MEXGB", ])
  }
  fails <- stats::aggregate(!is.na(error) ~ method + dataset + outcome, data = res, FUN = sum)
  log_tab("Real data: number of failed fits", fails)
  errs <- unique(res[!is.na(res$error), c("method", "dataset", "outcome", "error")])
  if (nrow(errs)) log_tab("Real data: error messages", errs)
}

## =============================================================================
## 4. Real survival data
## =============================================================================
sr <- read_raw("^sreal_.*\\.rds$")
if (!is.null(sr)) {
  res <- sr$res
  dsets <- intersect(CFG$real_surv$datasets, unique(res$dataset))
  Sc <- summarise_metric(res, "Cindex", c("method", "dataset")); Sc$col <- paste0("c_", Sc$dataset)
  Si <- summarise_metric(res, "IBS", c("method", "dataset")); Si$col <- paste0("i_", Si$dataset)
  utils::write.csv(rbind(Sc, Si), file.path(TAB, "sreal.csv"), row.names = FALSE)
  log_tab("Real survival data: C-index", Sc); log_tab("Real survival data: IBS", Si)
  M <- matrix("--", length(CFG$methods_surv), 2 * length(dsets),
              dimnames = list(CFG$methods_surv, c(paste0("c_", dsets), paste0("i_", dsets))))
  for (S in list(Sc, Si)) for (i in seq_len(nrow(S))) if (is.finite(S$mean[i]))
    M[S$method[i], S$col[i]] <- sprintf("%.3f\\,(%.3f)", S$mean[i], S$se[i])
  for (j in colnames(M)) {
    v <- suppressWarnings(as.numeric(sub("\\\\,.*", "", M[, j]))); if (all(is.na(v))) next
    b <- if (grepl("^c_", j)) max(v, na.rm = TRUE) else min(v, na.rm = TRUE)
    w <- which(!is.na(v) & v == b); M[w, j] <- paste0("\\textbf{", M[w, j], "}")
  }
  M <- M[rowSums(M != "--") > 0, , drop = FALSE]
  writeLines(c("\\begin{table}[htbp]\\centering\\small",
               paste0("\\caption{Real clustered time-to-event data: concordance index (higher better) and integrated Brier score (lower better), mean ($\\pm$ SE) over ", if (all(CFG$real_surv$cv_unit == "within")) paste0(real_design_txt(), " (cluster = herd for Insemination, study for BreastMulti and OvarianMulti)") else "cross-validation splits (Retinopathy and Readmission: cluster-level CV, new patients; BreastMulti and OvarianMulti: observation-level CV stratified by study)", ". ``--'' denotes a method not applicable in that dimensional regime.}"),
               "\\label{tab:sreal}", "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\begin{tabular}{l%s}", strrep("r", ncol(M))), "\\toprule",
               sprintf("& \\multicolumn{%d}{c}{C-index} & \\multicolumn{%d}{c}{Integrated Brier score}\\\\",
                       length(dsets), length(dsets)),
               sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}", length(dsets) + 1,
                       length(dsets) + 2, 2 * length(dsets) + 1),
               paste0("Method & ", paste(rep(dsets, 2), collapse = " & "), "\\\\"), "\\midrule",
               paste0(rownames(M), " & ", apply(M, 1, paste, collapse = " & "), "\\\\"),
               "\\bottomrule", "\\end{tabular}}", "\\end{table}"),
             file.path(TAB, "sreal.tex"))
  res$dataset <- factor(res$dataset, levels = dsets)
  box_fig(res, "Cindex", CFG$methods_surv, "dataset", "Test C-index", "fig_sreal_cindex.pdf", TRUE, 12, 4.5, ncol = 4)
  box_fig(res, "IBS", CFG$methods_surv, "dataset", "Integrated Brier score", "fig_sreal_ibs.pdf", FALSE, 12, 4.5, ncol = 4)
  sumf <- file.path(CFG$dir_raw, "sreal_data_summary.csv")
  if (file.exists(sumf)) log_tab("Real survival data: data sets", utils::read.csv(sumf))
  fails <- stats::aggregate(!is.na(error) ~ method + dataset, data = res, FUN = sum)
  log_tab("Real survival data: number of failed fits", fails)
  errs <- unique(res[!is.na(res$error), c("method", "dataset", "error")])
  if (nrow(errs)) log_tab("Real survival data: error messages", errs)
}

message("Tables written to ", TAB, "; figures to ", FIG, "; summary in ", out_txt)
