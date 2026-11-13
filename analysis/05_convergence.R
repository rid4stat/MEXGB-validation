## =============================================================================
## MEXGB paper 1 -- EM convergence of MEXGB (analogue of Figure 1 of the MEGB paper)
## Evolution of the -2 log-likelihood criterion across EM iterations for different
## n.minobsinnode, linear/nonlinear truth and p in {min, max of CFG$sim$p_grid}, for each
## outcome type, with the MEXGB 0.3 estimation settings of the main analysis (cross-fitted
## E-step, LMM start, colsample, tree method).
## Outputs (one figure per outcome):
##   results/figures/fig_mexgb_loglik_continuous.pdf   (family = "gaussian")
##   results/figures/fig_mexgb_loglik_binary.pdf       (family = "binomial", PQL working -2LL)
##   results/figures/fig_mexgb_loglik_survival.pdf     (family = "cox", PQL working -2LL)
##   results/tables/mexgb_loglik.csv                   (all curves)
## Every fit is cached in results/raw/loglik/ (an interrupted run resumes) and a failing fit
## is reported and skipped instead of stopping the script. Fits run in parallel on
## SLURM_CPUS_PER_TASK (Hali combine job) or CFG$N_CORES cores (Linux/macOS).
## =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(MEXGB); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
S <- CFG$sim
grid <- expand.grid(minobs = c(2, 5, 8), p = range(S$p_grid), dgp = S$dgp,
                    family = c("gaussian", "binomial", "cox"), stringsAsFactors = FALSE)
cache_dir <- file.path(CFG$dir_raw, "loglik"); dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

## calibrations once per (family, dgp)
cals <- list()
for (fam in c("binomial", "cox")) for (dgp in S$dgp) {
  sd1 <- if (fam == "cox") CFG$surv_sim$tau1 else S$tau1
  sd0 <- if (fam == "cox") CFG$surv_sim$tau0 else S$tau0
  cals[[paste(fam, dgp)]] <- MEXGB::calibrateSim(
    fam, target = if (fam == "binomial") S$target_auc else CFG$surv_sim$target_c,
    prevalence = S$prevalence, linear = dgp == "linear", random_sd_intercept = sd0,
    random_sd_slope = sd1, pilot_n = S$pilot_n, pilot_seed = CFG$SEED)
}

one_fit <- function(g) {
  gg <- grid[g, ]
  f <- file.path(cache_dir, sprintf("%s_%s_p%d_minobs%d.rds", gg$family, gg$dgp, gg$p, gg$minobs))
  if (file.exists(f)) return(readRDS(f))
  out <- tryCatch({
    set.seed(CFG$SEED + 77)          # same data for the three n.minobsinnode values
    cox <- gg$family == "cox"
    s <- MEXGB::simLong(n = if (cox) CFG$surv_sim$n else S$n, p = gg$p,
                        time_points = S$time_points, rho_W = S$rho_W, rho_Z = S$rho_Z,
                        random_sd_intercept = if (cox) CFG$surv_sim$tau0 else S$tau0,
                        random_sd_slope = if (cox) CFG$surv_sim$tau1 else S$tau1,
                        noise_sd = S$sigma, linear = gg$dgp == "linear", family = gg$family,
                        calibration = cals[[paste(gg$family, gg$dgp)]],
                        censoring = CFG$surv_sim$censoring,
                        weibull_shape = CFG$surv_sim$weibull_shape,
                        weibull_scale = CFG$surv_sim$weibull_scale)
    Y <- if (cox) cbind(s$Y, s$status) else s$Y
    t0 <- proc.time()[["elapsed"]]
    fit <- MEXGB::MEXGB(X = as.matrix(s[, grep("^Var", names(s))]), Y = Y, id = s$id,
                        Z = as.matrix(s[, c("RandomIntercept", "RandomSlope")]), time = s$time,
                        iter = CFG$em$iter, ntree = CFG$boost$ntree,
                        shrinkage = CFG$boost$shrinkage, interaction.depth = CFG$boost$depth,
                        n.minobsinnode = gg$minobs, lambda = CFG$xgb$lambda,
                        alpha = CFG$xgb$alpha, gamma = CFG$xgb$gamma,
                        subsample = CFG$xgb$subsample, colsample_bytree = CFG$colsample(gg$p),
                        family = gg$family, delta = CFG$em$delta,
                        crossfit = CFG$em$crossfit %||% 5, init = CFG$em$init %||% "lmm",
                        tree_method = CFG$xgb$tree_method %||% "auto", verbose = FALSE)
    data.frame(gg, iteration = seq_along(fit$LL), LL = fit$LL, converged = fit$converged,
               best_iter = fit$best_iter, seconds = proc.time()[["elapsed"]] - t0)
  }, error = function(e) {
    message("!! loglik fit failed (", gg$family, ", ", gg$dgp, ", p = ", gg$p, ", minobs = ",
            gg$minobs, "): ", conditionMessage(e)); NULL })
  if (!is.null(out)) saveRDS(out, f)
  out
}
cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (is.na(cores)) cores <- CFG$N_CORES
cores <- max(1L, min(cores, nrow(grid)))
## longest fits first (p = max) for a shorter makespan
ord <- order(-grid$p, grid$family)
res <- if (cores > 1 && .Platform$OS.type == "unix")
  parallel::mclapply(ord, one_fit, mc.cores = cores, mc.preschedule = FALSE) else lapply(ord, one_fit)
out <- do.call(rbind, Filter(function(x) is.data.frame(x), res))
if (is.null(out) || nrow(out) == 0) stop("no log-likelihood curve could be computed")
utils::write.csv(out, file.path(CFG$dir_tables, "mexgb_loglik.csv"), row.names = FALSE)

## one figure per outcome: panels = truth x dimension, colour = n.minobsinnode
outcomes <- list(gaussian = c("continuous", "-2 log-likelihood criterion"),
                 binomial = c("binary", "-2 PQL working log-likelihood"),
                 cox      = c("survival", "-2 PQL working log-likelihood (Poisson representation)"))
for (fam in names(outcomes)) {
  d <- out[out$family == fam, ]
  if (nrow(d) == 0) { message("no curves for ", fam); next }
  d$Truth <- factor(ifelse(d$dgp == "linear", "Linear", "Nonlinear"), levels = c("Linear", "Nonlinear"))
  d$Dimension <- factor(sprintf("p = %d", d$p), levels = sprintf("p = %d", sort(unique(d$p))))
  ## open circle = iteration returned by MEXGB (convergence, or best likelihood when the EM
  ## stopped on 'patience' / 'iter')
  sel <- d[d$iteration == d$best_iter, ]
  g <- ggplot(d, aes(iteration, LL, colour = factor(minobs))) +
    geom_line() + geom_point(size = 1) +
    geom_point(data = sel, shape = 21, size = 3.2, stroke = 0.9, fill = NA, show.legend = FALSE) +
    facet_wrap(vars(Truth, Dimension), ncol = 2, scales = "free_y",
               labeller = labeller(.multi_line = FALSE)) +   # own y-axis per panel
    scale_x_continuous(breaks = function(l) unique(floor(pretty(seq(1, max(1, l[2])))))) +
    labs(x = "EM iteration", y = outcomes[[fam]][2], colour = "n.minobsinnode",
         title = sprintf("MEXGB: evolution of the log-likelihood criterion, %s outcome", outcomes[[fam]][1]),
         caption = paste0("Open circle: iteration returned by MEXGB (converged, or best likelihood when stopped). ",
                          switch(fam,
                                 gaussian = "",
                                 binomial = "The PQL working likelihood is not monotone; with few binary\nmeasurements per subject the random-slope variance drifts upward over the iterations (as in glmer).",
                                 cox = "Dispersion phi estimated (multiplicative PQL working model)."))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0))
  fn <- file.path(CFG$dir_figures, sprintf("fig_mexgb_loglik_%s.pdf", outcomes[[fam]][1]))
  ggsave(fn, g, width = 8, height = 6)
  message("written ", fn)
}
message("done")
