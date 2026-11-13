## =============================================================================
## Study registry -- the SINGLE definition of every simulation scenario and real-data
## analysis. Used by the local scripts (01, 02, 03, 06, 07) AND by the Hali campaign
## (hali/R/*), so a cluster run executes exactly the same studies, seeds and settings as a
## local run: a study's results do not depend on where it was computed.
##
## A study = one raw results file (results/raw/<file>.rds) = one scenario or one
## (data set, outcome). It is evaluated by n_rep x K train/test splits; every split is
## reproducible from (study seed, rep, fold) alone, which is what lets the Hali array split
## a study into independent tasks.
## =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

new_study <- function(block, label, family, methods, n_rep, active, seed, tags, data,
                      raw_file, cv_unit = "cluster") {
  list(id = paste0(block, "_", label), block = block, label = label, family = family,
       methods = methods, n_rep = n_rep, cv_unit = cv_unit, active = active, seed = seed,
       tags = tags, data = data, raw_file = raw_file)
}

## ---- continuous simulation ---------------------------------------------------------
registry_csim <- function(cfg) {
  S <- cfg$sim
  main <- expand.grid(dgp = S$dgp, p = S$p_grid, stringsAsFactors = FALSE)
  main$rho_W <- S$rho_W; main$sigma <- S$sigma; main$study <- "main"
  sens <- expand.grid(rho_W = S$sens_rho, sigma = S$sens_sigma)
  sens <- sens[!(sens$rho_W == S$rho_W & sens$sigma == S$sigma), ]
  sens <- data.frame(dgp = "nonlinear", p = S$sens_p, rho_W = sens$rho_W,
                     sigma = sens$sigma, study = "sensitivity", stringsAsFactors = FALSE)
  scen <- rbind(main, sens)
  scen$label <- sprintf("%s_p%d_rho%.1f_sig%.1f", scen$dgp, scen$p, scen$rho_W, scen$sigma)
  lapply(seq_len(nrow(scen)), function(s) {
    sc <- scen[s, ]
    new_study("csim", sc$label, "gaussian", cfg$methods_cont, cfg$N_REP,
              paste0("Var", 1:S$n_active), cfg$SEED + 100000 * s,
              tags = list(dgp = sc$dgp, p = sc$p, rho_W = sc$rho_W, sigma = sc$sigma,
                          study = sc$study),
              data = list(type = "sim", family = "gaussian", dgp = sc$dgp, p = sc$p,
                          rho_W = sc$rho_W, sigma = sc$sigma),
              raw_file = paste0("csim_", sc$label, ".rds"), cv_unit = cfg$sim_cv_unit)
  })
}

## ---- binary simulation ---------------------------------------------------------------
calibrate_binary <- function(cfg, dgp, prevalence) {
  S <- cfg$sim
  if (is.null(S$target_auc)) return(list(signal = NULL, intercept = NULL, achieved = NA))
  MEXGB::calibrateSim("binomial", target = S$target_auc, prevalence = prevalence,
                      linear = (dgp == "linear"), rel_p = S$n_active,
                      time_points = S$time_points, rho_W = S$rho_W, rho_Z = S$rho_Z,
                      random_sd_intercept = S$tau0, random_sd_slope = S$tau1,
                      pilot_n = S$pilot_n, pilot_seed = cfg$SEED)
}

registry_bsim <- function(cfg, calibrate = TRUE) {
  S <- cfg$sim
  scen <- expand.grid(dgp = S$dgp, p = S$p_grid, stringsAsFactors = FALSE)
  scen$rho_W <- S$rho_W; scen$prevalence <- S$prevalence; scen$study <- "main"
  if (length(S$bin_prev_grid)) {
    scen <- rbind(scen, data.frame(dgp = "nonlinear", p = S$sens_p, rho_W = S$rho_W,
                                   prevalence = setdiff(S$bin_prev_grid, S$prevalence),
                                   study = "prevalence"))
  }
  if (isTRUE(S$run_sensitivity_binary)) {
    scen <- rbind(scen, data.frame(dgp = "nonlinear", p = S$sens_p,
                                   rho_W = setdiff(S$sens_rho, S$rho_W),
                                   prevalence = S$prevalence, study = "sensitivity"))
  }
  scen$label <- sprintf("%s_p%d_rho%.1f_prev%.2f", scen$dgp, scen$p, scen$rho_W, scen$prevalence)
  cal_keys <- unique(scen[, c("dgp", "prevalence")])
  cals <- if (calibrate) lapply(seq_len(nrow(cal_keys)), function(i)
    calibrate_binary(cfg, cal_keys$dgp[i], cal_keys$prevalence[i])) else
      rep(list(NULL), nrow(cal_keys))
  names(cals) <- paste(cal_keys$dgp, cal_keys$prevalence)
  lapply(seq_len(nrow(scen)), function(s) {
    sc <- scen[s, ]
    new_study("bsim", sc$label, "binomial", cfg$methods_bin, cfg$N_REP,
              paste0("Var", 1:S$n_active), cfg$SEED + 500000 + 100000 * s,
              tags = list(dgp = sc$dgp, p = sc$p, rho_W = sc$rho_W,
                          prevalence = sc$prevalence, study = sc$study),
              data = list(type = "sim", family = "binomial", dgp = sc$dgp, p = sc$p,
                          rho_W = sc$rho_W, prevalence = sc$prevalence,
                          cal = cals[[paste(sc$dgp, sc$prevalence)]]),
              raw_file = paste0("bsim_", sc$label, ".rds"), cv_unit = cfg$sim_cv_unit)
  })
}

## ---- survival simulation -------------------------------------------------------------
calibrate_survival <- function(cfg, dgp) {
  S <- cfg$surv_sim
  MEXGB::calibrateSim("cox", target = S$target_c, linear = (dgp == "linear"), rel_p = 6,
                      time_points = S$time_points, rho_W = S$rho_W, rho_Z = S$rho_Z,
                      random_sd_intercept = S$tau0, random_sd_slope = S$tau1,
                      pilot_n = S$pilot_n, pilot_seed = cfg$SEED)
}

registry_ssim <- function(cfg, calibrate = TRUE) {
  S <- cfg$surv_sim
  scen <- expand.grid(dgp = S$dgp, p = S$p_grid, stringsAsFactors = FALSE)
  scen$censoring <- S$censoring; scen$study <- "main"
  if (length(S$cens_grid)) {
    scen <- rbind(scen, data.frame(dgp = "nonlinear", p = S$sens_p,
                                   censoring = setdiff(S$cens_grid, S$censoring),
                                   study = "censoring"))
  }
  scen$label <- sprintf("%s_p%d_cens%.2f", scen$dgp, scen$p, scen$censoring)
  cals <- if (calibrate) lapply(S$dgp, function(g) calibrate_survival(cfg, g)) else
    rep(list(NULL), length(S$dgp))
  names(cals) <- S$dgp
  lapply(seq_len(nrow(scen)), function(s) {
    sc <- scen[s, ]
    new_study("ssim", sc$label, "cox", cfg$methods_surv, cfg$N_REP, paste0("Var", 1:6),
              cfg$SEED + 2000000 + 100000 * s,
              tags = list(dgp = sc$dgp, p = sc$p, censoring = sc$censoring, study = sc$study),
              data = list(type = "sim", family = "cox", dgp = sc$dgp, p = sc$p,
                          censoring = sc$censoring, cal = cals[[sc$dgp]]),
              raw_file = paste0("ssim_", sc$label, ".rds"), cv_unit = cfg$sim_cv_unit)
  })
}

## ---- real data: continuous + binary ---------------------------------------------------
real_datasets <- function(cfg, outcome = c("all", "continuous", "binary")) {
  outcome <- match.arg(outcome)
  ds <- switch(outcome, all = cfg$real$datasets, continuous = cfg$real$datasets_cont %||% cfg$real$datasets,
               binary = cfg$real$datasets_bin %||% cfg$real$datasets)
  if (isTRUE(cfg$SMOKE)) ds <- intersect(ds, c("BodyWeight", "Soybean", "PBC"))
  ds
}

registry_real <- function(cfg) {
  ds <- real_datasets(cfg)
  out <- list()
  for (k in seq_along(ds)) for (fam in c("gaussian", "binomial")) {
    outcome <- if (fam == "gaussian") "continuous" else "binary"
    if (!ds[k] %in% real_datasets(cfg, outcome)) next
    out[[length(out) + 1]] <- new_study(
      "real", paste(ds[k], outcome, sep = "_"), fam,
      if (fam == "gaussian") cfg$methods_cont else cfg$methods_bin, cfg$N_REP_REAL, NULL,
      cfg$SEED + 900000 + 10000 * k + (fam == "binomial") * 5000,
      tags = list(dataset = ds[k], outcome = outcome),
      data = list(type = "real", dataset = ds[k], outcome = outcome),
      raw_file = sprintf("real_%s_%s.rds", ds[k], outcome), cv_unit = cfg$real_cv_unit)
  }
  out
}

## ---- real data: survival ------------------------------------------------------------------
registry_sreal <- function(cfg) {
  ds <- cfg$real_surv$datasets
  lapply(seq_along(ds), function(k)
    new_study("sreal", ds[k], "cox", cfg$methods_surv, cfg$N_REP_REAL, NULL,
              cfg$SEED + 3000000 + 10000 * k,
              tags = list(dataset = ds[k], outcome = "survival"),
              data = list(type = "real", dataset = ds[k], outcome = "survival"),
              raw_file = sprintf("sreal_%s.rds", ds[k]),
              cv_unit = cfg$real_surv$cv_unit[[ds[k]]]))
}

ALL_BLOCKS <- c("csim", "bsim", "ssim", "real", "sreal")

#' Build the registry for the requested blocks. calibrate = FALSE skips the (seconds-long)
#' calibration, for task planning only.
build_registry <- function(cfg, blocks = ALL_BLOCKS, calibrate = TRUE) {
  reg <- list()
  if ("csim"  %in% blocks) reg <- c(reg, registry_csim(cfg))
  if ("bsim"  %in% blocks) reg <- c(reg, registry_bsim(cfg, calibrate))
  if ("ssim"  %in% blocks) reg <- c(reg, registry_ssim(cfg, calibrate))
  if ("real"  %in% blocks) reg <- c(reg, registry_real(cfg))
  if ("sreal" %in% blocks) reg <- c(reg, registry_sreal(cfg))
  names(reg) <- vapply(reg, `[[`, "", "id")
  reg
}

## calibration tables written next to the raw results (read by 04_tables_figures.R)
write_calibration_tables <- function(reg, dir) {
  b <- Filter(function(st) st$block == "bsim" && !is.null(st$data$cal), reg)
  if (length(b)) {
    tab <- unique(do.call(rbind, lapply(b, function(st) data.frame(
      dgp = st$data$dgp, prevalence = st$data$prevalence,
      target_auc = st$data$cal$target %||% NA_real_, kappa = st$data$cal$signal %||% NA_real_,
      c0 = st$data$cal$intercept %||% NA_real_, achieved_auc = st$data$cal$achieved %||% NA_real_))))
    utils::write.csv(tab, file.path(dir, "bsim_calibration.csv"), row.names = FALSE)
  }
  s <- Filter(function(st) st$block == "ssim" && !is.null(st$data$cal), reg)
  if (length(s)) {
    tab <- unique(do.call(rbind, lapply(s, function(st) data.frame(
      dgp = st$data$dgp, target_c = st$data$cal$target, kappa = st$data$cal$signal,
      achieved_c = st$data$cal$achieved))))
    utils::write.csv(tab, file.path(dir, "ssim_calibration.csv"), row.names = FALSE)
  }
}

## ---- data -----------------------------------------------------------------------------------
load_any_real <- function(name, cfg) {
  if (name %in% cfg$real$datasets) load_real(name, cfg) else load_real_surv(name, cfg)
}

#' data_fun(rep) for a study. Real data sets are read from real_dir/<name>.rds when the
#' Hali prep stage has saved them there, otherwise loaded directly.
study_data_fun <- function(st, cfg, real_dir = NULL) {
  dd <- st$data
  if (dd$type == "real") {
    f <- if (!is.null(real_dir)) file.path(real_dir, paste0(dd$dataset, ".rds")) else ""
    d <- if (nzchar(f) && file.exists(f)) readRDS(f) else load_any_real(dd$dataset, cfg)
    if (identical(dd$outcome, "binary"))   # supplied binary outcome if present, else median split
      d$y <- if (!is.null(d$y_binary)) d$y_binary else as.numeric(d$y > stats::median(d$y))
    return(const_data(d))
  }
  seed <- st$seed
  if (dd$family == "gaussian") {
    S <- cfg$sim
    function(rep) {
      set.seed(seed + rep)
      s <- MEXGB::simLong(n = S$n, p = dd$p, rel_p = S$n_active, time_points = S$time_points,
                          rho_W = dd$rho_W, rho_Z = S$rho_Z,
                          random_sd_intercept = S$tau0, random_sd_slope = S$tau1,
                          noise_sd = dd$sigma, linear = (dd$dgp == "linear"), family = "gaussian")
      make_data(s[, -1:-5], s$Y, s$id, s$time)
    }
  } else if (dd$family == "binomial") {
    S <- cfg$sim
    function(rep) {
      set.seed(seed + rep)
      s <- MEXGB::simLong(n = S$n, p = dd$p, rel_p = S$n_active, time_points = S$time_points,
                          rho_W = dd$rho_W, rho_Z = S$rho_Z,
                          random_sd_intercept = S$tau0, random_sd_slope = S$tau1,
                          linear = (dd$dgp == "linear"), family = "binomial",
                          prevalence = dd$prevalence,
                          calibration = if (is.null(dd$cal$signal)) NULL else dd$cal)
      make_data(s[, grep("^Var", names(s))], s$Y, s$id, s$time)
    }
  } else {
    S <- cfg$surv_sim
    function(rep) {
      set.seed(seed + rep)
      s <- MEXGB::simLong(n = S$n, p = dd$p, rel_p = 6, time_points = S$time_points,
                          rho_W = S$rho_W, rho_Z = S$rho_Z,
                          random_sd_intercept = S$tau0, random_sd_slope = S$tau1,
                          linear = (dd$dgp == "linear"), family = "cox",
                          censoring = dd$censoring,
                          weibull_shape = S$weibull_shape, weibull_scale = S$weibull_scale,
                          calibration = dd$cal)
      make_data(s[, grep("^Var", names(s))], s$Y, s$id, s$time, status = s$status)
    }
  }
}

#' Run a whole study locally (all reps x folds) -- what 01/02/03/06/07 call.
run_study <- function(st, cfg, real_dir = NULL) {
  run_cv(data_fun = study_data_fun(st, cfg, real_dir), methods = st$methods,
         family = st$family, cfg = cfg, n_rep = st$n_rep, active = st$active,
         seed = st$seed, tags = st$tags, label = st$label, cv_unit = st$cv_unit)
}

#' Run ONE split (rep, fold) of a study for a subset of its methods -- what a Hali array
#' task calls. Identical to the corresponding split of run_study(): same data (seed + rep),
#' same folds (seed + 1000 rep), same per-method seed (seed + 1000 rep + fold).
run_split <- function(st, cfg, rep, fold, methods = st$methods, real_dir = NULL) {
  d <- study_data_fun(st, cfg, real_dir)(rep)
  folds <- make_folds(d, cfg$K, seed = st$seed + 1000 * rep, unit = st$cv_unit,
                      n_test = cfg$n_test_per_subject %||% 3,
                      test_frac = cfg$test_frac %||% 0)
  if (fold > length(folds)) return(NULL)
  run_fold(d, folds[[fold]], methods, st$family, cfg, active = st$active,
           seed = st$seed + 1000 * rep + fold,
           tags = c(st$tags, list(rep = rep, fold = fold)))
}

## ---- planning: approximate cost of one split (seconds), for ordering and sizing ----------
## v0.4: per-split seconds on one core at N_train = 700 (n = 50 x 14 training measurements).
## MEXGB, MEGB, XGBoost/XGBCox and GBM were measured (diagnostics/v04); the other methods are
## the v0.3 measurements at N_train = 180 scaled by (700/180)^1.3 ~ 5.9. Real data scale with
## their training size relative to 700. Only used to order tasks, to decide which studies are
## split by method, and for the plan's resource projection.
## MEXGB 0.5 timings measured in the build container (diagnostics/v05), seconds per split
COST_TABLE <- local({
  sc <- (700 / 180)^1.3
  list(
  gaussian = list(p6 = c(MEXGB = 200, XGBoost = 0.3, MEGB = 7, GBM = 0.2,
                         c(MERF = 1, REEMforest = 7, RF = 0.1, GPBoost = 1, mboost = 0.2, glmmLasso = 0.2, LMM = 0.1) * sc),
                  p170 = c(MEXGB = 400, XGBoost = 3, MEGB = 38, GBM = 1.4,
                           c(MERF = 5, REEMforest = 7, RF = 2.5, GPBoost = 1.5, mboost = 0.3, glmmLasso = 1.5, LMM = 0.1) * sc),
                  p2000 = c(MEXGB = 1800, XGBoost = 34, MEGB = 1015, GBM = 22,
                            c(MERF = 160, REEMforest = 135, RF = 13, GPBoost = 5, mboost = 3, glmmLasso = 5, LMM = 0.1) * sc)),
  binomial = list(p6 = c(MEXGB = 25, XGBoost = 0.2, GBM = 0.1,
                         c(GMEXGBoost = 3, `GLMM-PQL` = 0.1, glmmLasso = 0.1, RF = 0.1, mboost = 0.1, GPBoost = 3) * sc),
                  p170 = c(MEXGB = 80, XGBoost = 2.7, GBM = 1.4,
                           c(GMEXGBoost = 10, `GLMM-PQL` = 0.1, glmmLasso = 1, RF = 0.3, mboost = 0.1, GPBoost = 3) * sc),
                  p2000 = c(MEXGB = 500, XGBoost = 33, GBM = 18.5,
                            c(GMEXGBoost = 60, `GLMM-PQL` = 0.1, glmmLasso = 5, RF = 3, mboost = 2, GPBoost = 8) * sc)),
  cox = list(p6 = c(MEXGB = 40, XGBCox = 0.2, GBM = 0.1,
                    c(RSF = 0.2, mboost = 0.2, CoxLasso = 0.1, CoxFrailty = 0.5) * sc),
             p170 = c(MEXGB = 100, XGBCox = 2.6, GBM = 1.3,
                      c(RSF = 0.6, mboost = 0.2, CoxLasso = 2, CoxFrailty = 0.5) * sc),
             p2000 = c(MEXGB = 1000, XGBCox = 33, GBM = 18.4,
                       c(RSF = 7, mboost = 0.5, CoxLasso = 3, CoxFrailty = 0.5) * sc)))
})

## dimensions of the real data sets, for planning without loading them
REAL_DIMS <- list(BodyWeight = c(N = 176, p = 3), Soybean = c(N = 412, p = 4),
                  PBC = c(N = 751, p = 13), Insemination = c(N = 10476, p = 4),
                  OvarianDebulking = c(N = 2013, p = 2000),
                  Sprockett = c(N = 1815, p = 1879), Retinopathy = c(N = 394, p = 6),
                  Readmission = c(N = 861, p = 7), BreastMulti = c(N = 867, p = 2000),
                  OvarianMulti = c(N = 2000, p = 2000))

split_cost <- function(st, methods = st$methods, cfg = NULL) {
  if (st$data$type == "sim") { p <- st$data$p; scale <- 1 }
  else {
    dm <- REAL_DIMS[[st$data$dataset]] %||% c(N = 500, p = 100)
    p <- if (st$data$dataset %in% c("BreastMulti", "OvarianMulti") && !is.null(cfg))
      cfg$real_surv$n_genes else dm[["p"]]
    scale <- (0.7 * dm[["N"]] / 700)^1.3      # trees / forests grow slightly super-linearly in N
  }
  tab <- COST_TABLE[[st$family]]
  bucket <- if (p <= 50) "p6" else if (p <= 500) "p170" else "p2000"
  base <- tab[[bucket]]
  pf <- if (bucket == "p2000") p / 2000 else 1
  sum(base[intersect(methods, names(base))]) * scale * pf + 2   # + start-up
}

#' Task table: one task = one (study, rep, fold); studies whose estimated split cost
#' exceeds split_threshold seconds get one task per (study, rep, fold, METHOD) instead.
#' Ordered longest-first so the long tasks start early (shorter makespan).
build_task_table <- function(reg, cfg, split_threshold = 1800) {
  rows <- list()
  for (st in reg) {
    per_method <- split_cost(st, st$methods, cfg) > split_threshold
    msets <- if (per_method) as.list(st$methods) else list(st$methods)
    for (rep in seq_len(st$n_rep)) for (fold in seq_len(cfg$K)) for (ms in msets)
      rows[[length(rows) + 1]] <- data.frame(
        study_id = st$id, block = st$block, rep = rep, fold = fold,
        methods = paste(ms, collapse = ","),
        unit = if (per_method) ms[1] else "all",
        est_sec = round(split_cost(st, ms, cfg)), stringsAsFactors = FALSE)
  }
  tt <- do.call(rbind, rows)
  tt <- tt[order(-tt$est_sec, tt$study_id, tt$rep, tt$fold, tt$unit), ]
  tt$task_id <- seq_len(nrow(tt))
  rownames(tt) <- NULL
  tt[, c("task_id", setdiff(names(tt), "task_id"))]
}

#' file stem of a task's output
task_stem <- function(fits_dir, study_id, rep, fold, unit)
  file.path(fits_dir, study_id, sprintf("r%02d_f%02d_%s", rep, fold, unit))

## ---- local driver used by 01/02/03/06/07 ----------------------------------------------------
run_block_local <- function(block, cfg) {
  reg <- build_registry(cfg, block)
  write_calibration_tables(reg, cfg$dir_raw)
  res <- list()
  for (i in seq_along(reg)) {
    st <- reg[[i]]
    message(sprintf("[%s %d/%d] %s", block, i, length(reg), st$label))
    res[[st$id]] <- tryCatch(cached(file.path(cfg$dir_raw, st$raw_file), run_study(st, cfg)),
                             error = function(e) {
                               message("!! ", st$id, " skipped: ", conditionMessage(e)); NULL })
  }
  invisible(res)
}

## one-line description of every real data set of a block (written to results/raw)
write_real_summary <- function(datasets, cfg, file, real_dir = NULL) {
  rows <- lapply(datasets, function(ds) {
    f <- if (!is.null(real_dir)) file.path(real_dir, paste0(ds, ".rds")) else ""
    d <- tryCatch(if (nzchar(f) && file.exists(f)) readRDS(f) else load_any_real(ds, cfg),
                  error = function(e) NULL)
    if (is.null(d)) return(NULL)
    data.frame(dataset = ds, N = length(d$y), clusters = length(unique(d$id)), p = ncol(d$X),
               events = if (is.null(d$status)) NA else sum(d$status),
               description = attr(d, "description") %||% "")
  })
  out <- do.call(rbind, rows)
  if (!is.null(out)) utils::write.csv(out, file, row.names = FALSE)
  invisible(out)
}
