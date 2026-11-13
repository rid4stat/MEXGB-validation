## v0.4 design benchmark: n=50, t=20, 30% within-subject test; one cell per call
setwd("../../analysis")   # run from diagnostics/v04
suppressMessages({source("config.R"); source("R/methods.R"); source("R/cv.R"); library(MEXGB)})
a <- commandArgs(TRUE); fam <- a[1]; dgp <- a[2]; p <- as.integer(a[3]); nrep <- as.integer(a[4]); K <- as.integer(a[5])
tag <- Sys.getenv("TAG"); S <- CFG$sim; SS <- CFG$surv_sim; cfg <- CFG
if (nzchar(Sys.getenv("OLDBUDGET"))) { cfg$boost$depth <- 3; cfg$boost$minobs <- 5; cfg$xgb$tree_method <- "exact" }
if (nzchar(Sys.getenv("CF"))) cfg$em$crossfit <- as.integer(Sys.getenv("CF")); if (nzchar(Sys.getenv("INNER"))) cfg$em$inner <- as.integer(Sys.getenv("INNER")); if (nzchar(Sys.getenv("DISP"))) cfg$em$dispersion <- Sys.getenv("DISP"); if (nzchar(Sys.getenv("CALB"))) cfg$em$calibrate_binomial <- as.logical(Sys.getenv("CALB")); if (nzchar(Sys.getenv("CAL"))) cfg$em$calibrate <- as.logical(Sys.getenv("CAL"))
meth <- switch(fam, gaussian = c("MEXGB","MEGB","XGBoost","GBM"), binomial = c("MEXGB","XGBoost","GBM"), cox = c("MEXGB","XGBCox","GBM"))
if (nzchar(Sys.getenv("BM"))) meth <- strsplit(Sys.getenv("BM"), ",")[[1]]
cal <- switch(fam, gaussian = NULL,
  binomial = calibrateSim("binomial", target = S$target_auc, prevalence = S$prevalence, linear = dgp == "linear",
                          time_points = S$time_points, random_sd_intercept = S$tau0, random_sd_slope = S$tau1, pilot_n = 300, pilot_seed = 2026),
  cox = calibrateSim("cox", target = SS$target_c, linear = dgp == "linear", time_points = SS$time_points,
                     random_sd_intercept = SS$tau0, random_sd_slope = SS$tau1, pilot_n = 300, pilot_seed = 2026))
xgo <- function(tr, te, zb_tr, zb_te) {     # oracle: xgboost with the TRUE random effects as offset
  obj <- switch(fam, gaussian = "reg:squarederror", binomial = "binary:logistic", cox = "survival:cox")
  lab <- if (fam == "cox") ifelse(tr$status == 1, tr$y, -tr$y) else tr$y
  prm <- list(objective = obj, eta = cfg$boost$shrinkage, max_depth = cfg$boost$depth,
              min_child_weight = if (fam == "gaussian") cfg$boost$minobs else cfg$boost$minobs / 4,
              nthread = 1, tree_method = if (nrow(tr$X) <= 10000) "exact" else "hist", base_score = if (fam == "gaussian") mean(tr$y - zb_tr) else 0.5)
  m <- xgboost::xgb.train(prm, xgboost::xgb.DMatrix(tr$X, label = lab, base_margin = if (fam == "gaussian") NULL else zb_tr), 500, verbose = 0)
  f <- as.numeric(predict(m, te$X, outputmargin = TRUE))
  switch(fam, gaussian = f + zb_te, binomial = plogis(f + zb_te), cox = f + zb_te) }
met <- function(te, pr) switch(fam, gaussian = mse(te$y, pr), binomial = auc(te$y, pr), cox = cindex(te$y, te$status, pr))
out <- NULL
for (r in 1:nrep) {
  set.seed(1000 * r + p + 7)
  s <- if (fam == "cox") simLong(n = SS$n, p = p, time_points = SS$time_points, linear = dgp == "linear", family = "cox", calibration = cal,
                                 random_sd_intercept = SS$tau0, random_sd_slope = SS$tau1, censoring = SS$censoring)
       else simLong(n = S$n, p = p, time_points = S$time_points, linear = dgp == "linear", family = fam, calibration = cal,
                    random_sd_intercept = S$tau0, random_sd_slope = S$tau1, noise_sd = S$sigma)
  X <- as.matrix(s[, grep("^Var", names(s))]); d <- make_data(X, s$Y, s$id, s$time, status = if (fam == "cox") s$status else NULL)
  zb <- rowSums(as.matrix(s[, 4:5]) * attr(s, "b")[s$id, ])
  folds <- make_folds(d, K, seed = r, unit = "within", n_test = 1, test_frac = 0.3)
  for (k in 1:K) {
    te_rows <- folds[[k]]; tr <- subset_data(d, -te_rows); te <- subset_data(d, te_rows)
    for (m in meth) { set.seed(k); rr <- run_method(m, tr, te, fam, cfg)
      if (!is.na(rr$error)) message(m, ": ", rr$error)
      inf <- rr$info
      out <- rbind(out, data.frame(fam, dgp, p, rep = r, fold = k, method = m, metric = if (is.null(rr$pred)) NA else met(te, rr$pred), time = rr$time,
        D11 = if (!is.null(inf$D)) inf$D[1] else NA, D22 = if (!is.null(inf$D) && length(inf$D) > 1) inf$D[2] else NA,
        sigma = if (!is.null(inf$sigma)) inf$sigma else NA, n_iter = if (!is.null(inf$n_iter)) inf$n_iter else NA,
        best_iter = if (!is.null(inf$best_iter)) inf$best_iter else NA)) }
    if (!nzchar(Sys.getenv("BM"))) out <- rbind(out, data.frame(fam, dgp, p, rep = r, fold = k, method = "Oracle", metric = met(te, xgo(tr, te, zb[-te_rows], zb[te_rows])),
        time = NA, D11 = NA, D22 = NA, sigma = NA, n_iter = NA, best_iter = NA))
  }
}
saveRDS(out, sprintf("../diagnostics/v04/raw/b4_%s_%s_%d%s.rds", fam, dgp, p, tag))
ag <- aggregate(cbind(metric, time, D11, D22, sigma, n_iter, best_iter) ~ method, out, mean, na.action = na.pass)
print(ag[order(ag$metric, decreasing = fam != "gaussian"), ], digits = 3)
