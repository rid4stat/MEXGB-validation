## Before/after benchmark of the v0.3 fixes (one scenario per call). Run from diagnostics/:
##   Rscript benchmark_before_after.R <gaussian|binomial|cox> <linear|nonlinear> <p> <n_datasets> <splits_per_dataset>
##   optional env: BM=MEXGB_new,XGBoost_new,... (subset)  ROFF=10 (data-set seed offset)  TAG=_b (file suffix)
## "old" = MEXGB/MEGB 0.2 algorithm (crossfit = 0, init = "zero") + old colsample; "new" = v0.3 defaults.
setwd("../analysis")   # run from the diagnostics folder
suppressMessages({source("config.R"); source("R/methods.R"); source("R/cv.R"); library(MEXGB)})
a <- commandArgs(TRUE); fam <- a[1]; dgp <- a[2]; p <- as.integer(a[3]); nrep <- as.integer(a[4]); K <- as.integer(a[5])
new <- CFG
old <- CFG; old$em$crossfit <- 0; old$em$init <- "zero"; old$xgb$tree_method <- "hist"
old$colsample <- function(p) if (p >= 1000) 0.15 else if (p >= 100) 0.5 else 1
meth0 <- if (fam == "cox") list(c("MEXGB","old"),c("MEXGB","new"),c("XGBCox","new"),c("GBM","new")) else if (fam == "gaussian") list(c("MEXGB","old"),c("MEXGB","new"),c("MEGB","old"),c("MEGB","new"),
                                    c("XGBoost","old"),c("XGBoost","new"),c("GBM","new")) else
        list(c("MEXGB","old"),c("MEXGB","new"),c("XGBoost","old"),c("XGBoost","new"),c("GBM","new"))
meth <- if (nzchar(Sys.getenv("BM"))) Filter(function(m) paste0(m[1],"_",m[2]) %in% strsplit(Sys.getenv("BM"),",")[[1]], meth0) else meth0
cal <- if (fam == "binomial") calibrateSim("binomial", target=.85, prevalence=.5, linear=dgp=="linear",
                                           pilot_n=300, pilot_seed=2026) else if (fam == "cox")
  calibrateSim("cox", target=.75, linear=dgp=="linear", random_sd_intercept=sqrt(.5), random_sd_slope=sqrt(.5),
               pilot_n=300, pilot_seed=2026) else NULL
out <- NULL
for (r in 1:nrep) {
  set.seed(1000*(r+as.integer(Sys.getenv("ROFF","0"))) + p)
  dat <- if (fam == "cox") simLong(n=20, p=p, time_points=10, linear=dgp=="linear", family="cox", calibration=cal,
                                   random_sd_intercept=sqrt(.5), random_sd_slope=sqrt(.5), censoring=.3) else
    simLong(n=20, p=p, time_points=10, linear=dgp=="linear", family=fam, calibration=cal)
  X <- as.matrix(dat[, grep("^Var", names(dat))])
  d <- make_data(X, dat$Y, dat$id, dat$time, status = if (fam == "cox") dat$status else NULL)
  folds <- make_folds(d, K, seed=r, unit="within")
  for (k in 1:K) {
    tr <- subset_data(d, -folds[[k]]); te <- subset_data(d, folds[[k]])
    for (m in meth) {
      cfg <- if (m[2]=="old") old else new
      set.seed(k); rr <- run_method(m[1], tr, te, fam, cfg)
      if (!is.na(rr$error)) message(m, ": ", rr$error)
      met <- if (is.null(rr$pred)) NA else if (fam=="gaussian") mse(te$y, rr$pred) else if (fam=="cox") cindex(te$y, te$status, rr$pred) else auc(te$y, rr$pred)
      out <- rbind(out, data.frame(fam, dgp, p, rep=r, fold=k, method=paste0(m[1], if (!m[1] %in% c("GBM","XGBCox")) paste0("_",m[2]) else ""),
                                   metric=met, time=rr$time))
    }
  }
}
saveRDS(out, sprintf("../diagnostics/raw/bench_%s_%s_%d%s.rds", fam, dgp, p, Sys.getenv("TAG")))
agg <- aggregate(cbind(metric,time) ~ method, out, mean); print(agg[order(agg$metric, decreasing = fam!="gaussian"),], digits=4)
