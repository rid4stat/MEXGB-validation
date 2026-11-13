setwd("../../analysis"); suppressMessages({source("config.R"); source("R/methods.R"); source("R/cv.R"); library(MEXGB)})
dgp <- commandArgs(TRUE)[1]; p <- as.integer(commandArgs(TRUE)[2]); res <- NULL
for (r in 1:3) { set.seed(1000*r + p + 7); s <- simLong(n=50, p=p, time_points=20, linear=dgp=="linear")
  X <- as.matrix(s[, grep("^Var", names(s))]); d <- make_data(X, s$Y, s$id, s$time)
  te_rows <- make_folds(d, 1, seed=r, unit="within", n_test=1, test_frac=.3)[[1]]; tr <- subset_data(d, -te_rows); te <- subset_data(d, te_rows)
  row <- c()
  for (mo in c(2, 5)) for (dp in c(2, 3)) { cfg <- CFG; cfg$boost$minobs <- mo; cfg$boost$depth <- dp; cfg$xgb$tree_method <- "exact"
    row[sprintf("GBM_m%d_d%d", mo, dp)] <- mse(te$y, run_method("GBM", tr, te, "gaussian", cfg)$pred)
    row[sprintf("XGB_m%d_d%d", mo, dp)] <- mse(te$y, run_method("XGBoost", tr, te, "gaussian", cfg)$pred) }
  res <- rbind(res, row) }
cat(dgp, p, "\n"); print(round(colMeans(res), 2))
