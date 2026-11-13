setwd("../../analysis")   # run from diagnostics/v04
suppressMessages({source("config.R"); source("R/methods.R"); source("R/cv.R"); source("R/real_data.R"); source("R/real_survival.R"); source("R/registry.R"); library(MEXGB)})
ds <- commandArgs(TRUE)[1]; nfold <- as.integer(commandArgs(TRUE)[2] %||% "1")
fam <- if (ds %in% CFG$real$datasets_cont) "gaussian" else if (ds %in% CFG$real$datasets_bin) "binomial" else "cox"
d <- load_any_real(ds, CFG); if (fam == "binomial") d$y <- d$y_binary %||% d$y
meth <- switch(fam, gaussian = c("MEXGB","MEGB","XGBoost","GBM"), binomial = c("MEXGB","XGBoost","GBM"), cox = c("MEXGB","XGBCox","GBM"))
folds <- make_folds(d, nfold, seed = 2026, unit = "within", n_test = 1, test_frac = 0.3)
out <- NULL
for (k in seq_len(nfold)) { te_rows <- folds[[k]]; tr <- subset_data(d, -te_rows); te <- subset_data(d, te_rows)
  for (m in meth) { set.seed(k); r <- run_method(m, tr, te, fam, CFG)
    met <- if (is.null(r$pred)) NA else switch(fam, gaussian = mse(te$y, r$pred), binomial = auc(te$y, r$pred), cox = cindex(te$y, te$status, r$pred))
    out <- rbind(out, data.frame(dataset = ds, fold = k, method = m, metric = met, time = r$time, error = r$error,
      D = if (!is.null(r$info$D)) paste(round(r$info$D, 3), collapse = "/") else NA, iters = r$info$n_iter %||% NA)) } }
saveRDS(out, sprintf("../diagnostics/v04/real_check/rc_%s.rds", ds)); print(out)
