## =============================================================================
## Metrics, data containers and the repeated (blocked) cross-validation driver
## Families: "gaussian" (continuous), "binomial" (binary), "cox" (right-censored survival)
## =============================================================================

## ---- metrics -------------------------------------------------------------------
mse <- function(y, p) mean((y - p)^2)

auc <- function(y, p) {                       # Mann-Whitney / rank AUC
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier <- function(y, p) mean((p - y)^2)

class_rates <- function(y, p, thr = 0.5) {
  yhat <- as.integer(p > thr)
  c(sens = if (any(y == 1)) 100 * mean(yhat[y == 1] == 1) else NA_real_,
    fpr  = if (any(y == 0)) 100 * mean(yhat[y == 0] == 1) else NA_real_,
    acc  = 100 * mean(yhat == y))
}

## Harrell's concordance index of a risk score (higher score = higher hazard)
cindex <- function(time, status, risk) {
  if (sum(status) == 0) return(NA_real_)
  if (length(unique(risk)) < 2) return(0.5)          # constant risk score (e.g. empty lasso model)
  as.numeric(survival::concordance(survival::Surv(time, status) ~ risk,
                                   reverse = TRUE)$concordance)
}

## Breslow cumulative baseline hazard (same estimator as MEXGB)
breslow_H <- function(time, status, eta, at) {
  r <- exp(eta - max(eta))
  ut <- sort(unique(time[status == 1]))
  if (length(ut) == 0) return(rep(0, length(at)))
  d <- tabulate(match(time[status == 1], ut), nbins = length(ut))
  o <- order(time); ts <- time[o]
  cs <- rev(cumsum(rev(r[o])))
  R <- cs[findInterval(ut, ts, left.open = TRUE) + 1]
  H <- cumsum(d / R) / exp(max(eta)) * 1          # undo the stabilising shift
  c(0, H)[findInterval(at, ut) + 1]
}

## Integrated Brier score (Graf et al. 1999) with inverse-probability-of-censoring weights.
## Survival curves are S(t|x) = exp(-Lambda0(t) exp(lp)), Lambda0 = Breslow estimator from the
## training risk scores (the same calibration for every method). Grid: 50 points up to the
## 75th percentile of the test follow-up times.
ibs <- function(tr_time, tr_status, lp_train, te_time, te_status, lp_test, n_grid = 50) {
  if (sum(te_status) < 2 || is.null(lp_train)) return(NA_real_)
  tau <- min(stats::quantile(te_time, 0.75), max(tr_time[tr_status == 1]))
  grid <- seq(min(te_time), tau, length.out = n_grid)
  H0 <- breslow_H(tr_time, tr_status, lp_train - mean(lp_train), grid)
  S <- exp(-outer(exp(lp_test - mean(lp_train)), H0))           # n_test x grid
  km <- survival::survfit(survival::Surv(te_time, 1 - te_status) ~ 1)
  Gfun <- stats::stepfun(km$time, c(1, km$surv))
  G_t <- pmax(Gfun(grid), 1e-8)
  G_ti <- pmax(Gfun(te_time - 1e-10), 1e-8)                     # G(t_i-)
  bs <- vapply(seq_along(grid), function(k) {
    t <- grid[k]
    w_event <- (te_time <= t & te_status == 1) / G_ti
    w_alive <- (te_time > t) / G_t[k]
    mean(S[, k]^2 * w_event + (1 - S[, k])^2 * w_alive)
  }, numeric(1))
  if (length(grid) < 2 || diff(range(grid)) <= 0) return(NA_real_)
  sum(diff(grid) * (head(bs, -1) + tail(bs, -1)) / 2) / diff(range(grid))
}

## rbind data frames with different columns (missing columns filled with NA)
rbind_fill <- function(lst) {
  lst <- Filter(function(x) !is.null(x) && nrow(x) > 0, lst)
  if (length(lst) == 0) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) {
    for (cc in setdiff(cols, names(x))) x[[cc]] <- NA
    x[, cols, drop = FALSE]
  }))
}

## ---- data container ------------------------------------------------------------
## X: covariates; y: outcome (survival time for "cox"); status: event indicator (cox);
## id: subject/cluster; time: longitudinal time/occasion used in the random-effect design;
## Z: random-effect design, default cbind(1, time); Z = "intercept" for a random intercept.
make_data <- function(X, y, id, time, status = NULL, Z = NULL) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", seq_len(ncol(X)))
  time <- as.numeric(time)
  if (is.null(Z)) Z <- cbind(Intercept = 1, time = time)
  else if (identical(Z, "intercept")) Z <- cbind(Intercept = rep(1, length(time)))
  list(X = X, y = as.numeric(y), id = id, time = time,
       status = if (is.null(status)) NULL else as.numeric(status), Z = as.matrix(Z))
}

subset_data <- function(d, rows) {
  list(X = d$X[rows, , drop = FALSE], y = d$y[rows], id = d$id[rows],
       time = d$time[rows], status = if (is.null(d$status)) NULL else d$status[rows],
       Z = d$Z[rows, , drop = FALSE])
}

## folds / splits as lists of TEST ROWS
##  unit = "within":      K random train/test splits of the same data set; each test set is
##                        obtained by randomly picking, from EVERY subject i with n_i
##                        measurements, max(round(test_frac * n_i), n_test) of them, capped at
##                        n_i - 1 so that every subject keeps >= 1 measurement in training
##                        (v0.4: test_frac = 0.3 -> 6 of 20 in the simulation; the same rule
##                        for every real data set). Subjects with a single
##                        measurement stay in training. Test subjects are seen in training, so
##                        mixed models predict with their random effects f(x) + z'b_i.
##  unit = "cluster":     all observations of a subject/cluster in the same fold (unseen subjects)
##  unit = "observation": observations split at random within every cluster into K folds (seen
##                        clusters; used when there are only a few large clusters, e.g. studies)
make_folds <- function(d, K, seed, unit = c("within", "cluster", "observation"), n_test = 3,
                       test_frac = 0) {
  unit <- match.arg(unit)
  set.seed(seed)
  if (unit == "within") {
    rows_by_id <- split(seq_along(d$id), factor(d$id, levels = unique(d$id)))
    return(lapply(seq_len(K), function(k) sort(unlist(lapply(rows_by_id, function(r) {
      m <- min(max(round(test_frac * length(r)), n_test),
               length(r) - 1)                  # always keep >= 1 measurement in training
      if (m < 1) integer(0) else r[sample.int(length(r), m)]
    }), use.names = FALSE))))
  }
  if (unit == "cluster") {
    ids <- sample(unique(d$id))
    K <- min(K, length(ids))
    grp <- split(ids, rep_len(seq_len(K), length(ids)))
    lapply(grp, function(g) which(d$id %in% g))
  } else {
    fold <- integer(length(d$id))
    for (g in unique(d$id)) {
      r <- which(d$id == g)
      fold[r] <- sample(rep_len(seq_len(K), length(r)))
    }
    split(seq_along(fold), fold)
  }
}

## ---- one fold: fit every applicable method, return a long data.frame ---------------
run_fold <- function(d, te_rows, methods, family, cfg, active = NULL, seed = 1,
                     tags = list()) {
  tr <- subset_data(d, -te_rows)
  te <- subset_data(d, te_rows)
  p <- ncol(d$X)
  out <- list(); preds <- list()
  for (m in methods) {
    if (!method_applicable(m, p, length(tr$y), cfg, family)) next
    set.seed(seed)
    r <- run_method(m, tr, te, family, cfg)
    row <- data.frame(method = m, n_test = length(te$y), time = r$time,
                      error = r$error, stringsAsFactors = FALSE)
    if (family == "gaussian") {
      row$MSE <- if (is.null(r$pred)) NA_real_ else mse(te$y, r$pred)
    } else if (family == "binomial") {
      if (is.null(r$pred)) {
        row[c("AUC", "Brier", "Sens", "FPRc", "Acc")] <- NA_real_
      } else {
        cr <- class_rates(te$y, r$pred)
        row$AUC <- auc(te$y, r$pred); row$Brier <- brier(te$y, r$pred)
        row$Sens <- cr[["sens"]]; row$FPRc <- cr[["fpr"]]; row$Acc <- cr[["acc"]]
      }
    } else {
      if (is.null(r$pred)) {
        row[c("Cindex", "IBS")] <- NA_real_
      } else {
        row$Cindex <- cindex(te$y, te$status, r$pred)
        row$IBS <- tryCatch(ibs(tr$y, tr$status, r$lp_train, te$y, te$status, r$pred),
                            error = function(e) NA_real_)
      }
    }
    if (!is.null(active) && !is.null(r$imp)) {
      imp <- r$imp[r$imp > 0]
      sel <- names(sort(imp, decreasing = TRUE))[seq_len(min(length(active), length(imp)))]
      row$TPR <- 100 * length(intersect(sel, active)) / length(active)
      row$FPR <- if (p > length(active)) 100 * length(setdiff(sel, active)) / (p - length(active)) else NA_real_
    } else {
      row$TPR <- NA_real_; row$FPR <- NA_real_
    }
    if (!is.null(r$imp)) {
      ## keep the 20 most important features for selection-stability summaries
      top <- head(sort(r$imp[r$imp > 0], decreasing = TRUE), 20)
      row$top_features <- paste(names(top), collapse = ";")
    } else row$top_features <- NA_character_
    out[[m]] <- row
    if (!is.null(r$pred))
      preds[[m]] <- data.frame(method = m, id = te$id, time = te$time, y = te$y,
                               status = if (is.null(te$status)) NA else te$status,
                               pred = r$pred, stringsAsFactors = FALSE)
  }
  res <- rbind_fill(out)
  pr <- rbind_fill(preds)
  for (nm in names(tags)) {
    res[[nm]] <- tags[[nm]]
    if (!is.null(pr)) pr[[nm]] <- tags[[nm]]
  }
  rownames(res) <- NULL
  attr(res, "pred") <- pr
  res
}

## ---- repeated K-fold CV on one data set -----------------------------------------------
## data_fun(rep) returns the data container for repetition 'rep' (a fresh simulated data
## set in the simulation, the same data set for real data).
run_cv <- function(data_fun, methods, family, cfg, n_rep, active = NULL, seed = cfg$SEED,
                   tags = list(), label = "", cv_unit = "cluster") {
  tasks <- expand.grid(fold = seq_len(cfg$K), rep = seq_len(n_rep))
  ## parallel workers start empty: they source the analysis code into their global env
  src <- normalizePath(list.files("R", pattern = "\\.R$", full.names = TRUE))
  one <- function(k) {
    if (!exists(".mexgb_code_loaded", envir = globalenv())) {
      for (f in src) sys.source(f, envir = globalenv())
      assign(".mexgb_code_loaded", TRUE, envir = globalenv())
    }
    rep <- tasks$rep[k]; fold <- tasks$fold[k]
    d <- data_fun(rep)
    folds <- make_folds(d, cfg$K, seed = seed + 1000 * rep, unit = cv_unit,
                        n_test = cfg$n_test_per_subject %||% 3,
                        test_frac = cfg$test_frac %||% 0)
    if (fold > length(folds)) return(NULL)
    run_fold(d, folds[[fold]], methods, family, cfg, active = active,
             seed = seed + 1000 * rep + fold,
             tags = c(tags, list(rep = rep, fold = fold)))
  }
  t0 <- Sys.time()
  if (cfg$N_CORES > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    res <- future.apply::future_lapply(seq_len(nrow(tasks)), one, future.seed = TRUE,
                                       future.packages = c("MEXGB"))
  } else {
    res <- lapply(seq_len(nrow(tasks)), one)
  }
  preds <- rbind_fill(lapply(res, function(r) attr(r, "pred")))
  res <- rbind_fill(res)
  rownames(res) <- NULL; if (!is.null(preds)) rownames(preds) <- NULL
  attr(res, "pred") <- preds
  err <- res[!is.na(res$error), , drop = FALSE]
  if (nrow(err) > 0) {
    tab <- table(err$method)
    message(sprintf("  [%s] failed fits: %s", label,
                    paste(sprintf("%s=%d/%d", names(tab), as.integer(tab),
                                  as.integer(table(res$method)[names(tab)])), collapse = ", ")))
    message("    first error per method: ",
            paste(sprintf("%s: %s", names(tab),
                          err$error[match(names(tab), err$method)]), collapse = " | "))
  }
  message(sprintf("  [%s] %d splits done in %.1f min", label, nrow(tasks),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

## data_fun for a fixed (real) data set; force() binds the data in the closure so that the
## parallel workers receive it
const_data <- function(d) { force(d); function(rep) d }

## load-or-run helper so that interrupted runs resume where they stopped
cached <- function(file, expr) {
  if (file.exists(file)) { message("  (cached) ", file); return(readRDS(file)) }
  res <- expr
  saveRDS(res, file)
  res
}

setup_parallel <- function(cfg) {
  if (cfg$N_CORES > 1 && requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = cfg$N_CORES)
    message("Running on ", cfg$N_CORES, " parallel workers (future::multisession).")
  } else {
    message("Running sequentially.")
  }
  options(future.globals.maxSize = 4 * 1024^3)
}
