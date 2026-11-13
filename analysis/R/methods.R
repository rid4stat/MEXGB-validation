## =============================================================================
## Fitting wrappers: every method is fitted with its own published R package.
##
## Each fit_<method>(tr, te, family, cfg) receives
##   tr, te : lists with X (matrix, named columns), y (numeric; 0/1 if binary; observed
##            time if survival), status (survival only), id, time, Z (random-effect design)
## and returns list(pred = predictions for te on the response scale (probability for
## binary; log relative hazard / risk score for survival), imp = named importance vector
## over colnames(X) or NULL, and for survival lp_train = training risk scores, used to
## estimate the Breslow baseline for the integrated Brier score).
## Every mixed-effect method predicts with f(x) + z'b_i, where b_i is the subject's predicted
## random effect when the subject was seen in training (simulation design: three test
## measurements of every subject; multi-study survival data) and 0 for an unseen subject
## (subject-level CV of the real data sets).
## =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

full_imp <- function(v, features) {
  out <- stats::setNames(rep(0, length(features)), features)
  if (is.null(v) || length(v) == 0) return(out)
  v <- v[names(v) %in% features]
  out[names(v)] <- abs(as.numeric(v))
  out[!is.finite(out)] <- 0
  out
}

## standardise train/test with training moments (linear / penalised methods)
standardise <- function(tr, te) {
  mu <- colMeans(tr$X); s <- apply(tr$X, 2, stats::sd)
  keep <- is.finite(s) & s > 0
  Xtr <- scale(tr$X[, keep, drop = FALSE], mu[keep], s[keep])
  Xte <- scale(te$X[, keep, drop = FALSE], mu[keep], s[keep])
  list(Xtr = Xtr, Xte = Xte, keep = colnames(tr$X)[keep])
}

quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

## Add each test observation's predicted random effect to a fixed-effect prediction:
##   eta = f(x) + z' b_i   for subjects seen in training (b_i = row of re_mat for re_ids),
##   eta = f(x)            for unseen subjects (prior mean b_i = 0).
## Z: random-effect design of the test rows; only its first ncol(re_mat) columns are used
## (a random-intercept model uses the intercept column only).
add_re <- function(fixed, Z, id, re_mat, re_ids) {
  re_mat <- as.matrix(re_mat)
  k <- match(as.character(id), as.character(re_ids))
  ok <- !is.na(k)
  out <- as.numeric(fixed)
  if (any(ok))
    out[ok] <- out[ok] + rowSums(as.matrix(Z)[ok, seq_len(ncol(re_mat)), drop = FALSE] *
                                   re_mat[k[ok], , drop = FALSE])
  out
}

## ---------------------------------------------------------------------------
## MEXGB (proposed) -- package MEXGB (xgboost back-end)
## ---------------------------------------------------------------------------
fit_MEXGB <- function(tr, te, family, cfg) {
  p <- ncol(tr$X)
  Y <- if (family == "cox") cbind(tr$y, tr$status) else tr$y
  mx <- cfg$mexgb %||% list()
  fit <- MEXGB::MEXGB(X = tr$X, Y = Y, id = tr$id, Z = tr$Z, time = tr$time,
                      iter = cfg$em$iter, ntree = mx$ntree_max %||% cfg$boost$ntree,
                      linear = mx$linear %||% FALSE, lin_alpha = mx$lin_alpha %||% 1,
                      honest = mx$honest %||% TRUE, honest_d = mx$honest_d %||% TRUE,
                      lin_relax = mx$lin_relax %||% TRUE, honest_offset = mx$honest_offset %||% TRUE,
                      ntree_select = mx$ntree_select %||% "fixed",
                      shrinkage = mx$shrinkage %||% cfg$boost$shrinkage,
                      interaction.depth = cfg$boost$depth,
                      n.minobsinnode = cfg$boost$minobs,
                      lambda = mx$lambda %||% cfg$xgb$lambda, alpha = cfg$xgb$alpha,
                      gamma = cfg$xgb$gamma, subsample = mx$subsample %||% cfg$xgb$subsample,
                      colsample_bytree = mx$colsample %||% cfg$colsample(p),
                      family = family, delta = cfg$em$delta, nthread = 1,
                      crossfit = cfg$em$crossfit %||% 5, init = mx$init %||% "auto",
                      inner = cfg$em$inner %||% 50,
                      dispersion = cfg$em$dispersion %||% "auto",
                      calibrate = cfg$em$calibrate %||% TRUE,
                      calibrate_binomial = cfg$em$calibrate_binomial %||% TRUE,
                      tree_method = cfg$xgb$tree_method %||% "auto",
                      verbose = FALSE)
  if (family == "cox") {
    pred <- predict(fit, X = te$X, Z = te$Z, id = te$id, time = te$time, type = "link")
    lp_train <- predict(fit, X = tr$X, Z = tr$Z, id = tr$id, time = tr$time, type = "link")
  } else {
    pred <- predict(fit, X = te$X, Z = te$Z, id = te$id, time = te$time, type = "response")
    lp_train <- NULL
  }
  list(pred = pred, lp_train = lp_train, imp = full_imp(fit$importance, colnames(tr$X)),
       info = list(n_iter = fit$n_iter, converged = fit$converged, best_iter = fit$best_iter,
                   D = diag(as.matrix(fit$var_random_effects)), sigma = as.numeric(fit$sigma),
                   nrounds = fit$nrounds %||% NA, cal_slope = fit$cal_slope %||% NA,
                   linear = fit$linear, lin_decision = fit$lin_decision, screen_decision = fit$screen_decision,
                   n_screened = length(fit$screened), honest_loss = fit$honest_loss[is.finite(fit$honest_loss)], d_scale = fit$d_scale[is.finite(fit$d_scale)]))
}

## ---------------------------------------------------------------------------
## MEGB -- package MEGB (gbm back-end). MEGB 0.2 returns NULL when the EM does not
## reach the tolerance within 'iter' iterations; we then refit once with a 10x looser
## tolerance, and report NA if it still fails.
## ---------------------------------------------------------------------------
fit_MEGB <- function(tr, te, family, cfg) {
  stopifnot(family == "gaussian")
  ## MEGB >= 0.3: cross-fitted E-step + LMM starting values (same settings as MEXGB)
  extra <- if ("crossfit" %in% names(formals(MEGB::MEGB)))
    c(list(crossfit = cfg$em$crossfit %||% 5, init = cfg$em$init %||% "lmm"),
      if ("inner" %in% names(formals(MEGB::MEGB))) list(inner = cfg$em$inner %||% 50)) else list()
  run <- function(delta) quiet(do.call(MEGB::MEGB, c(list(
    X = tr$X, Y = as.matrix(tr$y), id = tr$id, Z = tr$Z, iter = cfg$em$iter,
    ntree = cfg$boost$ntree, time = tr$time, shrinkage = cfg$boost$shrinkage,
    interaction.depth = cfg$boost$depth, n.minobsinnode = cfg$boost$minobs, cv.folds = 0,
    delta = delta, verbose = FALSE), extra)))
  fit <- run(cfg$em$delta)
  retried <- FALSE
  if (is.null(fit)) { fit <- run(10 * cfg$em$delta); retried <- TRUE }
  if (is.null(fit)) stop("MEGB did not converge (MEGB() returned NULL)")
  ## population prediction for unseen subjects (fixed effect of the gbm fit)
  ## fixed effect from the gbm fit + the subject's BLUP (computed here rather than with
  ## MEGB 0.2's predict(), which matches random effects to subjects by position)
  f <- as.numeric(gbm::predict.gbm(fit$forest, newdata = data.frame(te$X),
                                   n.trees = cfg$boost$ntree))
  pred <- add_re(f, te$Z, te$id, fit$random_effects, fit$id_btilde)
  ri <- gbm::summary.gbm(fit$forest, plotit = FALSE)
  list(pred = pred, imp = full_imp(stats::setNames(ri$rel.inf, as.character(ri$var)),
                                   colnames(tr$X)),
       info = list(retried = retried))
}

## ---------------------------------------------------------------------------
## MERF and REEMforest -- package LongituRF (Capitaine et al. 2021), sto = "none"
## ---------------------------------------------------------------------------
fit_longituRF <- function(tr, te, cfg, which) {
  p <- ncol(tr$X)
  ntree <- if (p > 1000) cfg$real$forest_ntree_hd else cfg$forest$ntree
  f <- if (which == "MERF") LongituRF::MERF else LongituRF::REEMforest
  fit <- quiet(f(X = tr$X, Y = tr$y, id = tr$id, Z = tr$Z, iter = cfg$em$iter,
                 mtry = cfg$forest$mtry_reg(p), ntree = ntree, time = tr$time,
                 sto = "none", delta = cfg$em$delta))
  pred <- add_re(as.numeric(stats::predict(fit$forest, te$X)), te$Z, te$id,
                 fit$random_effects, fit$id_btilde)
  imp <- randomForest::importance(fit$forest, type = 1)[, 1]
  list(pred = pred, imp = full_imp(imp, colnames(tr$X)))
}
fit_MERF       <- function(tr, te, family, cfg) fit_longituRF(tr, te, cfg, "MERF")
fit_REEMforest <- function(tr, te, family, cfg) fit_longituRF(tr, te, cfg, "REEMforest")

## ---------------------------------------------------------------------------
## RF -- package randomForest (independence-assuming)
## ---------------------------------------------------------------------------
fit_RF <- function(tr, te, family, cfg) {
  p <- ncol(tr$X)
  ntree <- if (p > 1000) cfg$real$forest_ntree_hd else cfg$forest$ntree
  if (family == "gaussian") {
    fit <- randomForest::randomForest(x = tr$X, y = tr$y, ntree = ntree,
                                      mtry = cfg$forest$mtry_reg(p), importance = TRUE)
    pred <- as.numeric(stats::predict(fit, te$X))
  } else {
    fit <- randomForest::randomForest(x = tr$X, y = factor(tr$y, levels = c(0, 1)),
                                      ntree = ntree, mtry = cfg$forest$mtry_cls(p),
                                      importance = TRUE)
    pred <- as.numeric(stats::predict(fit, te$X, type = "prob")[, "1"])
  }
  list(pred = pred, imp = full_imp(randomForest::importance(fit, type = 1)[, 1],
                                   colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## GBM -- package gbm (independence-assuming first-order boosting)
## ---------------------------------------------------------------------------
fit_GBM <- function(tr, te, family, cfg) {
  dist <- switch(family, gaussian = "gaussian", binomial = "bernoulli", cox = "coxph")
  if (family == "cox") {
    dtr <- data.frame(stime = tr$y, sstat = tr$status, tr$X)
    form <- stats::as.formula(paste("survival::Surv(stime, sstat) ~",
                                    paste(colnames(dtr)[-(1:2)], collapse = " + ")))
  } else {
    dtr <- data.frame(y = tr$y, tr$X)
    form <- y ~ .
  }
  fit <- gbm::gbm(form, data = dtr, distribution = dist, n.trees = cfg$boost$ntree,
                  interaction.depth = cfg$boost$depth, shrinkage = cfg$boost$shrinkage,
                  n.minobsinnode = cfg$boost$minobs, bag.fraction = cfg$boost$bag.fraction,
                  train.fraction = 1, verbose = FALSE)
  type <- if (family == "cox") "link" else "response"
  pred <- as.numeric(gbm::predict.gbm(fit, newdata = data.frame(te$X),
                                      n.trees = cfg$boost$ntree, type = type))
  lp_train <- if (family == "cox")
    as.numeric(gbm::predict.gbm(fit, newdata = data.frame(tr$X), n.trees = cfg$boost$ntree,
                                type = "link")) else NULL
  ri <- gbm::summary.gbm(fit, plotit = FALSE)
  list(pred = pred, lp_train = lp_train,
       imp = full_imp(stats::setNames(ri$rel.inf, as.character(ri$var)), colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## GPBoost -- package gpboost (tree boosting + grouped random intercept/slope)
## ---------------------------------------------------------------------------
fit_GPBoost <- function(tr, te, family, cfg) {
  lik <- if (family == "gaussian") "gaussian" else "bernoulli_logit"
  if (isTRUE(cfg$gpb$random_slope)) {
    gp <- gpboost::GPModel(group_data = tr$id, group_rand_coef_data = tr$time,
                           ind_effect_group_rand_coef = 1, likelihood = lik)
  } else {
    gp <- gpboost::GPModel(group_data = tr$id, likelihood = lik)
  }
  dtrain <- gpboost::gpb.Dataset(data = tr$X, label = tr$y)
  params <- list(learning_rate = cfg$boost$shrinkage, max_depth = cfg$boost$depth,
                 num_leaves = 2^cfg$boost$depth, min_data_in_leaf = cfg$boost$minobs,
                 feature_fraction = cfg$gpb$feature_fraction, num_threads = 1,
                 verbose = -1)
  bst <- quiet(gpboost::gpb.train(data = dtrain, gp_model = gp, nrounds = cfg$boost$ntree,
                                  params = params, verbose = 0))
  pr <- if (isTRUE(cfg$gpb$random_slope)) {
    stats::predict(bst, data = te$X, group_data_pred = te$id,
                   group_rand_coef_data_pred = te$time, predict_var = FALSE,
                   pred_latent = FALSE)
  } else {
    stats::predict(bst, data = te$X, group_data_pred = te$id, predict_var = FALSE,
                   pred_latent = FALSE)
  }
  imp <- tryCatch({
    ti <- gpboost::gpb.importance(bst)
    stats::setNames(ti$Gain, ti$Feature)
  }, error = function(e) NULL)
  list(pred = as.numeric(pr$response_mean), imp = full_imp(imp, colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## mboost -- package mboost: glmboost (component-wise linear boosting), centred
## ---------------------------------------------------------------------------
fit_mboost <- function(tr, te, family, cfg) {
  st <- standardise(tr, te)
  ctrl <- mboost::boost_control(mstop = cfg$mboost$mstop, nu = cfg$mboost$nu)
  if (family == "cox") {
    fit <- quiet(mboost::glmboost(x = st$Xtr, y = survival::Surv(tr$y, tr$status),
                                  family = mboost::CoxPH(), center = TRUE, control = ctrl))
    pred <- as.numeric(stats::predict(fit, newdata = st$Xte, type = "link"))
    lp_train <- as.numeric(stats::predict(fit, newdata = st$Xtr, type = "link"))
    cf <- quiet(stats::coef(fit))
    return(list(pred = pred, lp_train = lp_train,
                imp = full_imp(cf[names(cf) != "(Intercept)"], colnames(tr$X))))
  }
  if (family == "gaussian") {
    fit <- mboost::glmboost(x = cbind("(Intercept)" = 1, st$Xtr), y = tr$y,
                            family = mboost::Gaussian(), center = TRUE, control = ctrl)
  } else {
    fit <- mboost::glmboost(x = cbind("(Intercept)" = 1, st$Xtr),
                            y = factor(tr$y, levels = c(0, 1)),
                            family = mboost::Binomial(), center = TRUE, control = ctrl)
  }
  pred <- as.numeric(stats::predict(fit, newdata = cbind("(Intercept)" = 1, st$Xte),
                                    type = "response"))
  cf <- quiet(stats::coef(fit))
  cf <- cf[names(cf) != "(Intercept)"]
  list(pred = pred, imp = full_imp(cf, colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## glmmLasso -- package glmmLasso (L1-penalised (G)LMM), random intercept + slope
## ---------------------------------------------------------------------------
fit_glmmLasso <- function(tr, te, family, cfg) {
  st <- standardise(tr, te)
  vars <- make.names(st$keep)
  dtr <- data.frame(y = tr$y, id = factor(tr$id), st$Xtr)
  names(dtr) <- c("y", "id", vars)
  form <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")))
  fam <- if (family == "gaussian") stats::gaussian(link = "identity") else stats::binomial(link = "logit")
  fit <- quiet(glmmLasso::glmmLasso(fix = form, rnd = list(id = ~1), data = dtr,
                                    lambda = cfg$glmmlasso$lambda, family = fam,
                                    control = list(print.iter = FALSE)))
  b <- fit$coefficients
  eta <- b["(Intercept)"] + as.numeric(st$Xte %*% b[vars])
  ## random intercepts: glmmLasso names them "id<level>"
  rf <- fit$ranef
  if (!is.null(rf) && length(rf))
    eta <- add_re(eta, te$Z, te$id, matrix(as.numeric(rf), ncol = 1),
                  sub("^id", "", names(rf)))
  pred <- if (family == "gaussian") eta else stats::plogis(eta)
  imp <- stats::setNames(b[vars], st$keep)
  list(pred = as.numeric(pred), imp = full_imp(imp, colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## LMM -- package nlme: lme, REML, random intercept + slope
## ---------------------------------------------------------------------------
fit_LMM <- function(tr, te, family, cfg) {
  stopifnot(family == "gaussian")
  vars <- make.names(colnames(tr$X))
  dtr <- data.frame(y = tr$y, id = factor(tr$id), tr$X); names(dtr) <- c("y", "id", vars)
  dte <- data.frame(te$X); names(dte) <- vars
  form <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")))
  fit <- eval(bquote(nlme::lme(fixed = .(form), random = ~ 1 | id, data = dtr,
                               method = "REML",
                               control = nlme::lmeControl(returnObject = TRUE))))
  rf <- nlme::ranef(fit)
  pred <- add_re(as.numeric(stats::predict(fit, newdata = dte, level = 0)), te$Z, te$id,
                 as.matrix(rf), rownames(rf))
  tt <- summary(fit)$tTable
  tval <- stats::setNames(tt[vars, "t-value"], colnames(tr$X))
  list(pred = pred, imp = full_imp(tval, colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## GLMM-PQL -- package MASS: glmmPQL, logit link, random intercept + slope
## ---------------------------------------------------------------------------
`fit_GLMM-PQL` <- function(tr, te, family, cfg) {
  stopifnot(family == "binomial")
  vars <- make.names(colnames(tr$X))
  dtr <- data.frame(y = tr$y, id = factor(tr$id), tr$X); names(dtr) <- c("y", "id", vars)
  dte <- data.frame(te$X); names(dte) <- vars
  form <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")))
  fit <- quiet(eval(bquote(MASS::glmmPQL(fixed = .(form), random = ~ 1 | id,
                                          family = stats::binomial, data = dtr,
                                          verbose = FALSE))))
  rf <- nlme::ranef(fit)
  eta <- add_re(as.numeric(stats::predict(fit, newdata = dte, level = 0)), te$Z, te$id,
                as.matrix(rf), rownames(rf))
  tt <- summary(fit)$tTable
  tval <- stats::setNames(tt[vars, "t-value"], colnames(tr$X))
  list(pred = stats::plogis(eta), imp = full_imp(tval, colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## GMEXGBoost -- Asadi et al. (2026), BMC Med Res Methodol 26:24.
## Faithful port of the authors' gmexgb() (github.com/faas34188/GMEXGBoost) using the
## packages they use: xgboost (fixed effect) and lme4::glmer (random effects).
## Algorithm (authors' code): (1) GLM on the covariates gives an initial linear
## predictor eta0 (kept fixed); (2) iterate: XGBoost regression of eta0 - Z b on X;
## glmer(y ~ (1 + time | id), offset = f(X)); b <- ranef; stop when the largest
## relative change in b is < toll or after itmax iterations. New subjects:
## eta = glmer intercept + f(X).
## Changes needed to run it here, documented: xgb.train() replaces the legacy
## xgboost(data =, label =) call (removed in xgboost >= 3.0); when p > glm_init_max_p
## the initial GLM is a ridge-penalised logistic GLM (an unpenalised GLM is not estimable
## when p is large: it separates); the initial
## linear predictor is bounded to +-eta_bound; if glmer fails it is refitted with nAGQ = 0,
## then with a random intercept only; if it fails after the first iteration, the last
## successful fit is kept and the iterations stop.
## ---------------------------------------------------------------------------
gmexgb_fit <- function(y, X, id, time, cfg) {
  group <- factor(id)
  lev <- levels(group)
  q <- if (isTRUE(cfg$gmex$random_slope)) 2 else 1
  Zi <- if (q == 2) cbind(1, time) else matrix(1, length(y), 1)
  bi <- matrix(0, q, length(lev), dimnames = list(NULL, lev))

  ## STEP 2: GLM initialisation of eta (authors: unpenalised GLM on the covariates).
  ## (adaptation) when p > glm_init_max_p an unpenalised GLM is not estimable (separation),
  ## so a ridge-penalised logistic GLM (glmnet, alpha = 0, lambda by 5-fold CV) is used.
  p <- ncol(X)
  if (p <= cfg$gmex$glm_init_max_p) {
    glm0 <- suppressWarnings(stats::glm(y ~ ., data = data.frame(y = y, X), family = stats::binomial))
    eta.est <- glm0$family$linkfun(glm0$fitted.values)
  } else {
    sx <- scale(X); sx[!is.finite(sx)] <- 0
    rg <- glmnet::cv.glmnet(sx, y, family = "binomial", alpha = 0, nfolds = 5)
    eta.est <- as.numeric(stats::predict(rg, newx = sx, s = "lambda.min", type = "link"))
  }
  ## (adaptation) bound the GLM linear predictor: under (quasi-)separation it diverges and
  ## the glmer offset then makes PIRLS fail
  eta.est <- pmin(pmax(eta.est, -cfg$gmex$eta_bound), cfg$gmex$eta_bound)

  params <- list(objective = "reg:squarederror", eta = cfg$boost$shrinkage,
                 max_depth = cfg$boost$depth, subsample = cfg$xgb$subsample,
                 lambda = cfg$xgb$lambda, alpha = cfg$xgb$alpha,
                 colsample_bytree = cfg$colsample(p), nthread = 1)
  rand_form <- if (q == 2) y ~ (1 + time | group) else y ~ (1 | group)
  gdat <- data.frame(y = y, group = group, time = time)

  it <- 1; converged <- FALSE
  while (!converged && it < cfg$gmex$itmax) {
    target <- eta.est - rowSums(Zi * t(bi[, as.character(group), drop = FALSE]))
    xgb <- xgboost::xgb.train(params = c(params, list(base_score = mean(target))),
                              data = xgboost::xgb.DMatrix(X, label = target),
                              nrounds = cfg$boost$ntree, verbose = 0)
    fx <- as.numeric(stats::predict(xgb, X))
    fit_glmer <- function(form, nagq = 1) tryCatch(suppressMessages(suppressWarnings(
      lme4::glmer(form, data = gdat, family = stats::binomial, offset = fx, nAGQ = nagq))),
      error = function(e) NULL)
    gfit <- fit_glmer(rand_form)
    if (is.null(gfit)) gfit <- fit_glmer(rand_form, nagq = 0)   # PIRLS-only (more stable)
    if (is.null(gfit) && q == 2) {          # fall back to a random intercept
      q <- 1; Zi <- matrix(1, length(y), 1); bi <- bi[1, , drop = FALSE]
      rand_form <- y ~ (1 | group)
      gfit <- fit_glmer(rand_form)
      if (is.null(gfit)) gfit <- fit_glmer(rand_form, nagq = 0)
    }
    if (is.null(gfit)) {
      if (exists("gfit_last")) { gfit <- gfit_last; break }  # keep the last successful fit
      stop("glmer failed in the first GMEXGBoost iteration")
    }
    gfit_last <- gfit
    re <- lme4::ranef(gfit)$group
    bi.old <- bi
    bi <- t(as.matrix(re[lev, seq_len(q), drop = FALSE]))
    colnames(bi) <- lev
    diff.t <- abs(bi.old - bi)
    n.diff <- max(diff.t)
    n.old <- abs(bi.old[which(diff.t == n.diff, arr.ind = TRUE)[1, , drop = FALSE]])
    converged <- is.finite(n.diff / n.old) && (n.diff / n.old < cfg$gmex$toll)
    it <- it + 1
  }
  list(xgb = xgb, glmer = gfit, intercept = unname(lme4::fixef(gfit)[1]),
       bi = bi, n.iteration = it, converged = converged)
}

fit_GMEXGBoost <- function(tr, te, family, cfg) {
  stopifnot(family == "binomial")
  fit <- gmexgb_fit(tr$y, tr$X, tr$id, tr$time, cfg)
  ## glmer intercept + XGBoost + the subject's random effects (b = 0 for unseen subjects);
  ## fit$bi is q x levels with q = 1 (intercept) or 2 (intercept + slope on time)
  eta <- add_re(fit$intercept + as.numeric(stats::predict(fit$xgb, te$X)),
                cbind(1, te$time), te$id, t(fit$bi), colnames(fit$bi))
  imp <- tryCatch({
    ti <- as.data.frame(xgboost::xgb.importance(model = fit$xgb))
    stats::setNames(ti$Gain, ti$Feature)
  }, error = function(e) NULL)
  list(pred = stats::plogis(eta), imp = full_imp(imp, colnames(tr$X)),
       info = list(n_iter = fit$n.iteration, converged = fit$converged))
}

## ---------------------------------------------------------------------------
## Survival competitors
## ---------------------------------------------------------------------------
## XGBoost / XGB-Cox -- package xgboost with EXACTLY the booster settings MEXGB uses
## (same trees, depth, learning rate, min_child_weight, L1/L2, gamma, row and column
## subsampling) but no random effects: the D -> 0 limit of MEXGB. The difference between
## the two is therefore the contribution of the random effects / EM loop alone.
##   gaussian: reg:squarederror   binomial: binary:logistic   cox: survival:cox
fit_XGBoost <- function(tr, te, family, cfg) {
  p <- ncol(tr$X)
  obj <- switch(family, gaussian = "reg:squarederror", binomial = "binary:logistic",
                cox = "survival:cox")
  mcw <- if (family == "gaussian") cfg$boost$minobs else cfg$boost$minobs / 4  # as in MEXGB
  params <- list(objective = obj, eta = cfg$boost$shrinkage,
                 max_depth = cfg$boost$depth, min_child_weight = mcw,
                 lambda = cfg$xgb$lambda, alpha = cfg$xgb$alpha, gamma = cfg$xgb$gamma,
                 subsample = cfg$xgb$subsample, colsample_bytree = cfg$colsample(p),
                 nthread = 1)
  tm <- cfg$xgb$tree_method %||% "auto"           # same tree algorithm as MEXGB
  params$tree_method <- if (identical(tm, "auto")) (if (nrow(tr$X) <= 10000) "exact" else "hist") else tm
  if (family == "gaussian") params$base_score <- mean(tr$y)
  if (family == "binomial") params$base_score <- min(max(mean(tr$y), 1e-4), 1 - 1e-4)
  if (utils::packageVersion("xgboost") >= "2.0.0") params$seed <- 1
  lab <- if (family == "cox") ifelse(tr$status == 1, tr$y, -tr$y) else tr$y
  fit <- xgboost::xgb.train(params = params, data = xgboost::xgb.DMatrix(tr$X, label = lab),
                            nrounds = cfg$boost$ntree, verbose = 0)
  imp <- tryCatch({
    ti <- as.data.frame(xgboost::xgb.importance(model = fit))
    stats::setNames(ti$Gain, ti$Feature)
  }, error = function(e) NULL)
  if (family == "cox") {
    pred <- as.numeric(stats::predict(fit, te$X, outputmargin = TRUE))
    lp_train <- as.numeric(stats::predict(fit, tr$X, outputmargin = TRUE))
  } else {   # identity for gaussian, probability for binomial
    pred <- as.numeric(stats::predict(fit, te$X)); lp_train <- NULL
  }
  list(pred = pred, lp_train = lp_train, imp = full_imp(imp, colnames(tr$X)))
}
fit_XGBCox <- function(tr, te, family, cfg) {
  stopifnot(family == "cox")
  fit_XGBoost(tr, te, family, cfg)
}

## RSF -- package randomForestSRC (random survival forest, independence-assuming).
## Risk score = predicted ensemble mortality; importance = VIMP (negative values set to 0).
fit_RSF <- function(tr, te, family, cfg) {
  stopifnot(family == "cox")
  p <- ncol(tr$X)
  ntree <- if (p > 1000) cfg$real$forest_ntree_hd else cfg$forest$ntree
  dtr <- data.frame(stime = tr$y, sstat = tr$status, tr$X)
  ## one core: randomForestSRC's OpenMP VIMP is otherwise thread-schedule dependent
  ## (non-reproducible feature rankings); predictions are unaffected either way
  old_opt <- options(rf.cores = 1L, mc.cores = 1L); on.exit(options(old_opt), add = TRUE)
  fit <- randomForestSRC::rfsrc(Surv(stime, sstat) ~ ., data = dtr, ntree = ntree,
                                mtry = cfg$forest$mtry_cls(p), nodesize = cfg$rsf$nodesize,
                                importance = cfg$rsf$importance, seed = -1)
  pr <- predict(fit, newdata = data.frame(te$X))$predicted
  vimp <- fit$importance
  list(pred = log(pmax(pr, 1e-8)), lp_train = log(pmax(fit$predicted, 1e-8)),
       imp = full_imp(pmax(vimp, 0), colnames(tr$X)))
}

## Cox-Lasso -- package glmnet (L1-penalised Cox model, no random effects); lambda chosen
## by internal cross-validation with clusters kept together when there are >= 5 clusters.
fit_CoxLasso <- function(tr, te, family, cfg) {
  stopifnot(family == "cox")
  st <- standardise(tr, te)
  ncl <- length(unique(tr$id))
  nf <- cfg$coxlasso$nfolds
  foldid <- if (ncl >= nf) {
    cl <- sample(unique(tr$id)); fm <- rep_len(seq_len(nf), ncl)
    fm[match(tr$id, cl)]
  } else sample(rep_len(seq_len(nf), length(tr$y)))
  fit <- glmnet::cv.glmnet(x = st$Xtr, y = survival::Surv(tr$y, tr$status), family = "cox",
                           foldid = foldid, alpha = 1, nlambda = cfg$coxlasso$nlambda,
                           lambda.min.ratio = cfg$coxlasso$lambda.min.ratio)
  s <- cfg$coxlasso$s
  b <- as.numeric(stats::coef(fit, s = s)); names(b) <- st$keep
  list(pred = as.numeric(stats::predict(fit, newx = st$Xte, s = s, type = "link")),
       lp_train = as.numeric(stats::predict(fit, newx = st$Xtr, s = s, type = "link")),
       imp = full_imp(b, colnames(tr$X)))
}

## CoxFrailty -- package coxme (Cox model with Gaussian random effects; Ripatti-Palmgren),
## random intercept (+ random slope on 'time' when Z has two columns); low dimension only.
fit_CoxFrailty <- function(tr, te, family, cfg) {
  stopifnot(family == "cox")
  vars <- make.names(colnames(tr$X))
  dtr <- data.frame(stime = tr$y, sstat = tr$status, id = factor(tr$id), tm = tr$time, tr$X)
  names(dtr) <- c("stime", "sstat", "id", "tm", vars)
  slope <- ncol(tr$Z) == 2 && isTRUE(cfg$coxfrailty$random_slope)
  re <- if (slope) "(1 + tm | id)" else "(1 | id)"
  form <- stats::as.formula(paste("survival::Surv(stime, sstat) ~",
                                  paste(vars, collapse = " + "), "+", re))
  fit <- quiet(eval(bquote(coxme::coxme(.(form), data = dtr))))
  beta <- coxme::fixef(fit)[vars]
  rf <- coxme::ranef(fit)$id
  rf <- as.matrix(rf); if (ncol(rf) == 1 && is.null(colnames(rf))) colnames(rf) <- "Intercept"
  lp_of <- function(X, id, tm) {
    lp <- as.numeric(X %*% beta)
    k <- match(as.character(id), rownames(rf))
    ok <- !is.na(k)
    lp[ok] <- lp[ok] + rf[k[ok], 1] + (if (ncol(rf) == 2) rf[k[ok], 2] * tm[ok] else 0)
    lp
  }
  tval <- beta / sqrt(diag(as.matrix(stats::vcov(fit))))
  list(pred = lp_of(te$X, te$id, te$time), lp_train = lp_of(tr$X, tr$id, tr$time),
       imp = full_imp(stats::setNames(tval, colnames(tr$X)), colnames(tr$X)))
}

## ---------------------------------------------------------------------------
## applicability rules and dispatcher
## ---------------------------------------------------------------------------
method_applicable <- function(m, p, n_train, cfg, family = NULL) {
  if (m %in% c("LMM", "GLMM-PQL", "CoxFrailty")) return(p <= cfg$p_max_lmm)
  if (m == "glmmLasso") return(p <= cfg$p_max_glmmlasso && p < n_train)
  TRUE
}

FITTERS <- list(MEXGB = fit_MEXGB, MEGB = fit_MEGB, MERF = fit_MERF,
                REEMforest = fit_REEMforest, RF = fit_RF, GBM = fit_GBM,
                GPBoost = fit_GPBoost, mboost = fit_mboost, glmmLasso = fit_glmmLasso,
                LMM = fit_LMM, `GLMM-PQL` = `fit_GLMM-PQL`, GMEXGBoost = fit_GMEXGBoost,
                XGBoost = fit_XGBoost, XGBCox = fit_XGBCox, RSF = fit_RSF, CoxLasso = fit_CoxLasso,
                CoxFrailty = fit_CoxFrailty)

run_method <- function(m, tr, te, family, cfg) {
  f <- FITTERS[[m]]
  if (is.null(f)) stop("unknown method: ", m)
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(f(tr, te, family, cfg), error = function(e) e)
  ct <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error") || is.null(res$pred) || length(res$pred) != length(te$y) ||
      any(!is.finite(res$pred))) {
    msg <- if (inherits(res, "error")) conditionMessage(res) else "invalid predictions"
    return(list(pred = NULL, imp = NULL, time = NA_real_, error = msg))
  }
  res$time <- ct
  res$error <- NA_character_
  res
}
