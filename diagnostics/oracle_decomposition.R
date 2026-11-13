## Oracle decomposition of the test MSE (XGBoost vs MEXGB vs oracle random effects), nonlinear design.
## Rscript oracle_decomposition.R <p> <seeds, e.g. 1,2,3>
suppressMessages(library(MEXGB)); a <- commandArgs(TRUE); p <- as.integer(a[1]); seeds <- as.integer(strsplit(a[2],",")[[1]])
mse <- function(a,b) mean((a-b)^2)
xg <- function(X,y) xgboost::xgb.train(list(objective="reg:squarederror", eta=.05, max_depth=3, min_child_weight=5, nthread=1, tree_method="exact", base_score=mean(y)), xgboost::xgb.DMatrix(X,label=y), 500, verbose=0)
for (s in seeds) {
set.seed(s); d <- simLong(n=20, p=p, time_points=10, linear=FALSE)
X <- as.matrix(d[, grep("^Var", names(d))]); Z <- as.matrix(d[,4:5]); b <- attr(d,"b"); zb <- rowSums(Z*b[d$id,]); f <- attr(d,"f")
te <- unlist(lapply(split(seq_along(d$id), d$id), function(r) sample(r,3))); tr <- setdiff(seq_along(d$id), te)
m1 <- xg(X[tr,], d$Y[tr]); f1 <- predict(m1, X[te,])
m2 <- xg(X[tr,], d$Y[tr]-zb[tr]); f2 <- predict(m2, X[te,])
fit <- MEXGB(X[tr,], d$Y[tr], d$id[tr], Z[tr,], time=d$time[tr], ntree=500, interaction.depth=3, verbose=FALSE)
fm <- MEXGB:::fixedMEXGB(fit, X[te,]); bh <- fit$random_effects; zbh <- rowSums(Z[te,]*bh[match(d$id[te], fit$id_btilde),])
cat(sprintf("seed %d p %d | XGB %.2f | XGB+trueZb %.2f | oracleMEXGB(f on Y-Zb)+Zb %.2f | MEXGB %.2f | MEXGB f + true Zb %.2f | f-err: XGB %.2f MEXGBf %.2f | Zb-err MEXGB %.2f (var Zb %.2f) | D22 %.2f D11 %.2f sig %.2f\n",
  s, p, mse(d$Y[te],f1), mse(d$Y[te],f1+zb[te]), mse(d$Y[te],f2+zb[te]), mse(d$Y[te],fm+zbh), mse(d$Y[te],fm+zb[te]),
  mse(f[te],f1-mean(f1-f[te])), mse(f[te],fm), mse(zb[te],zbh), var(zb[te]), fit$var_random_effects[2,2], fit$var_random_effects[1,1], fit$sigma))
}
