fs <- list.files("raw", full.names = TRUE, "^bench_.*\\.rds$"); d <- do.call(rbind, lapply(fs, function(f) cbind(readRDS(f), src=basename(f))))
d$cell <- paste(c(gaussian="cont", binomial="bin", cox="surv")[as.character(d$fam)], d$dgp, d$p)
d <- d[!grepl("_phi", d$src), ]          # sensitivity runs with an estimated binary dispersion
d$method[d$method == "XGBCox"] <- "XGBoost_new"   # survival: the XGB-Cox baseline
ord <- c("MEXGB_old","MEXGB_new","MEGB_old","MEGB_new","XGBoost_old","XGBoost_new","GBM")
nm <- tapply(d$metric, list(d$cell, d$method), length)[, ord]; print(nm)
## common splits only (the main files, all methods present)
main <- d[!grepl("_b\\.rds|_c\\.rds", d$src), ]
m <- round(tapply(main$metric, list(main$cell, main$method), mean)[, ord], 3)
## nonlinear 2000: 6 splits for the 'new' methods
x <- d[d$cell=="cont nonlinear 2000" & d$method %in% c("MEXGB_new","MEGB_new","XGBoost_new","GBM"),]
cat("\nnonlinear p=2000, 6 splits (3 data sets):\n"); print(round(tapply(x$metric, x$method, mean),2))
wins <- function(a,b) { y <- main[main$method %in% c(a,b),]; w <- reshape(y[,c("cell","src","rep","fold","method","metric")], idvar=c("cell","src","rep","fold"), timevar="method", direction="wide")
  s <- ifelse(grepl("bin|surv", w$cell), 1, -1); tapply(sign(s*(w[[paste0("metric.",a)]]-w[[paste0("metric.",b)]]))>0, w$cell, function(z) sprintf("%d/%d", sum(z), length(z))) }
out <- data.frame(cell=rownames(m), m, MEXGBnew_beats_XGB=wins("MEXGB_new","XGBoost_new")[rownames(m)], MEXGBnew_beats_GBM=wins("MEXGB_new","GBM")[rownames(m)], check.names=FALSE)
print(out, row.names=FALSE); write.csv(out, "benchmark_before_after.csv", row.names=FALSE)
