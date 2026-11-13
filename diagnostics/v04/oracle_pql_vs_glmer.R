setwd("../../analysis"); source("config.R"); suppressPackageStartupMessages(library(MEXGB)); S <- CFG$sim
for (dgp in c("linear","nonlinear")) for (sd in 1:3) {
cal <- calibrateSim("binomial", target = S$target_auc, prevalence = S$prevalence, linear = dgp=="linear", time_points = 20, random_sd_intercept = S$tau0, random_sd_slope = S$tau1, pilot_n = 300, pilot_seed = 2026)
set.seed(sd); s <- simLong(n = 50, p = 6, time_points = 20, linear = dgp=="linear", family = "binomial", calibration = cal)
Z <- as.matrix(s[,4:5]); bt <- attr(s,"b"); fix <- attr(s,"eta") - rowSums(Z * bt[s$id,])
il <- split(seq_along(s$id), s$id); D <- diag(2); b <- matrix(0, 50, 2); tr <- c()
for (it in 1:40) { eta <- fix + rowSums(Z * b[s$id,]); mu <- plogis(eta); w <- pmax(mu*(1-mu),1e-6); yt <- eta + (s$Y-mu)/w
  u <- MEXGB:::pql_em_additive(yt - fix, Z, D, 0, w, il); b <- u$b; D <- u$D; if (it %in% c(1,5,10,20,40)) tr <- c(tr, sprintf("it%d %.2f", it, D[2,2])) }
m <- suppressMessages(lme4::glmer(Y ~ (1 + time | id), data = s, family = binomial, offset = fix))
cat(dgp, sd, "| PQL D22:", paste(tr, collapse=" "), "| glmer D22 = ", round(lme4::VarCorr(m)$id[2,2],2), "| var(true b1) =", round(var(bt[,2]),2), "\n") }
