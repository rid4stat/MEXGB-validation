## =============================================================================
## MEXGB paper 1 -- install everything needed (run once)
## Run from the 'analysis' folder:   Rscript 00_install.R
## =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("xgboost", "MEGB", "gbm", "randomForest", "gpboost", "mboost", "glmmLasso",
          "nlme", "MASS", "lme4", "ggplot2", "future", "future.apply", "remotes",
          "BiocManager", "data.table", "jsonlite",
          ## survival
          "survival", "randomForestSRC", "glmnet", "coxme")
miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)

## LongituRF (MERF, REEMforest) is archived on CRAN: install from the CRAN archive mirror
if (!requireNamespace("LongituRF", quietly = TRUE)) {
  for (d in c("latex2exp", "mvtnorm", "rpart"))
    if (!requireNamespace(d, quietly = TRUE)) install.packages(d)
  ok <- tryCatch({ remotes::install_github("cran/LongituRF"); TRUE }, error = function(e) FALSE)
  if (!ok) install.packages("https://cran.r-project.org/src/contrib/Archive/LongituRF/LongituRF_0.9.tar.gz",
                            repos = NULL, type = "source")
}

## Sprockett data: shipped as data/SprockettData.csv (no download needed)
## Readmission data: shipped as data/readmission.rda (from CRAN package frailtypack; installing
## frailtypack is optional)
## Multi-study gene-expression survival data (Bioconductor data packages)
for (b in c("Biobase", "breastCancerVDX", "breastCancerTRANSBIG", "breastCancerMAINZ",
            "breastCancerUNT", "curatedOvarianData"))
  if (!requireNamespace(b, quietly = TRUE)) BiocManager::install(b, update = FALSE, ask = FALSE)

## the MEXGB package itself (source folder shipped next to 'analysis')
install.packages("../MEXGB", repos = NULL, type = "source")
## MEGB 0.3 (cross-fitted E-step; same fix as MEXGB 0.3) -- replaces CRAN's MEGB 0.2
install.packages("../MEGB", repos = NULL, type = "source")

pk <- c("MEXGB", "xgboost", "MEGB", "gbm", "randomForest", "LongituRF", "gpboost", "mboost",
        "glmmLasso", "nlme", "MASS", "lme4", "ggplot2", "future.apply", "survival",
        "randomForestSRC", "glmnet", "coxme", "Biobase", "breastCancerVDX",
        "breastCancerTRANSBIG", "breastCancerMAINZ", "breastCancerUNT", "curatedOvarianData")
print(data.frame(package = pk,
                 version = vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                                           error = function(e) "NOT INSTALLED"), "")))
