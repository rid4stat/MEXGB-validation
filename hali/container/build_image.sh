#!/usr/bin/env bash
# =============================================================================
# hali/container/build_image.sh -- build the analysis image
# =============================================================================
#   cd hali/container
#   ./build_image.sh                      # uses ../../MEXGB_0.5.0.tar.gz
#   ./build_image.sh /path/to/MEXGB_0.5.0.tar.gz
#
#   A. --fakeroot build from mexgb.def: everything baked in. Preferred.
#   B. If fakeroot is disabled, pull the base image and install the R packages into
#      a HOST library bind-mounted at run time (no privilege needed).
#
# Run on the LOGIN node -- compute nodes usually have no outbound network.
# gpboost and xgboost compile large C++ sources: allow ~30-60 min and >= 8 GB RAM.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
TARBALL="${1:-$HERE/../../MEXGB_0.5.0.tar.gz}"
[[ -f "$TARBALL" ]] || { echo "MEXGB tarball not found: $TARBALL" >&2; exit 1; }
MEGB_TARBALL="${2:-$HERE/../../MEGB_0.4.tar.gz}"   # MEGB 0.4 (cross-fitted + converged E-step), replaces CRAN 0.2
[[ -f "$MEGB_TARBALL" ]] || { echo "MEGB tarball not found: $MEGB_TARBALL" >&2; exit 1; }
command -v apptainer >/dev/null 2>&1 || { echo "apptainer not on PATH. Run: module load apptainer/1.4.0-rc.1" >&2; exit 1; }

export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRATCH:-/tmp}/${USER:-$(id -un)}/apptainer_cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${SCRATCH:-/tmp}/${USER:-$(id -un)}/apptainer_tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
echo "cache: $APPTAINER_CACHEDIR"; echo "tmp  : $APPTAINER_TMPDIR"

SIF="$HERE/mexgb.sif"
cp -f "$TARBALL" "$HERE/MEXGB_0.5.0.tar.gz"
cp -f "$MEGB_TARBALL" "$HERE/MEGB_0.4.tar.gz"

echo; echo "=== Route A: --fakeroot build (self-contained image) ==="
if apptainer build --fakeroot "$SIF" mexgb.def; then
  echo "baked" > "$HERE/image_kind"
  rm -f "$HERE/MEXGB_0.5.0.tar.gz" "$HERE/MEGB_0.4.tar.gz"
  echo; echo "BUILT: $SIF"
  apptainer exec "$SIF" Rscript -e 'library(MEXGB); v1 <- as.character(packageVersion("MEXGB")); v2 <- as.character(packageVersion("MEGB")); cf <- all(c("crossfit","init","dispersion","inner","calibrate","linear","ntree_select","screen","estep","zfixed","honest","honest_d") %in% names(formals(MEXGB::MEXGB))) && "inner" %in% names(formals(MEGB::MEGB)); cat(R.version.string, "\n"); cat("MEXGB", v1, "| MEGB", v2, "(crossfit:", if (cf) "yes" else "NO", ")\n"); if (!cf || package_version(v1) < "0.5.0" || package_version(v2) < "0.4") stop("STALE BUILD: MEXGB >= 0.5.0 and MEGB >= 0.4 are required (re-copy the v0.5 bundle)")'
  exit 0
fi

echo; echo "=== fakeroot build failed; Route B: base image + host-side R library ==="
apptainer pull --force "$SIF" docker://rocker/r-ver:4.4.3
HOSTLIB="${R_LIBS_USER:-$HOME/R/hali/4.4.3}"
mkdir -p "$HOSTLIB"
echo "host library: $HOSTLIB"
RX=(apptainer exec --bind "$HOSTLIB:$HOSTLIB" --bind "$HERE:$HERE" --env R_LIBS_USER="$HOSTLIB" "$SIF")

"${RX[@]}" Rscript -e "
  .libPaths(c('$HOSTLIB', .libPaths()))
  install.packages(c('xgboost','MEGB','gbm','randomForest','gpboost','mboost','glmmLasso',
    'nlme','MASS','lme4','survival','randomForestSRC','glmnet','coxme','ggplot2',
    'data.table','jsonlite','future','future.apply','remotes','BiocManager','latex2exp',
    'mvtnorm','rpart'), lib='$HOSTLIB', Ncpus=max(1L, parallel::detectCores()))
  ok <- tryCatch({remotes::install_github('cran/LongituRF', lib='$HOSTLIB', upgrade='never'); TRUE},
                 error=function(e) FALSE)
  if (!ok) install.packages('https://cran.r-project.org/src/contrib/Archive/LongituRF/LongituRF_0.9.tar.gz',
                            repos=NULL, type='source', lib='$HOSTLIB')
  BiocManager::install(c('Biobase','breastCancerVDX','breastCancerTRANSBIG','breastCancerMAINZ',
                         'breastCancerUNT','curatedOvarianData'), lib='$HOSTLIB', update=FALSE, ask=FALSE)
"
"${RX[@]}" R CMD INSTALL --library="$HOSTLIB" "$HERE/MEXGB_0.5.0.tar.gz"
"${RX[@]}" R CMD INSTALL --library="$HOSTLIB" "$HERE/MEGB_0.4.tar.gz"
echo "hostlib" > "$HERE/image_kind"
"${RX[@]}" Rscript -e "
  .libPaths(c('$HOSTLIB', .libPaths()))
  req <- c('MEXGB','xgboost','MEGB','gbm','randomForest','LongituRF','gpboost','mboost',
           'glmmLasso','lme4','survival','randomForestSRC','glmnet','coxme','ggplot2')
  miss <- req[!req %in% rownames(installed.packages())]
  if (length(miss)) stop('MISSING: ', paste(miss, collapse=', '))
  cat(R.version.string, '\n')
  v1 <- as.character(packageVersion('MEXGB')); v2 <- as.character(packageVersion('MEGB'))
  cf <- all(c('crossfit','dispersion','inner','calibrate','linear','ntree_select','screen','estep','zfixed','honest','honest_d') %in% names(formals(MEXGB::MEXGB))) && 'inner' %in% names(formals(MEGB::MEGB)) && 'crossfit' %in% names(formals(MEGB::MEGB))
  cat('MEXGB', v1, '| MEGB', v2, '(crossfit:', if (cf) 'yes' else 'NO', ')\n')
  if (!cf || package_version(v1) < '0.5.0' || package_version(v2) < '0.4')
    stop('STALE BUILD: MEXGB >= 0.5.0 and MEGB >= 0.4 are required (re-copy the v0.5 bundle)')
  opt <- c('Biobase','breastCancerVDX','curatedOvarianData')
  mo <- opt[!opt %in% rownames(installed.packages())]
  if (length(mo)) cat('NOTE: missing (multi-study survival data skipped):', paste(mo, collapse=', '), '\n')
"
rm -f "$HERE/MEXGB_0.5.0.tar.gz" "$HERE/MEGB_0.4.tar.gz"
echo; echo "BUILT (host-library mode): $SIF"
echo "Keep R_LIBS_USER=$HOSTLIB in hali/config.sh"
