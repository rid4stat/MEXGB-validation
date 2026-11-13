## Diagnostic harness: run selected methods / MEXGB variants on the first splits of
## benchmark studies, with exactly the data and splits of the Hali campaign.
## Usage: Rscript harness.R <block> <label_regex> <reps> <folds> <methods,comma> <out.csv>
args <- commandArgs(TRUE)
setwd("/home/claude/out/v05/MEXGB_paper1/analysis")
suppressPackageStartupMessages({source("config.R"); source("R/methods.R"); source("R/cv.R")
  source("R/real_data.R"); source("R/real_survival.R"); source("R/registry.R"); library(MEXGB)}); if (nzchar(Sys.getenv("MEXGB_DEV"))) pkgload::load_all("/home/claude/out/v05/MEXGB_paper1/MEXGB", quiet = TRUE)
source("../diagnostics/v05/variants.R")
block <- args[1]; rx <- args[2]; reps <- eval(parse(text = args[3])); folds <- eval(parse(text = args[4]))
methods <- strsplit(args[5], ",")[[1]]; out <- args[6]
reg <- switch(block, csim = registry_csim(CFG), bsim = registry_bsim(CFG), ssim = registry_ssim(CFG),
              real = registry_real(CFG), sreal = registry_sreal(CFG))
reg <- Filter(function(st) grepl(rx, st$label), reg)
res <- list()
for (st in reg) for (r in reps) {
  d <- study_data_fun(st, CFG)(r)
  fl <- make_folds(d, CFG$K, seed = st$seed + 1000 * r, unit = st$cv_unit,
                   n_test = CFG$n_test_per_subject %||% 3, test_frac = CFG$test_frac %||% 0)
  for (k in folds) {
    x <- run_fold(d, fl[[k]], methods, st$family, CFG, active = st$active,
                  seed = st$seed + 1000 * r + k, tags = c(st$tags, list(label = st$label, rep = r, fold = k)))
    x$top_features <- NULL
    res[[length(res) + 1]] <- x
    cat(sprintf("%s rep %d fold %d done\n", st$label, r, k)); flush.console()
    write.csv(rbind_fill(res), out, row.names = FALSE)
  }
}
