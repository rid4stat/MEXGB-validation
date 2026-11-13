## =============================================================================
## MEXGB -- descriptive summary of the real data sets (paper data-description table)
## Outputs: results/tables/tab_real_data.csv, tab_real_data.tex
## Run from the 'analysis' folder:   Rscript 08_data_description.R
## =============================================================================
source("config.R")
source("R/methods.R"); source("R/cv.R"); source("R/real_data.R")
source("R/real_survival.R"); source("R/registry.R"); source("R/describe_real.R")

tab <- describe_real_data(CFG)
if (is.null(tab)) stop("no real data set could be loaded")
utils::write.csv(tab, file.path(CFG$dir_tables, "tab_real_data.csv"), row.names = FALSE)
print(tab[, c("dataset", "N", "clusters", "per_cluster", "p", "outcome", "events",
              "binary_prevalence", "test_per_split")], right = FALSE)

esc <- function(x) gsub("%", "\\\\%", gsub("_", "\\\\_", ifelse(is.na(x), "--", x)))
lines <- with(tab, sprintf("%s & %d & %d & %s & %d & %s & %s \\\\", dataset, N, clusters,
                           per_cluster, p, esc(random_effects),
                           esc(ifelse(is.na(events), binary_prevalence, events))))
writeLines(c("\\begin{table}[htbp]\\centering\\small",
             "\\caption{Real data sets. $N$: observations; $n$: clusters (subjects, plots, patients or studies); $n_i$: observations per cluster, median [range]; $p$: predictors. Last column: event rate of the binary outcome (continuous data sets) or number (\\%) of events (survival data sets). Each test set holds max(15\\%, 3) observations of every cluster, at least one observation of every cluster being kept for training.}",
             "\\label{tab:realdata}",
             "\\begin{tabular}{lrrlrll}", "\\toprule",
             "Data set & $N$ & $n$ & $n_i$ & $p$ & Random effects & Binary rate / events\\\\", "\\midrule",
             lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
           file.path(CFG$dir_tables, "tab_real_data.tex"))
message("written ", file.path(CFG$dir_tables, "tab_real_data.tex"))
