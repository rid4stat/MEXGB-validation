#!/usr/bin/env Rscript
# =============================================================================
# hali/R/04_outstanding.R -- write the ids of tasks that are not yet complete
# =============================================================================
#   Rscript 04_outstanding.R <out_file>
# Used by './submit.sh fit' on a resubmission, so only outstanding tasks are queued.
# =============================================================================
source(file.path(Sys.getenv("HALI_ROOT", "../.."), "hali", "R", "hali_env.R"))
out <- commandArgs(trailingOnly = TRUE)[1]
tt <- readRDS(file.path(PREP_DIR, "tasks.rds"))
stems <- task_stem(FITS_DIR, tt$study_id, tt$rep, tt$fold, tt$unit)
todo <- tt$task_id[!file.exists(paste0(stems, ".done"))]
writeLines(as.character(todo), out)
cat(length(todo), "\n")
