# =============================================================================
# hali/R/hali_env.R -- bootstrap for every Hali stage
# =============================================================================
# Loads the SHARED analysis code from ../analysis (config.R, R/*.R) -- the same
# files the local scripts use. There is deliberately no cluster-only copy of the
# estimation code: every study, seed and setting comes from analysis/R/registry.R,
# so a split computed on Hali is identical to the same split computed locally.
# =============================================================================

HALI_ROOT <- Sys.getenv("HALI_ROOT", unset = normalizePath(file.path("..", "..")))
ANALYSIS  <- file.path(HALI_ROOT, "analysis")
if (!dir.exists(ANALYSIS))
  stop("analysis directory not found: ", ANALYSIS,
       "\nSet HALI_ROOT to the directory containing MEXGB/, analysis/ and hali/.", call. = FALSE)
setwd(ANALYSIS)                       # analysis/config.R uses paths relative to analysis/

source("config.R")                    # reads the MEXGB_* overrides exported by hali/config.sh
for (f in c("R/methods.R", "R/cv.R", "R/real_data.R", "R/real_survival.R", "R/registry.R"))
  source(f)
suppressPackageStartupMessages(library(MEXGB))

OUT_DIR  <- Sys.getenv("OUT_DIR", unset = file.path(HALI_ROOT, "output_local"))
PREP_DIR <- file.path(OUT_DIR, "prep")
FITS_DIR <- file.path(OUT_DIR, "fits")
LOG_DIR  <- file.path(OUT_DIR, "logs")
REAL_DIR <- file.path(PREP_DIR, "real")
for (d in c(PREP_DIR, FITS_DIR, file.path(LOG_DIR, "tasks"), REAL_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

BLOCKS <- strsplit(Sys.getenv("BLOCKS", paste(ALL_BLOCKS, collapse = ",")), ",")[[1]]
BLOCKS <- intersect(trimws(BLOCKS), ALL_BLOCKS)
SPLIT_THRESHOLD <- as.numeric(Sys.getenv("SPLIT_THRESHOLD", "1800"))

## array tasks are single-core; the local parallel driver is not used on Hali
CFG$N_CORES <- 1L

# ---- small utilities ---------------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
  flush.console()
}

peak_rss_gb <- function() {                    # high-water mark of this process (Linux)
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  l <- grep("^VmHWM:", st, value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", l)) / 1024^2
}

current_task_id <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  id <- if (length(a) >= 1) a[1] else Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (!nzchar(id)) stop("no task id: pass it as an argument or run as a SLURM array", call. = FALSE)
  as.integer(id)
}

in_container <- function() nzchar(Sys.getenv("APPTAINER_CONTAINER")) ||
  nzchar(Sys.getenv("SINGULARITY_CONTAINER")) || file.exists("/.singularity.d")

atomic_save <- function(x, file) {             # never leave a half-written file behind
  tmp <- paste0(file, ".tmp")
  saveRDS(x, tmp)
  file.rename(tmp, file)
}
