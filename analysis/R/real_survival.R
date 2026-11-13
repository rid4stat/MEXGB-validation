## =============================================================================
## Real clustered time-to-event data sets (frailty / mixed-effects survival setting)
##
## Low-dimensional (cluster = patient; subject-level CV, i.e. prediction for new patients)
##   Retinopathy  survival::retinopathy (Huster et al. 1989; DRS): 197 high-risk diabetic
##                patients, one laser-treated and one control eye each; time to visual loss.
##                Covariates: laser type, treated eye, eye side, age at diagnosis, diabetes
##                type, risk score. Random effect: patient intercept (shared frailty).
##   Readmission  frailtypack::readmission (Gonzalez et al. 2005): 403 colorectal-cancer
##                patients, 861 rehospitalisation gap times (recurrent events).
##                Covariates: chemotherapy, sex, Dukes stage, Charlson index, episode number.
##                Random effect: patient intercept (optionally + slope on episode number).
## High-dimensional (cluster = study; observation-level CV stratified by study)
##   BreastMulti  Bioconductor breastCancerVDX / TRANSBIG / MAINZ / UNT: Affymetrix hgu133a
##                expression of ~860 node-negative breast tumours from 4 cohorts with distant-
##                metastasis-free survival (DMFS). Random effect: cohort intercept.
##   OvarianMulti Bioconductor curatedOvarianData (Ganzfried et al. 2013): gene-level
##                expression of ovarian tumours from the studies meeting the package's
##                patient-selection rules, overall survival. Random effect: study intercept.
## Multi-study expression: each gene is z-scored within study; genes are ranked by their
## average within-study SD rank and the CFG$real_surv$n_genes top genes are kept; follow-up
## is administratively censored at CFG$real_surv$max_years.
## =============================================================================

admin_censor <- function(time, status, max_t) {
  status[time > max_t] <- 0
  list(time = pmin(time, max_t), status = status)
}

prep_Retinopathy <- function(cfg) {
  r <- survival::retinopathy
  X <- cbind(laser_argon = as.numeric(r$laser == "argon"),
             treated = as.numeric(r$trt == 1),
             eye_right = as.numeric(r$eye == "right"),
             age = r$age,
             type_adult = as.numeric(r$type == "adult"),
             risk = r$risk)
  d <- make_data(X, r$futime / 12, r$id, rep(0, nrow(r)), status = r$status, Z = "intercept")
  attr(d, "name") <- "Retinopathy"
  attr(d, "description") <- sprintf("Retinopathy: N=%d eyes, n=%d patients, p=%d, events=%d (time in years)",
                                    length(d$y), length(unique(d$id)), ncol(X), sum(d$status))
  d
}

load_readmission <- function(cfg) {
  if (requireNamespace("frailtypack", quietly = TRUE)) {
    e <- new.env(); utils::data("readmission", package = "frailtypack", envir = e)
    return(e$readmission)
  }
  f <- file.path(cfg$dir_data, "readmission.rda")
  if (!file.exists(f)) stop("install 'frailtypack' or place readmission.rda in ", cfg$dir_data)
  e <- new.env(); load(f, envir = e); e$readmission
}

prep_Readmission <- function(cfg) {
  r <- load_readmission(cfg)
  r <- r[r$time > 0, ]
  X <- cbind(chemo = as.numeric(r$chemo == "Treated"),
             female = as.numeric(r$sex == "Female"),
             dukes_C = as.numeric(r$dukes == "C"),
             dukes_D = as.numeric(r$dukes == "D"),
             charlson_1_2 = as.numeric(r$charlson == "1-2"),
             charlson_3 = as.numeric(r$charlson == "3"),
             episode = r$enum)
  tm <- (r$enum - 1) / (max(r$enum) - 1)
  Z <- if (isTRUE(cfg$real_surv$readmission_random_slope)) NULL else "intercept"
  d <- make_data(X, r$time / 365.25, r$id, tm, status = r$event, Z = Z)
  attr(d, "name") <- "Readmission"
  attr(d, "description") <- sprintf("Readmission: N=%d gap times, n=%d patients, p=%d, events=%d (time in years)",
                                    length(d$y), length(unique(d$id)), ncol(X), sum(d$status))
  d
}

## ---- multi-study expression helpers -------------------------------------------------
load_rda_object <- function(file) {
  e <- new.env(); nm <- load(file, envir = e); get(nm[1], envir = e)
}

## keep the n most variable genes (average within-study SD rank), z-score within study
merge_studies <- function(mats, times, status, n_genes, max_years) {
  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) < 50) stop("fewer than 50 genes common to all studies")
  rk <- sapply(mats, function(m) {
    s <- apply(m[common, , drop = FALSE], 1, stats::sd, na.rm = TRUE)
    rank(-s, ties.method = "average")
  })
  keep <- common[order(rowMeans(rk))][seq_len(min(n_genes, length(common)))]
  Xs <- lapply(mats, function(m) {
    z <- t(m[keep, , drop = FALSE])
    z <- scale(z)
    z[!is.finite(z)] <- 0
    z
  })
  X <- do.call(rbind, Xs)
  colnames(X) <- make.names(keep, unique = TRUE)
  id <- rep(seq_along(mats), vapply(mats, ncol, 1L))
  ac <- admin_censor(unlist(times) / 365.25, unlist(status), max_years)
  list(X = X, time = ac$time, status = ac$status, id = id, n_common = length(common))
}

prep_BreastMulti <- function(cfg) {
  rs <- cfg$real_surv
  mats <- list(); times <- list(); stat <- list()
  for (co in rs$breast_cohorts) {
    obj <- tolower(co); pkg <- paste0("breastCancer", co)
    eset <- if (requireNamespace(pkg, quietly = TRUE)) {
      e <- new.env(); utils::data(list = obj, package = pkg, envir = e); get(obj, envir = e)
    } else if (!is.null(rs$breast_local_dir)) {
      f <- list.files(rs$breast_local_dir, pattern = paste0("^", obj, "\\.rda$"),
                      recursive = TRUE, full.names = TRUE)
      if (!length(f)) stop("no ", obj, ".rda in ", rs$breast_local_dir)
      load_rda_object(f[1])
    } else stop("install Bioconductor package ", pkg, " (BiocManager::install('", pkg, "'))")
    pd <- Biobase::pData(eset)
    ok <- !is.na(pd$t.dmfs) & !is.na(pd$e.dmfs) & pd$t.dmfs > 0
    mats[[co]] <- Biobase::exprs(eset)[, ok, drop = FALSE]
    times[[co]] <- pd$t.dmfs[ok]; stat[[co]] <- pd$e.dmfs[ok]
  }
  m <- merge_studies(mats, times, stat, rs$n_genes, rs$max_years)
  d <- make_data(m$X, m$time, m$id, rep(0, length(m$time)), status = m$status, Z = "intercept")
  d$cluster_names <- names(mats)
  attr(d, "name") <- "BreastMulti"
  attr(d, "description") <- sprintf("BreastMulti: N=%d patients, %d cohorts (%s), p=%d of %d common probes, events=%d (DMFS, years, censored at %d y)",
                                    length(d$y), length(mats),
                                    paste(sprintf("%s=%d", names(mats), vapply(mats, ncol, 1L)), collapse = ", "),
                                    ncol(m$X), m$n_common, sum(d$status), rs$max_years)
  d
}

prep_OvarianMulti <- function(cfg) {
  rs <- cfg$real_surv
  have_pkg <- requireNamespace("curatedOvarianData", quietly = TRUE)
  if (have_pkg) {
    items <- utils::data(package = "curatedOvarianData")$results[, "Item"]
    items <- grep("_eset$", items, value = TRUE)
    get_eset <- function(nm) { e <- new.env(); utils::data(list = nm, package = "curatedOvarianData", envir = e); get(nm, envir = e) }
    conf <- system.file("extdata", "patientselection.config", package = "curatedOvarianData")
  } else if (!is.null(rs$ovarian_local_dir)) {
    files <- list.files(rs$ovarian_local_dir, pattern = "_eset\\.rda$", full.names = TRUE)
    items <- sub("\\.rda$", "", basename(files))
    get_eset <- function(nm) load_rda_object(files[match(nm, items)])
    conf <- file.path(rs$ovarian_local_dir, "patientselection.config")
  } else stop("install Bioconductor package curatedOvarianData (BiocManager::install('curatedOvarianData'))")

  ## samples flagged as duplicates / outliers by the package's own selection config
  remove_samples <- character(0)
  if (nzchar(conf) && file.exists(conf)) {
    e <- new.env()
    try(sys.source(conf, envir = e), silent = TRUE)
    if (exists("remove.samples", envir = e)) remove_samples <- get("remove.samples", envir = e)
  }
  ## mRNA expression studies only (drop TCGA RNA-seq duplicate and miRNA arrays)
  items <- items[!grepl("RNASeq|mirna", items)]

  mats <- list(); times <- list(); stat <- list()
  for (nm in items) {
    es <- get_eset(nm)
    pd <- Biobase::pData(es)
    if (!all(c("days_to_death", "vital_status") %in% names(pd))) next
    ok <- !is.na(pd$days_to_death) & !is.na(pd$vital_status) & pd$days_to_death > 0
    if ("sample_type" %in% names(pd)) ok <- ok & grepl("^tumor$", pd$sample_type)
    ok <- ok & !(paste0(nm, ":", Biobase::sampleNames(es)) %in% remove_samples)
    ev <- as.numeric(pd$vital_status == "deceased")
    if (sum(ok) < rs$ovarian_min_n || sum(ev[ok]) < rs$ovarian_min_events) next
    m <- Biobase::exprs(es)[, ok, drop = FALSE]
    m <- m[!grepl("///", rownames(m), fixed = TRUE), , drop = FALSE]
    if (nrow(m) < 1000) next
    mats[[nm]] <- m; times[[nm]] <- pd$days_to_death[ok]; stat[[nm]] <- ev[ok]
  }
  if (length(mats) < 2) stop("fewer than two eligible ovarian studies")
  mm <- merge_studies(mats, times, stat, rs$n_genes, rs$max_years)
  d <- make_data(mm$X, mm$time, mm$id, rep(0, length(mm$time)), status = mm$status, Z = "intercept")
  d$cluster_names <- names(mats)
  attr(d, "name") <- "OvarianMulti"
  attr(d, "description") <- sprintf("OvarianMulti: N=%d patients, %d studies, p=%d of %d common genes, events=%d (OS, years, censored at %d y)",
                                    length(d$y), length(mats), ncol(mm$X), mm$n_common,
                                    sum(d$status), rs$max_years)
  d
}

## ---- Insemination (low-dimensional; cluster = herd) ----------------------------------------
## parfm::insem (Duchateau and Janssen 2008): time (days) from calving to first insemination of
## 10513 dairy cows in 181 herds; event = insemination, censored otherwise (culled, sold, end of
## study). Herds with >= min_cows cows. Predictors (p = 4): milk ureum and protein concentration,
## parity, heifer (first parity). Random intercept per herd (shared frailty).
prep_Insemination <- function(cfg) {
  ic <- cfg$real_surv$insem
  e <- new.env(); load(ic$file, envir = e); x <- get(ls(e)[1], envir = e)
  h <- table(x$Herd); x <- x[x$Herd %in% names(h)[h >= ic$min_cows], , drop = FALSE]
  x <- x[x$Time > 0, , drop = FALSE]
  X <- cbind(ureum = x$Ureum, protein = x$Protein, parity = x$Parity, heifer = x$Heifer)
  d <- make_data(X, x$Time / 30.4375, as.integer(factor(x$Herd)), rep(0, nrow(x)),
                 status = x$Status, Z = "intercept")
  attr(d, "name") <- "Insemination"
  attr(d, "description") <- sprintf("Insemination: N=%d cows, n=%d herds (>= %d cows), p=%d, events=%d (time in months)",
                                    length(d$y), length(unique(d$id)), ic$min_cows, ncol(X), sum(d$status))
  d
}

## ---- OvarianDebulking (high-dimensional, binary; cluster = study) ----------------------------
## curatedOvarianData: tumour expression of the mRNA studies that record the surgical debulking
## status (optimal / suboptimal; Riester et al. 2014 predicted it across these studies);
## outcome = suboptimal debulking (1). Studies with >= ovarian_min_n patients and >= 1000 genes;
## duplicate / outlier samples of the package's selection config removed. Genes standardised
## within study, the n_genes most variable (average within-study SD rank) common genes kept.
## Random intercept per study.
prep_OvarianDebulking <- function(cfg) {
  rs <- cfg$real_surv
  if (!requireNamespace("curatedOvarianData", quietly = TRUE) && is.null(rs$ovarian_local_dir))
    stop("install Bioconductor package curatedOvarianData (BiocManager::install('curatedOvarianData'))")
  if (requireNamespace("curatedOvarianData", quietly = TRUE)) {
    items <- grep("_eset$", utils::data(package = "curatedOvarianData")$results[, "Item"], value = TRUE)
    get_eset <- function(nm) { e <- new.env(); utils::data(list = nm, package = "curatedOvarianData", envir = e); get(nm, envir = e) }
    conf <- system.file("extdata", "patientselection.config", package = "curatedOvarianData")
  } else {
    files <- list.files(rs$ovarian_local_dir, pattern = "_eset\\.rda$", full.names = TRUE)
    items <- sub("\\.rda$", "", basename(files))
    get_eset <- function(nm) load_rda_object(files[match(nm, items)])
    conf <- file.path(rs$ovarian_local_dir, "patientselection.config")
  }
  remove_samples <- character(0)
  if (nzchar(conf) && file.exists(conf)) {
    e <- new.env(); try(sys.source(conf, envir = e), silent = TRUE)
    if (exists("remove.samples", envir = e)) remove_samples <- get("remove.samples", envir = e)
  }
  items <- items[!grepl("RNASeq|mirna", items)]
  mats <- list(); ys <- list()
  for (nm in items) {
    es <- get_eset(nm); pd <- Biobase::pData(es)
    if (!"debulking" %in% names(pd)) next
    ok <- pd$debulking %in% c("optimal", "suboptimal")
    if ("sample_type" %in% names(pd)) ok <- ok & grepl("^tumor$", pd$sample_type)
    ok <- ok & !(paste0(nm, ":", Biobase::sampleNames(es)) %in% remove_samples)
    if (sum(ok) < rs$ovarian_min_n) next
    m <- Biobase::exprs(es)[, ok, drop = FALSE]
    m <- m[!grepl("///", rownames(m), fixed = TRUE), , drop = FALSE]
    if (nrow(m) < 1000) next
    mats[[nm]] <- m; ys[[nm]] <- as.numeric(pd$debulking[ok] == "suboptimal")
  }
  if (length(mats) < 2) stop("fewer than two eligible ovarian studies with debulking status")
  mm <- merge_studies(mats, lapply(ys, function(y) rep(1, length(y))), lapply(ys, function(y) rep(0, length(y))),
                      rs$n_genes, Inf)
  y <- unlist(ys)
  d <- make_data(mm$X, y, mm$id, rep(0, length(y)), Z = "intercept")
  d$y_binary <- y
  d$cluster_names <- names(mats)
  attr(d, "name") <- "OvarianDebulking"
  attr(d, "description") <- sprintf("OvarianDebulking: N=%d patients, %d studies, p=%d of %d common genes, suboptimal debulking %.1f%%",
                                    length(y), length(mats), ncol(mm$X), mm$n_common, 100 * mean(y))
  d
}

load_real_surv <- function(name, cfg) {
  switch(name,
         Insemination = prep_Insemination(cfg),
         Retinopathy = prep_Retinopathy(cfg),
         Readmission = prep_Readmission(cfg),
         BreastMulti = prep_BreastMulti(cfg),
         OvarianMulti = prep_OvarianMulti(cfg),
         stop("unknown survival data set ", name))
}
