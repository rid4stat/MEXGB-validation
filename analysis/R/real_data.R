## =============================================================================
## Real longitudinal data sets (paper Section 7)
##   BodyWeight : nlme::BodyWeight  (16 rats, 11 occasions, 3 diets)
##   Soybean    : nlme::Soybean     (48 plots, 8-10 occasions)
##   Sprockett  : Tsimane gut microbiome cohort (Sprockett et al. 2020, Nat Commun 11:3772),
##                16S rRNA ASV profiles, read from data/SprockettData.csv
## Every loader returns make_data(X, y, id, time) with time rescaled to [0, 1], plus
## attributes 'name' and 'description'.
## =============================================================================

rescale01 <- function(t) (t - min(t)) / (max(t) - min(t))

prep_BodyWeight <- function(cfg) {
  bw <- as.data.frame(nlme::BodyWeight)
  id <- as.integer(as.character(bw$Rat))
  X <- cbind(Time = bw$Time,
             Diet2 = as.numeric(bw$Diet == "2"),
             Diet3 = as.numeric(bw$Diet == "3"))
  d <- make_data(X, bw$weight, id, rescale01(bw$Time))
  attr(d, "name") <- "BodyWeight"
  attr(d, "description") <- sprintf("BodyWeight: N=%d, n=%d rats, p=%d", length(d$y),
                                    length(unique(id)), ncol(X))
  d
}

prep_Soybean <- function(cfg) {
  sb <- as.data.frame(nlme::Soybean)
  id <- as.integer(factor(as.character(sb$Plot)))
  X <- cbind(Time = sb$Time,
             VarietyP = as.numeric(sb$Variety == "P"),
             Year1989 = as.numeric(sb$Year == "1989"),
             Year1990 = as.numeric(sb$Year == "1990"))
  d <- make_data(X, sb$weight, id, rescale01(sb$Time))
  attr(d, "name") <- "Soybean"
  attr(d, "description") <- sprintf("Soybean: N=%d, n=%d plots, p=%d", length(d$y),
                                    length(unique(id)), ncol(X))
  d
}

## ---- Sprockett microbiome cohort ---------------------------------------------------
## Source: data/SprockettData.csv (one row per sample): sample metadata, 2087 ASV columns
## ("ASV_<n>_<taxon>", stored as log(1 + relative abundance)), and the derived columns
## 'shannon' (Shannon diversity of the sample) and 'shannon_binary' (median split).
## Cohort: feces samples of subjects with >= min_visits samples -> 193 subjects, 1815 samples.
## Features: ASV columns centred within sample (log-ratio to the sample mean), then the
## lowest-variance decile of ASVs removed (2087 -> 1878).
## Longitudinal time (random slope):
##   time_var = "followup": months since the subject's first sample (Age_Months), / max -> [0,1]
##   time_var = "visit"   : within-subject visit index ordered by age, rescaled to [0,1]
prep_Sprockett <- function(cfg) {
  sc <- cfg$real$sprockett
  if (!file.exists(sc$file)) stop("Sprockett file not found: ", sc$file)
  raw <- utils::read.csv(sc$file, check.names = FALSE, stringsAsFactors = FALSE)
  raw <- raw[raw$Sample_Type == sc$sample_type, , drop = FALSE]
  nv <- table(raw$Subject_ID)
  raw <- raw[raw$Subject_ID %in% names(nv)[nv >= sc$min_visits], , drop = FALSE]
  raw <- raw[order(raw$Subject_ID, raw$Age_Months, raw$collection_date), , drop = FALSE]

  subj <- as.character(raw$Subject_ID)
  visit <- stats::ave(seq_along(subj), subj, FUN = seq_along)
  time <- if (identical(sc$time_var, "visit")) {
    (visit - 1) / (max(visit) - 1)
  } else {
    fu <- raw$Age_Months - stats::ave(raw$Age_Months, subj, FUN = min)
    fu / max(fu)
  }

  asv <- grep("^ASV_", names(raw))
  X <- as.matrix(raw[, asv])
  if (isTRUE(sc$centre)) X <- X - rowMeans(X)
  v <- apply(X, 2, stats::var)
  n_drop <- floor(sc$drop_low_var_frac * ncol(X))
  keep <- order(v, decreasing = TRUE)[seq_len(ncol(X) - n_drop)]
  X <- X[, sort(keep), drop = FALSE]
  colnames(X) <- make.names(colnames(X), unique = TRUE)

  id <- as.integer(factor(subj))
  d <- make_data(X, raw$shannon, id, time)
  d$y_binary <- as.numeric(raw$shannon_binary)
  attr(d, "name") <- "Sprockett"
  attr(d, "description") <- sprintf("Sprockett: N=%d samples, n=%d subjects, p=%d ASVs, time=%s",
                                    length(d$y), length(unique(id)), ncol(X), sc$time_var)
  d
}

## ---- PBC (low-dimensional, binary) ---------------------------------------------------------
## Mayo Clinic trial in primary biliary cirrhosis (JM::pbc2; Murtaugh et al. 1994): repeated
## clinical visits of 312 patients randomised to D-penicillamine or placebo. Patients with >=
## min_visits (10) visits; binary outcome = hepatomegaly (enlarged liver) at the visit; visits
## with a missing outcome or predictor are dropped. Predictors (p = 13): treatment, age, sex,
## log serum bilirubin, albumin, log alkaline phosphatase, log SGOT, platelets, prothrombin
## time, histologic stage, oedema (2 indicators) and years since enrolment. Random intercept +
## slope on time (years since enrolment, rescaled to [0, 1]).
prep_PBC <- function(cfg) {
  pc <- cfg$real$pbc
  e <- new.env(); load(pc$file, envir = e); pb <- e$pbc2
  vars <- c(pc$outcome, "drug", "age", "sex", "year", "serBilir", "albumin", "alkaline", "SGOT",
            "platelets", "prothrombin", "histologic", "edema")
  pb <- pb[stats::complete.cases(pb[, vars]), , drop = FALSE]     # complete visits first,
  nv <- table(pb$id); pb <- pb[pb$id %in% names(nv)[nv >= pc$min_visits], , drop = FALSE]  # then >= 10
  pb <- pb[order(pb$id, pb$year), , drop = FALSE]
  X <- cbind(D_penicil = as.numeric(pb$drug == "D-penicil"), age = pb$age,
             female = as.numeric(pb$sex == "female"), log_bilirubin = log(pb$serBilir),
             albumin = pb$albumin, log_alkaline = log(pb$alkaline), log_SGOT = log(pb$SGOT),
             platelets = pb$platelets, prothrombin = pb$prothrombin, histologic = pb$histologic,
             edema_nodiur = as.numeric(pb$edema == "edema no diuretics"),
             edema_diur = as.numeric(pb$edema == "edema despite diuretics"), years = pb$year)
  y <- as.numeric(pb[[pc$outcome]] == "Yes")
  id <- as.integer(factor(pb$id))
  d <- make_data(X, y, id, pb$year / max(pb$year))
  d$y_binary <- y
  attr(d, "name") <- "PBC"
  attr(d, "description") <- sprintf("PBC: N=%d visits, n=%d patients (>= %d visits), p=%d, %s rate %.1f%%",
                                    length(y), length(unique(id)), pc$min_visits, ncol(X), pc$outcome, 100 * mean(y))
  d
}

## outcome transformation of the continuous data sets (CFG$real$transform)
transform_y <- function(d, name, cfg) {
  tr <- (cfg$real$transform %||% c())[name]
  if (!is.na(tr) && tr == "log") {
    d$y <- log(d$y)
    attr(d, "description") <- paste0(attr(d, "description"), ", outcome log-transformed")
  }
  d
}

load_real <- function(name, cfg) {
  d <- switch(name,
              BodyWeight = prep_BodyWeight(cfg),
              Soybean = prep_Soybean(cfg),
              Sprockett = prep_Sprockett(cfg),
              PBC = prep_PBC(cfg),
              OvarianDebulking = prep_OvarianDebulking(cfg),
              stop("unknown data set ", name))
  if (name %in% (cfg$real$datasets_cont %||% character(0))) d <- transform_y(d, name, cfg)
  d
}
