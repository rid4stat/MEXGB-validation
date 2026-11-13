cells <- sub("^b4_(.*)_v04f\\.rds$", "\\1", list.files("raw", full.names = FALSE, "^b4_.*_v04f\\.rds$"))
out <- NULL
for (c in cells) { n <- readRDS(sprintf("b4_%s_v04f.rds", c)); fam <- n$fam[1]
  o <- if (file.exists(sprintf("b4_%s_v031.rds", c))) readRDS(sprintf("b4_%s_v031.rds", c)) else if (file.exists(sprintf("b4_%s.rds", c))) readRDS(sprintf("b4_%s.rds", c)) else NULL
  m <- tapply(n$metric, n$method, mean); hib <- fam != "gaussian"
  old <- if (!is.null(o)) mean(o$metric[o$method == "MEXGB"]) else NA
  base <- if (fam == "cox") "XGBCox" else "XGBoost"
  w <- function(a, b) { x <- n$metric[n$method == a]; y <- n$metric[n$method == b]; sprintf("%d/%d", sum(if (hib) x > y else x < y), length(x)) }
  out <- rbind(out, data.frame(cell = c, MEXGB_031 = round(old, 3), MEXGB_04 = round(m["MEXGB"], 3), MEGB_04 = round(if ("MEGB" %in% names(m)) m["MEGB"] else NA, 3),
    XGB = round(m[base], 3), GBM = round(m["GBM"], 3), Oracle = round(m["Oracle"], 3), beats_XGB = w("MEXGB", base), beats_GBM = w("MEXGB", "GBM"),
    D11 = round(mean(n$D11[n$method == "MEXGB"]), 2), D22 = round(mean(n$D22[n$method == "MEXGB"]), 2), sec = round(mean(n$time[n$method == "MEXGB"])))) }
print(out, row.names = FALSE)
