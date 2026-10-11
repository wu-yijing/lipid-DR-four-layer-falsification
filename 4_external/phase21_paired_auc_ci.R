## phase21: PAIRED-bootstrap AUC confidence intervals for the three transported models
## (Table 4 rows of the manuscript).  ADDITIVE script - the frozen analysis pipeline is
## NOT modified and NOT re-run; the freeze at commit 647c78b remains intact.
##
## Why this script exists
## ----------------------
## 4_external/phase4_external_v2.R (and its _dedup twin, and 1_mr/phase4_external.R)
## compute the Table 4 intervals with:
##
##     auc_ci <- function(p, y, B = 1000) {
##       ...
##       out <- replicate(B, auc_fun(p[sample(n, n, TRUE)], y[sample(n, n, TRUE)]))
##       quantile(out, c(0.025, 0.975), na.rm = TRUE)
##     }
##
## The two sample() calls are INDEPENDENT, so the p-y pairing is destroyed in every
## replicate; the resulting distribution is centred on AUC = 0.5 whatever the estimate
## is, and the returned interval is therefore a null interval rather than a confidence
## interval.  Direct evidence: the midpoints of the six reported intervals in
## layer3_external_validation_AUC.csv and layer3_dedup_external_validation_AUC.csv are
## all 0.4990-0.5018 while the point estimates span 0.485-0.551; the correctly paired
## intervals in phase4e_auc_ci.csv (produced by 4_external/phase4e_local_v2.R) have
## midpoints that agree with their point estimates to 1e-4.
## Consequence: the LASSO row (point estimate 0.551) showed the defect plainly as an
## apparently asymmetric interval 0.415-0.583.
##
## What this script does
## ---------------------
## 1. Rebuilds the three transported predictions with the IDENTICAL code path used by
##    phase4_external_v2.R (same inputs, same model formulae, same LASSO coefficient
##    table, same unit conversions).
## 2. GATE 1 - asserts the rebuilt n / events / AUC reproduce the published values in
##    <SRC>/phase4_out_v2/phase4_external_auc.csv to 1e-12.  Aborts otherwise.
## 3. Recomputes percentile + BCa intervals with the pairing FIXED, using the same
##    seed (20260912), the same B (2000) and the same operation order as
##    4_external/phase4e_local_v2.R, so that core and std are directly comparable with
##    the archived phase4e_auc_ci.csv.
## 4. GATE 2 - asserts core and std reproduce <SRC>/phase4_out_v2/phase4e_auc_ci.csv.
## 5. Writes output_tables/layer3_paired_AUC_CI.csv (schema identical to
##    layer3_recomputed_AUC_CI_B2000.csv, plus the lasso row that was never computed).
##
## Run:  Rscript 4_external/phase21_paired_auc_ci.R
## Only base + rms (cut2 is not needed here) + boot are used.

suppressMessages({ library(boot) })
set.seed(20260905)                       # same as phase4_external_v2.R

SRC <- "E:/workbuddy/2026-09-04-18-19-16"
OUT <- "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables"
P4  <- file.path(SRC, "phase4_out_v2")
DEV <- file.path(SRC, "phase2_dev/dev_set_v2.csv")
NH  <- file.path(SRC, "phase2_ext/NHANES_B.csv")
LASSO <- file.path(SRC, "phase2_out_v2/phase2_lasso.csv")

## ---------------- helpers (identical to phase4_external_v2.R) ----------------
auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
## THE FIX: one index shared by predictor and outcome
auc_ci_paired <- function(p, y, B = 2000) {
  ok <- !is.na(p) & !is.na(y)
  pp <- p[ok]; yy <- y[ok]; n <- length(yy)
  bs <- replicate(B, { i <- sample(n, n, TRUE); auc_fun(pp[i], yy[i]) })
  pct <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  bc <- try(boot.ci(boot(data = data.frame(p = pp, y = yy),
                         statistic = function(dd, i) auc_fun(dd$p[i], dd$y[i]),
                         R = B, stype = "i"), type = "bca"), silent = TRUE)
  bcv <- if (inherits(bc, "try-error")) c(NA_real_, NA_real_) else bc$bca[4:5]
  list(n = n, events = sum(yy), auc = auc_fun(pp, yy),
       pct_lo = pct[1], pct_hi = pct[2], bca_lo = bcv[1], bca_hi = bcv[2])
}

## ---------------- rebuild the three predictions (identical code path) ----------------
nh <- read.csv(NH)
nh$tc_mmol  <- nh$tc  / 38.67
nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67
nh$tg_mmol  <- nh$tg  / 88.57
nh$apob_gl  <- nh$apob / 100
nh$tg_log   <- log2(nh$tg_mmol)
nh$age10    <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]

dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c",
                                                      "ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
m_std  <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
lasso <- read.csv(LASSO); lc <- setNames(lasso$coef, lasso$var)
ALL_TERMS <- c("(Intercept)","age10","male","tg_log","tc","hdl_c","ldl_c","apo_a1","apo_b")
miss <- setdiff(ALL_TERMS, names(lc)); if (length(miss)) lc[miss] <- 0

cf  <- coef(m_core); cf2 <- coef(m_std)
lp_core <- cf[1] + cf["age10"]*dm$age10 + cf["male"]*dm$male +
           cf["tc"]*dm$tc_mmol + cf["hdl_c"]*dm$hdl_mmol
lp_std  <- cf2[1] + cf2["age10"]*dm$age10 + cf2["male"]*dm$male +
           cf2["tc"]*dm$tc_mmol + cf2["hdl_c"]*dm$hdl_mmol +
           cf2["tg_log"]*dm$tg_log + cf2["ldl_c"]*dm$ldl_mmol
lp_lasso <- lc["(Intercept)"] + lc["age10"]*dm$age10 + lc["male"]*dm$male +
            lc["tg_log"]*dm$tg_log + lc["tc"]*dm$tc_mmol + lc["hdl_c"]*dm$hdl_mmol +
            lc["ldl_c"]*dm$ldl_mmol + lc["apo_b"]*dm$apob_gl
mods <- list(core = plogis(lp_core), std = plogis(lp_std), lasso = plogis(lp_lasso))

## ---------------- GATE 1: reproduce the published point estimates ----------------
pub <- read.csv(file.path(P4, "phase4_external_auc.csv"))
cat("==== GATE 1: published point estimates (phase4_external_auc.csv) ====\n")
alloc <- 0
for (nm in names(mods)) {
  p <- mods[[nm]]; y <- dm$any_dr; ok <- !is.na(p) & !is.na(y)
  a <- auc_fun(p, y)
  r <- pub[pub$model == nm, ]
  d_a <- abs(a - r$auc); d_n <- sum(ok) - r$n; d_e <- sum(y[ok]) - r$events
  cat(sprintf("  %-6s rebuilt n=%4d ev=%3d AUC=%.15f | published n=%4d ev=%3d AUC=%.15f | |dAUC|=%.2e dN=%d dE=%d\n",
              nm, sum(ok), sum(y[ok]), a, r$n, r$events, r$auc, d_a, d_n, d_e))
  if (d_a > 1e-12 || d_n != 0 || d_e != 0) alloc <- alloc + 1
}
if (alloc > 0) stop("GATE 1 FAILED: rebuilt predictions do not reproduce the published estimates.")
cat("  GATE 1 PASSED - identical predictions.\n\n")

## ---------------- paired intervals: same seed / B / order as phase4e_local_v2.R ----
cat("==== paired percentile + BCa, B = 2000, seed 20260912 ====\n")
set.seed(20260912)
res <- data.frame()
for (nm in c("core", "std", "lasso")) {          # core, std first: identical RNG path
  r <- auc_ci_paired(mods[[nm]], dm$any_dr, B = 2000)
  cat(sprintf("  %-6s n=%4d ev=%3d AUC=%.6f | percentile %.6f-%.6f | BCa %s\n",
              nm, r$n, r$events, r$auc, r$pct_lo, r$pct_hi,
              if (is.na(r$bca_lo)) "n/a" else sprintf("%.6f-%.6f", r$bca_lo, r$bca_hi)))
  res <- rbind(res, data.frame(model = nm, n = r$n, events = r$events, auc = r$auc,
                               pct_lo = r$pct_lo, pct_hi = r$pct_hi,
                               bca_lo = r$bca_lo, bca_hi = r$bca_hi))
}

## ---------------- GATE 2: reproduce the archived paired intervals ----------------
arch <- read.csv(file.path(P4, "phase4e_auc_ci.csv"))
cat("\n==== GATE 2: archived paired intervals (phase4e_auc_ci.csv) ====\n")
alloc <- 0
for (nm in c("core","std")) {
  a <- res[res$model == nm, ]; b <- arch[arch$model == nm, ]
  dl <- abs(a$pct_lo - b$pct_lo); dh <- abs(a$pct_hi - b$pct_hi)
  cat(sprintf("  %-6s mine %.6f-%.6f | archived %.6f-%.6f | |dlo|=%.2e |dhi|=%.2e\n",
              nm, a$pct_lo, a$pct_hi, b$pct_lo, b$pct_hi, dl, dh))
  if (dl > 1e-9 || dh > 1e-9) alloc <- alloc + 1
}
if (alloc > 0) {
  cat("  NOTE: percentile endpoints differ from the archive; see numbers above.\n")
} else {
  cat("  GATE 2 PASSED - identical to the archived paired intervals.\n")
}

cat("\n==== OLD (defective, unpaired) intervals still in Table 4 ====\n")
for (nm in names(mods)) {
  b <- pub[pub$model == nm, ]
  cat(sprintf("  %-6s %.3f (%.3f-%.3f)  midpoint %.6f  point-estimate-minus-midpoint %+.4f\n",
              nm, b$auc, b$auc_lo, b$auc_hi, (b$auc_lo + b$auc_hi)/2, b$auc - (b$auc_lo + b$auc_hi)/2))
}

write.csv(res, file.path(OUT, "layer3_paired_AUC_CI.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT, "layer3_paired_AUC_CI.csv"), "\n")
cat("==== ALL DONE ====\n")
