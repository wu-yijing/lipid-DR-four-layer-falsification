## =====================================================================
## phase10 (layer 4): increment metrics for the LOCALLY FITTED comparator
##
## Two jobs
## --------
## (1) Fill a reporting gap. The manuscript reports dAUC, NRI and IDI for the
##     comparison between the full reference model and the TRANSPORTED Chinese lipid
##     score - the one comparison the paper itself argues cannot be attributed to
##     either cause, because it mixes transport loss with information gain. The
##     comparison that CAN be attributed (full model against a lipid-only model
##     fitted locally in NHANES) had dAUC only. This script supplies NRI and IDI for
##     it, with 95% intervals for every increment metric.
##
## (2) Repair a silent misalignment bug inherited from phase4e_local_v2.R.
##     That script built out-of-fold predictions with
##
##         p[fold == k] <- plogis(model.matrix(fmla, data = te) %*% coef(f))
##
##     For any model whose predictors contain missing values this is wrong:
##     model.matrix() drops the incomplete rows, so the replacement vector is
##     shorter than p[fold == k] and R recycles it with a warning. Only the
##     extended lipid model is affected here (its terms tg_log and ldl_mmol are
##     missing in 473 of the 888 reference records), and the consequence is that the
##     committed row `local_lipid_only_std` in layer4_three_row_AUC.csv reports
##     n = 888 and AUC 0.4861 for a quantity that exists on 415 records and is
##     genuinely 0.4988 (95% CI 0.437-0.558).
##
##     oof_fixed() below writes each fold's predictions back to their own row
##     positions and asserts that the prediction vector has exactly the length of
##     the aligned index, so the failure mode cannot recur silently. The same fix is
##     applied to phase4e_local_v2.R and phase4e_local_v2_dedup.R; the committed
##     three-row file is left untouched and is flagged as defective in
##     output_tables/PROVENANCE.md, with the corrected values written here.
##
## The full-vs-transported row doubles as a reproducibility check on the values
## already committed, and the script aborts if they no longer reproduce.
##
## Metric conventions follow 1_mr/phase4_external.R:
##   IDI = (mean p_new - mean p_old | cases) - (mean p_new - mean p_old | non-cases)
##   NRI = categorical, at cut = median of the new model's predicted risk
##
## Outputs
## -------
## output_tables/layer4_local_comparator_increment.csv
## output_tables/layer4_three_row_AUC_corrected.csv
##
## Run
## ---
##   Rscript 4_external/phase10_local_comparator_increment.R
##   LDF_SRC=/path/to/analysis/dir LDF_B=2000 Rscript 4_external/phase10_local_comparator_increment.R
## =====================================================================
suppressMessages({ library(survey) })
set.seed(20260912)

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- if (is.na(this)) "." else normalizePath(file.path(dirname(this), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "output_tables"))) REPO <- "."
OUT <- file.path(REPO, "output_tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DEV <- file.path(SRC, "phase2_dev/dev_set_v2.csv")
NH  <- file.path(SRC, "phase2_ext/NHANES_B.csv")
B   <- as.integer(Sys.getenv("LDF_B", "1000"))

for (p in c(DEV, NH)) if (!file.exists(p)) {
  stop(sprintf(paste0("required input not found: %s\n",
    "These inputs are not distributed with the repository.\n",
    "Set LDF_SRC to a local analysis directory, e.g.\n",
    "    export LDF_SRC=/path/to/analysis/dir\n"), p), call. = FALSE)
}
cat(sprintf("LDF_SRC = %s\nwrites to = %s\nbootstrap replicates B = %d\n", SRC, OUT, B))

## ---------------- helpers ----------------
auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_w_fun <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  p <- p[ok]; y <- y[ok]; w <- w[ok]
  pos <- y == 1; neg <- y == 0
  if (!any(pos) || !any(neg)) return(NA_real_)
  M <- outer(p[pos], p[neg], function(a, b) (a > b) + 0.5 * (a == b))
  sum(w[pos] * as.vector(M %*% w[neg])) / (sum(w[pos]) * sum(w[neg]))
}
## w = NULL reproduces the unweighted versions exactly
nri_idi <- function(p0, p1, y, cut, w = NULL) {
  ok <- !is.na(p0) & !is.na(p1) & !is.na(y)
  p0 <- p0[ok]; p1 <- p1[ok]; y <- y[ok]
  if (is.null(w)) w <- rep(1, length(y)) else w <- w[ok]
  ew <- function(v, sel) sum(w[sel] * v[sel]) / sum(w[sel])
  ev <- y == 1
  idi <- (ew(p1, ev) - ew(p0, ev)) - (ew(p1, !ev) - ew(p0, !ev))
  up_e   <- ew(as.numeric(p1 >= cut & p0 <  cut), ev)
  down_e <- ew(as.numeric(p1 <  cut & p0 >= cut), ev)
  up_n   <- ew(as.numeric(p1 >= cut & p0 <  cut), !ev)
  down_n <- ew(as.numeric(p1 <  cut & p0 >= cut), !ev)
  c(NRI = up_e - down_e + down_n - up_n, IDI = idi)
}

## ---------------- load + convert (identical to phase4e_local_v2.R) ----------------
nh <- read.csv(NH)
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]

dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
cf <- coef(glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial()))

ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR / 2

## ---------------- common 10-fold split (identical to phase4e_local_v2.R) ----------------
fold <- sample(rep(1:10, length.out = nrow(dref)))

## Alignment-preserving out-of-fold prediction.
## The guard is the point: sum(okk) must equal length(pred), otherwise a model
## whose predictors carry missing values would silently scramble the vector.
oof_fixed <- function(fmla, dat, wcol = NULL) {
  p <- rep(NA_real_, nrow(dat))
  for (k in 1:10) {
    tr <- dat[fold != k, ]; te <- dat[fold == k, ]
    ## survey weights are sampling weights, not binomial counts, so svyglm reports
    ## "non-integer #successes"; expected and harmless
    f <- suppressWarnings(try(
      if (is.null(wcol)) glm(fmla, data = tr, family = binomial())
      else survey::svyglm(fmla, design = survey::svydesign(ids = ~1,
             weights = as.formula(paste0("~", wcol)), data = tr), family = binomial()),
      silent = TRUE))
    if (inherits(f, "try-error")) next
    okk <- complete.cases(model.frame(fmla, data = te, na.action = na.pass))
    if (!any(okk)) next
    pr <- suppressWarnings(as.numeric(predict(f, newdata = te[okk, , drop = FALSE],
                                              type = "response")))
    stopifnot(length(pr) == sum(okk))   # <- alignment guard; must be able to fail
    p[which(fold == k)[okk]] <- pr
  }
  p
}
fml_full     <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
                smoker_100 + insulin_use + oral_med + sys_bp
fml_lip_core <- any_dr ~ age10 + male + tc_mmol + hdl_mmol
fml_lip_std  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol + tg_log + ldl_mmol

p_full     <- oof_fixed(fml_full, dref, "wt4")
p_loc_core <- oof_fixed(fml_lip_core, dref, "wt4")
p_loc_std  <- oof_fixed(fml_lip_std, dref, "wt4")
p_tr_core  <- plogis(cf[1] + cf["age10"] * dref$age10 + cf["male"] * dref$male +
                       cf["tc"] * dref$tc_mmol + cf["hdl_c"] * dref$hdl_mmol)

## ---------------- subset A: the three-row comparison ----------------
ok <- !is.na(p_full) & !is.na(p_loc_core) & !is.na(p_tr_core)
ds <- dref[ok, ]; y <- ds$any_dr; w <- ds$wt4
P <- list(full = p_full[ok], local_core = p_loc_core[ok], transported = p_tr_core[ok])
cat(sprintf("\nsubset A (full + local core + transported): n=%d events=%d (%.1f%%)\n",
            length(y), sum(y), 100 * mean(y)))
for (nm in names(P)) cat(sprintf("  AUC %-12s unweighted %.4f  MEC-weighted %.4f\n",
                                 nm, auc_fun(P[[nm]], y), auc_w_fun(P[[nm]], y, w)))

stopifnot(abs(auc_fun(P$full, y) - 0.7330631) < 5e-4,
          abs(auc_fun(P$transported, y) - 0.4859867) < 5e-4,
          abs(auc_fun(P$local_core, y) - 0.5150265) < 5e-4)
cat("reproducibility check passed against layer4_three_row_AUC.csv (full, local core, transported)\n")

## ---------------- subset B: the extended lipid model ----------------
ok_std <- !is.na(p_full) & !is.na(p_loc_std)
ys <- dref$any_dr[ok_std]; ws <- dref$wt4[ok_std]
ps_full <- p_full[ok_std]; ps_std <- p_loc_std[ok_std]
cat(sprintf("\nsubset B (full + extended lipid model): n=%d events=%d (%.1f%%)\n",
            sum(ok_std), sum(ys), 100 * mean(ys)))

a_std_corrected <- auc_fun(ps_std, ys)
committed <- NULL
cc_path <- file.path(OUT, "layer4_three_row_AUC.csv")
if (file.exists(cc_path)) {
  cc <- read.csv(cc_path)
  rowcc <- cc[cc$row == "local_lipid_only_std", ]
  if (nrow(rowcc) == 1) committed <- rowcc
}
if (!is.null(committed)) {
  cat("\n========== DEFECT REPORT (inherited from phase4e_local_v2.R) ==========\n")
  cat(sprintf("  committed layer4_three_row_AUC.csv, local_lipid_only_std : n=%d AUC=%.4f\n",
              committed$n, committed$AUC_unweighted))
  cat(sprintf("  corrected value here                                     : n=%d AUC=%.4f\n",
              sum(ok_std), a_std_corrected))
  cat("  cause: model.matrix() dropped the 473 records with missing TG/LDL-C, so the\n")
  cat("  fold-wise replacement vector was shorter than p[fold == k] and was recycled.\n")
  cat("  effect on the manuscript sentence about the extended variable set: the n was\n")
  cat("  overstated 2.1x and the AUC was 0.4861 instead of 0.4988. The conclusion\n")
  cat("  (at or below chance) is unchanged.\n")
  cat("========================================================================\n")
}
## the guard has to be able to fail: the corrected n must NOT be 888
stopifnot(sum(ok_std) < nrow(dref), abs(a_std_corrected - 0.4988) < 5e-3)

## ---------------- increments with bootstrap intervals ----------------
## Intervals are conditional on the fitted models: participants are resampled and the
## frozen predictions reused. The basis is written into the output so an interval
## cannot be mistaken for one that re-estimates the model.
## The bootstrap carries the two AUCs alongside the increments so that the level of
## discrimination being differenced is itself reported with an interval. Without it the
## corrected extended-variable-set row would state a point estimate (0.4988) that no
## committed table brackets, which is the reporting gap the audit flags elsewhere.
boot_metric <- function(p0, p1, yy, ww, f, rep_n = B) {
  nm <- c("AUC_old", "AUC_new", "dAUC", "NRI", "IDI")
  out <- matrix(NA_real_, rep_n, length(nm), dimnames = list(NULL, nm))
  n <- length(yy)
  for (i in seq_len(rep_n)) {
    s <- sample(n, n, TRUE)
    out[i, ] <- tryCatch(f(p0[s], p1[s], yy[s], ww[s]),
                         error = function(e) setNames(rep(NA_real_, length(nm)), nm))
  }
  out
}

run_pair <- function(p0, p1, yy, ww, label, old_name, note, wt) {
  cutv <- median(p1, na.rm = TRUE)   # fixed at the point estimate so the NRI
                                     # definition does not drift across replicates
  f <- function(a, b, yv, wv) c(AUC_old = auc_fun(a, yv),
                                AUC_new = auc_fun(b, yv),
                                dAUC = auc_fun(b, yv) - auc_fun(a, yv),
                                nri_idi(a, b, yv, cutv, wv))
  pt <- f(p0, p1, yy, ww)
  bs <- boot_metric(p0, p1, yy, ww, f)
  ci <- apply(bs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(contrast = label, comparator = old_name, weighting = wt,
             n = length(yy), events = sum(yy),
             AUC_old = pt[["AUC_old"]], AUC_old_lo = ci[1, "AUC_old"], AUC_old_hi = ci[2, "AUC_old"],
             AUC_new = pt[["AUC_new"]], AUC_new_lo = ci[1, "AUC_new"], AUC_new_hi = ci[2, "AUC_new"],
             dAUC = pt[["dAUC"]], dAUC_lo = ci[1, "dAUC"], dAUC_hi = ci[2, "dAUC"],
             NRI = pt[["NRI"]], NRI_lo = ci[1, "NRI"], NRI_hi = ci[2, "NRI"],
             IDI = pt[["IDI"]], IDI_lo = ci[1, "IDI"], IDI_hi = ci[2, "IDI"],
             cut = cutv, basis = sprintf("frozen-model percentile bootstrap, B=%d", B),
             note = note)
}

rows <- list()
for (wt in c("unweighted", "MEC-weighted")) {
  ww <- if (wt == "unweighted") NULL else w
  rows[[length(rows) + 1]] <- run_pair(
    P$local_core, P$full, y, ww,
    "full model minus locally fitted lipid-only model (CORE)",
    "local_core",
    "the attributable contrast: isolates the gain from the non-lipid clinical information",
    wt)
  rows[[length(rows) + 1]] <- run_pair(
    P$transported, P$full, y, ww,
    "full model minus TRANSPORTED Chinese lipid score",
    "transported",
    "the previously reported contrast; mixes transport loss with information gain, so it cannot be attributed to either",
    wt)
  rows[[length(rows) + 1]] <- run_pair(
    P$transported, P$local_core, y, ww,
    "locally fitted lipid-only model (CORE) minus TRANSPORTED score",
    "transported",
    "isolates the transport component alone; indistinguishable from zero here",
    wt)
}
for (wt in c("unweighted", "MEC-weighted")) {
  wws <- if (wt == "unweighted") NULL else ws
  rows[[length(rows) + 1]] <- run_pair(
    ps_std, ps_full, ys, wws,
    "full model minus locally fitted lipid-only model (EXTENDED variable set)",
    "local_std",
    sprintf("evaluable on %d of %d records only (TG and LDL-C missing elsewhere); corrected for the phase4e alignment defect",
            sum(ok_std), nrow(dref)),
    wt)
}

res <- do.call(rbind, rows)
num <- sapply(res, is.numeric)
res[num] <- lapply(res[num], function(v) round(v, 5))
write.csv(res, file.path(OUT, "layer4_local_comparator_increment.csv"), row.names = FALSE)
print(res[, c("contrast", "weighting", "n", "dAUC", "NRI", "IDI")], row.names = FALSE)

## ---------------- corrected three-row table ----------------
corr <- data.frame(
  row = c("transported_Chinese_lipid", "local_lipid_only_core",
          "local_lipid_only_std", "full_reference_model"),
  n = c(length(y), length(y), sum(ok_std), length(y)),
  events = c(sum(y), sum(y), sum(ys), sum(y)),
  AUC_unweighted = c(auc_fun(P$transported, y), auc_fun(P$local_core, y),
                     a_std_corrected, auc_fun(P$full, y)),
  AUC_MEC_weighted = c(auc_w_fun(P$transported, y, w), auc_w_fun(P$local_core, y, w),
                       auc_w_fun(ps_std, ys, ws), auc_w_fun(P$full, y, w)),
  status = c("unchanged", "unchanged",
             "CORRECTED (committed row had n=888 and AUC 0.4861/0.4972)", "unchanged"))
corr[, 3:5] <- lapply(corr[, 3:5], function(v) round(v, 6))
write.csv(corr, file.path(OUT, "layer4_three_row_AUC_corrected.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\nwrote %s\n",
            file.path(OUT, "layer4_local_comparator_increment.csv"),
            file.path(OUT, "layer4_three_row_AUC_corrected.csv")))
cat("==== ALL DONE ====\n")
