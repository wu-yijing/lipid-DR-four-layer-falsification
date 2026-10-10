## =====================================================================
## phase11 (layer 4): incremental value of lipids ON TOP OF the full
## reference model
##
## Why this exists
## ---------------
## The manuscript reports the increment in the opposite direction only: full
## reference model minus a lipid-only model (dAUC 0.218). That answers "how much do
## the clinical variables add over lipids", which is not the question a bedside tool
## designer asks. The question that matters clinically is the reverse one, and it is
## the one the paper's exclusionary conclusion actually needs:
##
##     given the full reference model, do lipids add anything?
##
## If the answer is no, the conclusion upgrades from "lipids are not sufficient on
## their own" to "lipids add nothing on top of what is already measured".
##
## Conventions are copied verbatim from phase10_local_comparator_increment.R and
## phase4_external.R so the two directions are directly comparable:
##   * out-of-fold predictions, 10-fold, the same folds as the committed layer-4 run
##   * IDI = (mean p_new - mean p_old | cases) - (mean p_new - mean p_old | non-cases)
##   * NRI = categorical, cut = median of the NEW model's predicted risk at the point
##     estimate and held fixed across replicates
##   * intervals = frozen-model percentile bootstrap (participants resampled, fitted
##     predictions reused), B = 1000
##   * both unweighted and MEC-weighted (WTMEC2YR/2)
##
## Subsets. Core lipids (TC, HDL-C) are non-fasting and available throughout; the
## extended set adds TG and LDL-C, which exist only in the fasting subsample. Each
## subset refits BOTH models on its own rows so that every increment is paired.
##
## Outputs
## -------
## output_tables/layer4_lipid_increment_over_full_model.csv
##
## Run
## ---
##   LDF_SRC=/path/to/analysis/dir Rscript 4_external/phase11_lipid_increment_over_full_model.R
## =====================================================================
suppressMessages({ library(survey) })
set.seed(20260912)

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- if (is.na(this)) "." else normalizePath(file.path(dirname(this), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "output_tables"))) REPO <- "."
OUT <- file.path(REPO, "output_tables")
B <- as.integer(Sys.getenv("LDF_B", "1000"))

NH <- file.path(SRC, "phase2_ext/NHANES_B.csv")
DEV <- file.path(SRC, "phase2_dev/dev_set_v2.csv")
if (!file.exists(NH)) stop("NHANES_B.csv not found; set LDF_SRC", call. = FALSE)
if (!file.exists(DEV)) stop("dev_set_v2.csv not found; set LDF_SRC", call. = FALSE)
cat(sprintf("LDF_SRC = %s\nwrites to = %s\nB = %d\n", SRC, OUT, B))

## ---------------- helpers (identical to phase10) ----------------
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
brier_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); mean((p[ok] - y[ok])^2)
}
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
oof_fixed <- function(fmla, dat, fold, wcol = NULL) {
  p <- rep(NA_real_, nrow(dat))
  for (k in sort(unique(fold))) {
    tr <- dat[fold != k, ]; te <- dat[fold == k, ]
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
    stopifnot(length(pr) == sum(okk))
    p[which(fold == k)[okk]] <- pr
  }
  p
}
boot_metric <- function(p0, p1, yy, ww, f, rep_n) {
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

## ---------------- load ----------------
nh <- read.csv(NH)
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
nh$nonhdl_mmol <- nh$tc_mmol - nh$hdl_mmol
dm <- nh[nh$diabetes == 1, ]

## the Chinese transportable core coefficients, exactly as in phase4e/phase10
dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(
  dev[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
cf <- coef(glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial()))

ref_vars <- c("age10", "male", "hba1c_pct", "dm_duration", "egfr_ckdepi2021", "bmi",
              "smoker_100", "insulin_use", "oral_med", "sys_bp", "any_dr", "WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR / 2
dref$p_transported <- plogis(cf[1] + cf["age10"] * dref$age10 + cf["male"] * dref$male +
                               cf["tc"] * dref$tc_mmol + cf["hdl_c"] * dref$hdl_mmol)

fold <- sample(rep(1:10, length.out = nrow(dref)))

cat(sprintf("\nn(diabetes) = %d | n(full reference model, complete) = %d | events = %d\n",
            nrow(dm), nrow(dref), sum(dref$any_dr)))

## ---------------- models ----------------
fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_plus_core <- update(fml_full, . ~ . + tc_mmol + hdl_mmol)
fml_plus_std  <- update(fml_full, . ~ . + tc_mmol + hdl_mmol + tg_log + ldl_mmol)
fml_plus_apob <- update(fml_plus_std, . ~ . + apob)
fml_plus_transported <- update(fml_full, . ~ . + p_transported)

p_full <- oof_fixed(fml_full, dref, fold, "wt4")
## anchor check: must reproduce the committed layer-4 full model exactly
stopifnot(abs(auc_fun(p_full, dref$any_dr) - 0.7330631) < 5e-4)
cat(sprintf("anchor OK: full model OOF AUC = %.4f (weighted %.4f)\n",
            auc_fun(p_full, dref$any_dr), auc_w_fun(p_full, dref$any_dr, dref$wt4)))

## ---------------- run one subset ----------------
run_subset <- function(label, extra_vars, fml_new, note) {
  key <- c("age10", "male", "hba1c_pct", "dm_duration", "egfr_ckdepi2021", "bmi",
           "smoker_100", "insulin_use", "oral_med", "sys_bp", "any_dr", "wt4", extra_vars)
  keep <- complete.cases(dref[, key])
  sub <- dref[keep, ]
  fl <- fold[keep]
  p0 <- p_full[keep]
  p1 <- oof_fixed(fml_new, sub, fl, "wt4")
  ok <- !is.na(p0) & !is.na(p1)
  sub <- sub[ok, ]; p0 <- p0[ok]; p1 <- p1[ok]
  yy <- sub$any_dr; ww <- sub$wt4
  cat(sprintf("\n---- %s: n = %d, events = %d (%.1f%%) ----\n",
              label, length(yy), sum(yy), 100 * mean(yy)))
  out <- list()
  for (wt in c("unweighted", "MEC-weighted")) {
    wv <- if (wt == "unweighted") NULL else ww
    cutv <- median(p1, na.rm = TRUE)
    f <- function(a, b, yv, wvv) c(AUC_old = if (is.null(wvv)) auc_fun(a, yv) else auc_w_fun(a, yv, wvv),
                                   AUC_new = if (is.null(wvv)) auc_fun(b, yv) else auc_w_fun(b, yv, wvv),
                                   dAUC = (if (is.null(wvv)) auc_fun(b, yv) else auc_w_fun(b, yv, wvv)) -
                                          (if (is.null(wvv)) auc_fun(a, yv) else auc_w_fun(a, yv, wvv)),
                                   nri_idi(a, b, yv, cutv, wvv))
    pt <- f(p0, p1, yy, wv)
    bs <- boot_metric(p0, p1, yy, wv, f, B)
    ci <- apply(bs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    cat(sprintf("   %-12s AUC full %.4f -> full+lipid %.4f | dAUC %.4f (%.4f-%.4f) | NRI %.4f (%.4f-%.4f) | IDI %.4f (%.4f-%.4f) | Brier %.4f -> %.4f\n",
                wt, pt[["AUC_old"]], pt[["AUC_new"]], pt[["dAUC"]], ci[1, "dAUC"], ci[2, "dAUC"],
                pt[["NRI"]], ci[1, "NRI"], ci[2, "NRI"],
                pt[["IDI"]], ci[1, "IDI"], ci[2, "IDI"],
                brier_fun(p0, yy), brier_fun(p1, yy)))
    out[[length(out) + 1]] <- data.frame(
      contrast = label, added = paste(extra_vars, collapse = " + "), weighting = wt,
      n = length(yy), events = sum(yy),
      AUC_full = pt[["AUC_old"]], AUC_full_lo = ci[1, "AUC_old"], AUC_full_hi = ci[2, "AUC_old"],
      AUC_plus = pt[["AUC_new"]], AUC_plus_lo = ci[1, "AUC_new"], AUC_plus_hi = ci[2, "AUC_new"],
      dAUC = pt[["dAUC"]], dAUC_lo = ci[1, "dAUC"], dAUC_hi = ci[2, "dAUC"],
      NRI = pt[["NRI"]], NRI_lo = ci[1, "NRI"], NRI_hi = ci[2, "NRI"],
      IDI = pt[["IDI"]], IDI_lo = ci[1, "IDI"], IDI_hi = ci[2, "IDI"],
      Brier_full = brier_fun(p0, yy), Brier_plus = brier_fun(p1, yy),
      cut = cutv, basis = sprintf("frozen-model percentile bootstrap, B=%d", B), note = note)
  }
  out
}

rows <- c(
  run_subset("full model + core lipids (TC, HDL-C)", c("tc_mmol", "hdl_mmol"), fml_plus_core,
             "non-fasting lipids available throughout; same n as the committed full model"),
  run_subset("full model + extended lipids (TC, HDL-C, TG, LDL-C)",
             c("tc_mmol", "hdl_mmol", "tg_log", "ldl_mmol"), fml_plus_std,
             "evaluable on the fasting subsample only"),
  run_subset("full model + extended lipids + ApoB",
             c("tc_mmol", "hdl_mmol", "tg_log", "ldl_mmol", "apob"), fml_plus_apob,
             "exploratory; ApoB exists only in the fasting subsample"),
  run_subset("full model + TRANSPORTED Chinese lipid score (single term)",
             "p_transported", fml_plus_transported,
             "the transported score carries the same information as the lipids it is built from, so a model with the raw lipids should dominate it")
)

res <- do.call(rbind, rows)

## ---------------- negative control: the same number of useless predictors -------
## Adding any k predictors costs degrees of freedom, and out-of-fold AUC falls with
## overfitting when the event count is fixed. So a negative dAUC is not by itself
## evidence that lipids carry harmful or merely uninformative signal. The control
## replaces the lipids with k independent standard-normal columns and repeats the
## whole out-of-fold fit: if the lipid increment is indistinguishable from this
## distribution, the defensible reading is "no gain", not "harm".
run_noise_control <- function(label, k, extra_vars, reps = 20) {
  key <- c("age10", "male", "hba1c_pct", "dm_duration", "egfr_ckdepi2021", "bmi",
           "smoker_100", "insulin_use", "oral_med", "sys_bp", "any_dr", "wt4", extra_vars)
  keep <- complete.cases(dref[, key])
  sub <- dref[keep, ]; fl <- fold[keep]; p0 <- p_full[keep]
  yy <- sub$any_dr
  d0 <- auc_fun(p0, yy)
  out <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(20260912 + r)
    z <- as.data.frame(matrix(rnorm(nrow(sub) * k), nrow(sub), k))
    names(z) <- paste0("noise", seq_len(k))
    zz <- cbind(sub, z)
    fml <- as.formula(paste("any_dr ~", paste(c(all.vars(fml_full)[-1], names(z)), collapse = " + ")))
    pn <- oof_fixed(fml, zz, fl, "wt4")
    out[r] <- auc_fun(pn, yy) - d0
  }
  cat(sprintf("\n---- negative control: %s, %d noise columns, %d draws, n = %d ----\n",
              label, k, reps, nrow(sub)))
  cat(sprintf("   dAUC across draws: median %.4f, range %.4f to %.4f, 2.5th-97.5th %.4f to %.4f\n",
              median(out), min(out), max(out),
              quantile(out, .025), quantile(out, .975)))
  data.frame(contrast = paste0("NEGATIVE CONTROL: full model + ", k, " noise columns (", label, ")"),
             added = paste0(k, " standard-normal draw(s)"), weighting = "unweighted",
             n = nrow(sub), events = sum(yy),
             AUC_full = d0, AUC_full_lo = NA, AUC_full_hi = NA,
             AUC_plus = NA, AUC_plus_lo = NA, AUC_plus_hi = NA,
             dAUC = min(out), dAUC_lo = quantile(out, .025), dAUC_hi = max(out),
             NRI = NA, NRI_lo = NA, NRI_hi = NA, IDI = NA, IDI_lo = NA, IDI_hi = NA,
             Brier_full = NA, Brier_plus = NA, cut = NA,
             basis = sprintf("out-of-fold refit, %d draws; dAUC range across draws", reps),
             note = "reported in dAUC/dAUC_lo/dAUC_hi as min / 2.5th percentile / max of the draw distribution")
}
ctrl <- rbind(
  run_noise_control("core subset", 2L, c("tc_mmol", "hdl_mmol")),
  run_noise_control("extended subset", 4L, c("tc_mmol", "hdl_mmol", "tg_log", "ldl_mmol")),
  run_noise_control("ApoB subset", 5L, c("tc_mmol", "hdl_mmol", "tg_log", "ldl_mmol", "apob"))
)
res_all <- rbind(res, ctrl)
num <- sapply(res_all, is.numeric)
res_all[num] <- lapply(res_all[num], function(v) ifelse(is.na(v), NA, round(v, 5)))
f <- file.path(OUT, "layer4_lipid_increment_over_full_model.csv")
write.csv(res_all, f, row.names = FALSE)
write.csv(ctrl, file.path(OUT, "layer4_lipid_increment_negative_control.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\nwrote %s\n==== ALL DONE ====\n", f,
            file.path(OUT, "layer4_lipid_increment_negative_control.csv")))
