## AUTO-GENERATED v2 from phase7_out/phase7h_weighted_pms.R (corrected outcome)

# Phase 7h: 加权外验 + Brier 全套 + pmsampsize（R2-2 / R2-3 / R2-4 / R2-7 相关联项）
# 1) 三个转移版本(core/std/LASSO)的 NHANES survey 加权 AUC/校准/Brier（WTMEC2YR/2, ids=~1）
# 2) 完整参考模型（加权 svyglm）10折 OOF + 血脂-only 的加权对比（新 fold 种子，已在日志注明）
# 3) pmsampsize（Riley）：中国血脂模型 与 NHANES 完整模型
suppressMessages({ library(survey); library(pmsampsize); library(rms) })
options(survey.lonely.psu = "adjust")
P2ext <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext"
P2dev <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
P2out <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
OUT   <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out"

## ---------------- helpers ----------------
aucw <- function(p, y, w) {   # weighted Mann-Whitney AUC
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  p <- p[ok]; y <- y[ok]; w <- w[ok]
  i1 <- which(y == 1); i0 <- which(y == 0)
  W1 <- sum(w[i1]); W0 <- sum(w[i0])
  if (W1 == 0 || W0 == 0) return(NA_real_)
  num <- 0
  for (i in i1) {
    num <- num + w[i] * sum(w[i0] * (p[i] > p[i0]) + 0.5 * w[i0] * (p[i] == p[i0]))
  }
  num / (W1 * W0)
}
brierw <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  sum(w[ok] * (p[ok] - y[ok])^2) / sum(w[ok])
}
calibw <- function(p, y, w) {  # survey-weighted recalibration
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  lp <- qlogis(pmin(pmax(p[ok], 1e-6), 1 - 1e-6))
  dd <- data.frame(y = y[ok], lp = lp, w = w[ok])
  f <- svyglm(y ~ lp, design = svydesign(ids = ~1, weights = ~w, data = dd),
              family = binomial())
  c(intercept = as.numeric(coef(f)[1]), slope = as.numeric(coef(f)[2]))
}
auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
brier_fun <- function(p, y) { ok <- !is.na(p) & !is.na(y); mean((p[ok] - y[ok])^2) }

## ---------------- data & transferred models (replicate phase4) ----------------
nh <- read.csv(file.path(P2ext, "NHANES_B.csv"))
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$apob_gl <- nh$apob / 100
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]
dm$wt4 <- dm$WTMEC2YR / 2
cat(sprintf("NHANES diabetes n=%d events=%d (%.1f%%) | wt4 range [%.0f, %.0f]\n",
            nrow(dm), sum(dm$any_dr), 100 * mean(dm$any_dr),
            min(dm$wt4, na.rm = TRUE), max(dm$wt4, na.rm = TRUE)))
cat(sprintf("weight missing: %d\n", sum(is.na(dm$wt4))))

dev <- read.csv(file.path(P2dev, "dev_set_v2.csv"))
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
m_std  <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
lasso <- read.csv(file.path("E:/workbuddy/2026-09-04-18-19-16/phase2_out_v2", "phase2_lasso.csv")); lc <- setNames(lasso$coef, lasso$var)
ALL_TERMS <- c("(Intercept)","age10","male","tg_log","tc","hdl_c","ldl_c","apo_a1","apo_b")
miss <- setdiff(ALL_TERMS, names(lc))
if (length(miss)) { cat("LASSO dropped:", paste(miss, collapse=", "), "-> set to 0\n"); lc[miss] <- 0 }

cf <- coef(m_core)
lp_core <- cf[1] + cf["age10"]*dm$age10 + cf["male"]*dm$male + cf["tc"]*dm$tc_mmol + cf["hdl_c"]*dm$hdl_mmol
cf2 <- coef(m_std)
lp_std <- cf2[1] + cf2["age10"]*dm$age10 + cf2["male"]*dm$male + cf2["tc"]*dm$tc_mmol +
          cf2["hdl_c"]*dm$hdl_mmol + cf2["tg_log"]*dm$tg_log + cf2["ldl_c"]*dm$ldl_mmol
lp_lasso <- lc["(Intercept)"] + lc["age10"]*dm$age10 + lc["male"]*dm$male + lc["tg_log"]*dm$tg_log +
            lc["tc"]*dm$tc_mmol + lc["hdl_c"]*dm$hdl_mmol + lc["ldl_c"]*dm$ldl_mmol + lc["apo_b"]*dm$apob_gl
mods <- list(core = plogis(lp_core), std = plogis(lp_std), lasso = plogis(lp_lasso))

## ---------------- 1) weighted external validation, 3 transferred versions ----------------
cat("\n==== weighted external validation (WTMEC2YR/2, ids=~1) ====\n")
extw <- data.frame()
for (nm in names(mods)) {
  p <- mods[[nm]]; y <- dm$any_dr; w <- dm$wt4
  a_w <- aucw(p, y, w); cs_w <- calibw(p, y, w); br_w <- brierw(p, y, w)
  a_u <- auc_fun(p, y); br_u <- brier_fun(p, y)
  extw <- rbind(extw, data.frame(model = nm,
    auc_unweighted = a_u, auc_weighted = a_w,
    brier_unweighted = br_u, brier_weighted = br_w,
    calib_slope_w = cs_w["slope"], calib_intercept_w = cs_w["intercept"]))
  cat(sprintf("%-6s AUC u=%.3f w=%.3f | Brier u=%.4f w=%.4f | slope_w=%.3f int_w=%.3f\n",
      nm, a_u, a_w, br_u, br_w, cs_w["slope"], cs_w["intercept"]))
}
write.csv(extw, file.path(OUT, "v2_phase7h_weighted_ext.csv"), row.names = FALSE)

## ---------------- 2) full reference model OOF (survey-weighted) + lipid-only ----------------
cat("\n==== full reference model, survey-weighted 10-fold OOF (new fold seed 20260906) ====\n")
ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR / 2
fmla <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
        smoker_100 + insulin_use + oral_med + sys_bp
fit_fun <- function(dat) svyglm(fmla, design = svydesign(ids = ~1, weights = ~wt4, data = dat),
                                family = binomial())
set.seed(20260906)   # new fold seed (phase4 used 20260905 stream); documented in log
fold <- sample(rep(1:10, length.out = nrow(dref)))
p_ref_oof <- rep(NA_real_, nrow(dref))
for (k in 1:10) {
  tr <- dref[fold != k, ]; te <- dref[fold == k, ]
  fk <- fit_fun(tr)
  p_ref_oof[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% as.numeric(coef(fk))))
}
lp_lip <- cf[1] + cf["age10"]*dref$age10 + cf["male"]*dref$male +
          cf["tc"]*dref$tc_mmol + cf["hdl_c"]*dref$hdl_mmol
p_lip <- plogis(lp_lip)
ok <- !is.na(p_ref_oof) & !is.na(p_lip)
cat(sprintf("comparable subset n=%d events=%d\n", sum(ok), sum(dref$any_dr[ok])))
w_ok <- dref$wt4[ok]; y_ok <- dref$any_dr[ok]
full_cmp <- data.frame(
  metric = c("AUC_full_oof_w", "AUC_lipid_only_w", "Brier_full_oof_w", "Brier_lipid_only_w",
             "AUC_full_oof_unw", "AUC_lipid_only_unw", "Brier_full_oof_unw", "Brier_lipid_only_unw"),
  value = c(aucw(p_ref_oof[ok], y_ok, w_ok), aucw(p_lip[ok], y_ok, w_ok),
            brierw(p_ref_oof[ok], y_ok, w_ok), brierw(p_lip[ok], y_ok, w_ok),
            auc_fun(p_ref_oof[ok], y_ok), auc_fun(p_lip[ok], y_ok),
            brier_fun(p_ref_oof[ok], y_ok), brier_fun(p_lip[ok], y_ok)))
cs_ref <- calibw(p_ref_oof[ok], y_ok, w_ok)
cs_lip <- calibw(p_lip[ok], y_ok, w_ok)
cat("\n---- weighted full vs lipid-only (comparable subset) ----\n"); print(full_cmp)
cat(sprintf("weighted calib: full slope=%.3f int=%.3f | lipid slope=%.3f int=%.3f\n",
    cs_ref["slope"], cs_ref["intercept"], cs_lip["slope"], cs_lip["intercept"]))
write.csv(full_cmp, file.path(OUT, "v2_phase7h_full_vs_lipid_weighted.csv"), row.names = FALSE)
write.csv(data.frame(model = c("full_oof", "lipid_only"),
                     calib_slope_w = c(cs_ref["slope"], cs_lip["slope"]),
                     calib_intercept_w = c(cs_ref["intercept"], cs_lip["intercept"])),
          file.path(OUT, "v2_phase7h_calib_weighted.csv"), row.names = FALSE)

## ---------------- 3) pmsampsize (Riley) ----------------
cat("\n==== pmsampsize (Riley, binary) ====\n")
# (a) Chinese lipid model: analysis population p=275/1966, S=9, expected C=0.67 (optimism-corrected 0.669)
ps1 <- pmsampsize(type = "b", rsquared = 0.10, parameters = 9,
                  prevalence = 275 / 1966, cstatistic = 0.67, shrinkage = 0.9)
cat("\n-- (a) Chinese lipid model (p=0.140, S=9, C=0.67, s=0.90) --\n"); print(ps1)
# (b) NHANES full reference model: p=267/888, S=10, expected C=0.735
ps2 <- pmsampsize(type = "b", rsquared = 0.15, parameters = 10,
                  prevalence = 267 / 888, cstatistic = 0.735, shrinkage = 0.9)
cat("\n-- (b) NHANES full model (p=0.301, S=10, C=0.735, s=0.90) --\n"); print(ps2)
grab_pms <- function(ps) {
  data.frame(
    sample_size_required = ps$sample_size,
    final_shrinkage = ps$final_shrinkage,
    final_EPP = ps$EPP,
    n_events_required = ps$events,
    stringsAsFactors = FALSE)
}
pms_out <- rbind(
  cbind(model = "chinese_lipid", prevalence = 275/1966, n_parameters = 9,
        expected_C = 0.67, grab_pms(ps1)),
  cbind(model = "nhanes_full", prevalence = 267/888, n_parameters = 10,
        expected_C = 0.735, grab_pms(ps2)))
write.csv(pms_out, file.path(OUT, "v2_phase7h_pmsampsize.csv"), row.names = FALSE)

cat("\n==== ALL DONE (7h) ====\n")
