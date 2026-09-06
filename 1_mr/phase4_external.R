# Phase 4: external validation of Chinese bedside lipid DR score on NHANES (set B)
# 4.1 full reference model | 4.2 lipid-only score | 4.3 dAUC/NRI/IDI | 4.4 cross-ethnicity shift
# Pure ASCII; all metrics self-implemented (no pROC / PredictABEL dependency)
suppressMessages({ library(rms) })   # rms for cut2 only
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv"
NH  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv"
LASSO <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out/phase2_lasso.csv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "phase4_external.log"), split = TRUE)

has_survey <- requireNamespace("survey", quietly = TRUE)
cat("survey package available:", has_survey, "\n")

## ---------------- helpers ----------------
auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_ci <- function(p, y, B = 1000) {
  ok <- !is.na(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]; n <- length(y)
  out <- replicate(B, auc_fun(p[sample(n, n, TRUE)], y[sample(n, n, TRUE)]))
  quantile(out, c(0.025, 0.975), na.rm = TRUE)
}
brier_fun <- function(p, y) { ok <- !is.na(p) & !is.na(y); mean((p[ok] - y[ok])^2) }
calib_tbl <- function(p, y, g = 10) {
  ok <- !is.na(p) & !is.na(y)
  p <- p[ok]; y <- y[ok]
  grp <- cut2(p, g = g)
  data.frame(pred_mean = as.numeric(tapply(p, grp, mean)),
             obs_rate  = as.numeric(tapply(y, grp, mean)),
             n         = as.numeric(tapply(y, grp, length)))
}
calib_slope_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  lp <- qlogis(pmin(pmax(p[ok], 1e-6), 1 - 1e-6))
  f <- glm(y[ok] ~ lp, family = binomial())
  c(intercept = as.numeric(coef(f)[1]), slope = as.numeric(coef(f)[2]))
}
nri_idi <- function(p0, p1, y, cut) {   # p0=old(null/lipid), p1=new(full)
  ev <- y == 1
  idi <- (mean(p1[ev], na.rm = TRUE) - mean(p0[ev], na.rm = TRUE)) -
         (mean(p1[!ev], na.rm = TRUE) - mean(p0[!ev], na.rm = TRUE))
  up_e <- mean(p1[ev] >= cut & p0[ev] < cut); down_e <- mean(p1[ev] < cut & p0[ev] >= cut)
  up_n <- mean(p1[!ev] >= cut & p0[!ev] < cut); down_n <- mean(p1[!ev] < cut & p0[!ev] >= cut)
  c(NRI = up_e - down_e + down_n - up_n, IDI = idi)
}

## ---------------- load NHANES + unit conversion ----------------
cat("\n==== NHANES set B: load + unit conversion ====\n")
nh <- read.csv(NH)
# NHANES reports lipids in mg/dL; Chinese cohort in mmol/L (apo in g/L)
nh$tc_mmol  <- nh$tc  / 38.67
nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67
nh$tg_mmol  <- nh$tg  / 88.57
nh$apob_gl  <- nh$apob / 100      # ApoB mg/dL -> g/L
nh$tg_log   <- log2(nh$tg_mmol)
nh$age10    <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]
cat(sprintf("diabetes subgroup n=%d  any_dr=%d (%.1f%%)\n", nrow(dm), sum(dm$any_dr), 100*mean(dm$any_dr)))
av <- function(...) sum(complete.cases(dm[, c(...)]))
cat(sprintf("available: TC+HDL=%d | +TG+LDL=%d | +ApoB=%d\n",
    av("tc_mmol","hdl_mmol"),
    av("tc_mmol","hdl_mmol","tg_mmol","ldl_mmol"),
    av("tc_mmol","hdl_mmol","tg_mmol","ldl_mmol","apob_gl")))

## ---------------- fit transportable models on Chinese dev cohort ----------------
cat("\n==== fit models on Chinese development cohort ====\n")
dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
cat(sprintf("Chinese dev n=%d events=%d (%.1f%%)\n", nrow(d), sum(d$any_dr), 100*mean(d$any_dr)))
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
m_std  <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
cat("\n-- core (age+sex+TC+HDL) --\n"); print(summary(m_core)$coefficients)
cat("\n-- std (core+TG+LDL) --\n"); print(summary(m_std)$coefficients)
lasso <- read.csv(LASSO); lc <- setNames(lasso$coef, lasso$var)
cat("\n-- LASSO (phase 2, includes ApoB) --\n"); print(lasso)

# in-sample AUC on Chinese dev (for head-to-head reference)
auc_dev <- data.frame(
  model = c("core","std","lasso"),
  auc_chinese_dev = c(auc_fun(predict(m_core, type="response"), d$any_dr),
                      auc_fun(predict(m_std,  type="response"), d$any_dr),
                      NA))

## ---------------- 4.2 + 4.4 external validation ----------------
cat("\n==== 4.2/4.4 external validation on NHANES diabetes subgroup ====\n")
cf <- coef(m_core)
lp_core <- cf[1] + cf["age10"]*dm$age10 + cf["male"]*dm$male +
           cf["tc"]*dm$tc_mmol + cf["hdl_c"]*dm$hdl_mmol
cf2 <- coef(m_std)
lp_std <- cf2[1] + cf2["age10"]*dm$age10 + cf2["male"]*dm$male +
          cf2["tc"]*dm$tc_mmol + cf2["hdl_c"]*dm$hdl_mmol +
          cf2["tg_log"]*dm$tg_log + cf2["ldl_c"]*dm$ldl_mmol
lp_lasso <- lc["(Intercept)"] + lc["age10"]*dm$age10 + lc["male"]*dm$male +
            lc["tg_log"]*dm$tg_log + lc["tc"]*dm$tc_mmol + lc["hdl_c"]*dm$hdl_mmol +
            lc["ldl_c"]*dm$ldl_mmol + lc["apo_b"]*dm$apob_gl
mods <- list(core = plogis(lp_core), std = plogis(lp_std), lasso = plogis(lp_lasso))

ext <- data.frame()
for (nm in names(mods)) {
  p <- mods[[nm]]; y <- dm$any_dr
  ok <- !is.na(p) & !is.na(y)
  a <- auc_fun(p, y); ci <- auc_ci(p, y, B = 1000); cs <- calib_slope_fun(p, y)
  ext <- rbind(ext, data.frame(model = nm, n = sum(ok), events = sum(y[ok]),
    auc = a, auc_lo = ci[1], auc_hi = ci[2], brier = brier_fun(p, y),
    calib_slope = cs["slope"], calib_intercept = cs["intercept"],
    auc_chinese_dev = auc_dev$auc_chinese_dev[match(nm, auc_dev$model)]))
  cat(sprintf("%-6s n=%4d ev=%3d AUC=%.3f (%.3f-%.3f) Brier=%.4f slope=%.3f int=%.3f\n",
      nm, sum(ok), sum(y[ok]), a, ci[1], ci[2], brier_fun(p, y), cs["slope"], cs["intercept"]))
  write.csv(calib_tbl(p, y), file.path(OUT, paste0("phase4_calib_", nm, ".csv")), row.names = FALSE)
}
print(ext)
write.csv(ext, file.path(OUT, "phase4_external_auc.csv"), row.names = FALSE)

## ---------------- 4.1 full reference model (NHANES) ----------------
cat("\n==== 4.1 full reference model (NHANES diabetes subgroup) ====\n")
ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR / 2
cat(sprintf("reference n=%d events=%d (%.1f%%)\n", nrow(dref), sum(dref$any_dr), 100*mean(dref$any_dr)))
fmla <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
        smoker_100 + insulin_use + oral_med + sys_bp
fit_fun <- function(dat) {
  if (has_survey) {
    survey::svyglm(fmla, design = survey::svydesign(ids = ~1, weights = ~wt4, data = dat),
                   family = binomial())
  } else glm(fmla, data = dat, family = binomial())
}
fit_ref <- fit_fun(dref)
cat(if (has_survey) "(survey-weighted)\n" else "(unweighted - survey unavailable)\n")
b <- coef(fit_ref); se <- sqrt(diag(vcov(fit_ref)))
ref_tab <- data.frame(var = names(b), coef = as.numeric(b), OR = exp(as.numeric(b)),
                      p = 2*pnorm(-abs(as.numeric(b)/se)))
print(ref_tab)
write.csv(ref_tab, file.path(OUT, "phase4_reference_model.csv"), row.names = FALSE)

# out-of-fold predictions (fair comparison vs lipid-only)
fold <- sample(rep(1:10, length.out = nrow(dref)))
p_ref_oof <- rep(NA_real_, nrow(dref))
for (k in 1:10) {
  tr <- dref[fold != k, ]; te <- dref[fold == k, ]
  fk <- fit_fun(tr)
  p_ref_oof[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% as.numeric(coef(fk))))
}
auc_ref_oof <- auc_fun(p_ref_oof, dref$any_dr)
cat(sprintf("full reference model AUC (10-fold out-of-fold) = %.3f\n", auc_ref_oof))

## ---------------- 4.3 dAUC / NRI / IDI ----------------
cat("\n==== 4.3 lipid-only vs full model increment (NHANES) ====\n")
lp_lip <- cf[1] + cf["age10"]*dref$age10 + cf["male"]*dref$male +
          cf["tc"]*dref$tc_mmol + cf["hdl_c"]*dref$hdl_mmol
p_lip <- plogis(lp_lip)
ok <- !is.na(p_lip) & !is.na(p_ref_oof)
cat(sprintf("comparable subset n=%d events=%d\n", sum(ok), sum(dref$any_dr[ok])))
auc_lip <- auc_fun(p_lip[ok], dref$any_dr[ok])
auc_ref <- auc_fun(p_ref_oof[ok], dref$any_dr[ok])
ni <- nri_idi(p_lip[ok], p_ref_oof[ok], dref$any_dr[ok], cut = median(p_ref_oof[ok], na.rm = TRUE))
cat(sprintf("lipid-only AUC=%.3f | full-model AUC=%.3f | dAUC=%.3f\n", auc_lip, auc_ref, auc_ref - auc_lip))
cat(sprintf("NRI=%.4f  IDI=%.4f (cut=median full-model risk)\n", ni["NRI"], ni["IDI"]))
inc <- data.frame(metric = c("n","events","AUC_lipid_only","AUC_full_oof","dAUC","NRI","IDI"),
                  value = c(sum(ok), sum(dref$any_dr[ok]), auc_lip, auc_ref, auc_ref - auc_lip,
                            ni["NRI"], ni["IDI"]))
write.csv(inc, file.path(OUT, "phase4_increment.csv"), row.names = FALSE)
print(inc)

## ---------------- calibration plot ----------------
png(file.path(OUT, "phase4_calibration_plot.png"), width = 1000, height = 420, res = 120)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.6, 1))
for (nm in names(mods)) {
  ct <- calib_tbl(mods[[nm]], dm$any_dr)
  lim <- max(c(ct$pred_mean, ct$obs_rate), na.rm = TRUE) * 1.1
  plot(ct$pred_mean, ct$obs_rate, xlim = c(0, lim), ylim = c(0, lim),
       xlab = "Predicted risk", ylab = "Observed rate",
       main = paste0("NHANES external: ", nm), pch = 19, col = "#2C7FB8")
  abline(0, 1, lty = 2, col = "red", lwd = 1.5)
  lines(ct$pred_mean, ct$obs_rate, col = "#2C7FB8", lwd = 1.2)
}
dev.off()
cat("\ncalibration plot written\n")

cat("\n==== ALL DONE ====\n")
sink()
