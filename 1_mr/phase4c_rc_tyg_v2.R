## AUTO-GENERATED v2 from phase4c_rc_tyg.R (corrected outcome, option A)

# Improvement 1: remnant cholesterol (RC) and TyG analyses in both cohorts
# Tests whether emerging lipid metrics escape the "lipid dimension is weak" conclusion
suppressMessages({ library(rms) })
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---------------- NHANES ----------------
cat("==== NHANES diabetes subgroup ====\n")
nh <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv")
nh$tc_mmol  <- nh$tc  / 38.67
nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67
nh$tg_mmol  <- nh$tg  / 88.57
# RC = TC - LDL - HDL (mmol/L)
nh$rc_mmol <- nh$tc_mmol - nh$ldl_mmol - nh$hdl_mmol
# TyG = ln[TG(mg/dL) x FPG(mg/dL) / 2]
nh$tyg <- log(nh$tg * nh$glucose / 2)
nh$tyg_bmi <- nh$tyg * nh$bmi
nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]

mets <- c("tc_mmol", "hdl_mmol", "ldl_mmol", "tg_mmol", "rc_mmol", "tyg", "tyg_bmi")
mlab <- c("TC (mmol/L)", "HDL-C (mmol/L)", "LDL-C (mmol/L)", "TG (mmol/L)",
          "Remnant chol. (mmol/L)", "TyG index", "TyG x BMI")
res_n <- data.frame()
for (i in seq_along(mets)) {
  x <- mets[i]; nm <- mlab[i]
  dx <- dm[!is.na(dm[[x]]) & !is.na(dm$any_dr) & !is.na(dm$age10) & !is.na(dm$male), ]
  sd_x <- sd(dx[[x]], na.rm = TRUE)
  # crude
  f0 <- glm(any_dr ~ dx[[x]], data = dx, family = binomial())
  or0 <- exp(coef(f0)[2]); p0 <- summary(f0)$coefficients[2, 4]
  auc0 <- auc_fun(dx[[x]], dx$any_dr)
  # adjusted for age+sex
  f1 <- glm(any_dr ~ dx[[x]] + age10 + male, data = dx, family = binomial())
  or1 <- exp(coef(f1)[2]); p1 <- summary(f1)$coefficients[2, 4]
  auc1 <- auc_fun(predict(f1, type = "response"), dx$any_dr)
  # base (age+sex) AUC on same subset for delta
  fb <- glm(any_dr ~ age10 + male, data = dx, family = binomial())
  aucb <- auc_fun(predict(fb, type = "response"), dx$any_dr)
  res_n <- rbind(res_n, data.frame(cohort = "NHANES", metric = nm, n = nrow(dx),
    events = sum(dx$any_dr), OR_perSD_crude = or0, p_crude = p0,
    OR_perSD_adj = or1, p_adj = p1, AUC_single = auc0,
    AUC_base_age_sex = aucb, AUC_base_plus_metric = auc1,
    deltaAUC = auc1 - aucb))
  cat(sprintf("%-26s n=%4d ev=%3d  ORadj=%.2f (p=%.3g)  AUC=%.3f  dAUC=%.3f\n",
      nm, nrow(dx), sum(dx$any_dr), or1, p1, auc0, auc1 - aucb))
}
write.csv(res_n, file.path(OUT, "phase4c_nhanes_rc_tyg.csv"), row.names = FALSE)

## ---------------- Chinese cohort ----------------
cat("\n==== Chinese screening cohort (gradable) ====\n")
dev <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv")
g <- dev[dev$gradeable == 1, ]
g$rc_mmol <- g$tc - g$hdl_c - g$ldl_c     # mmol/L
g$age10 <- g$age / 10
cmap <- c(tc_mmol = "tc", hdl_mmol = "hdl_c", ldl_mmol = "ldl_c", tg_mmol = "tg", rc_mmol = "rc_mmol")
res_c <- data.frame()
for (i in seq_along(mets)) {
  x <- mets[i]; nm <- mlab[i]
  if (x == "tyg" || x == "tyg_bmi") next    # no glucose in Chinese cohort
  xg <- if (x %in% names(cmap)) cmap[[x]] else x
  dx <- g[!is.na(g[[xg]]) & !is.na(g$any_dr) & !is.na(g$age10) & !is.na(g$male), ]
  f0 <- glm(any_dr ~ dx[[xg]], data = dx, family = binomial())
  or0 <- exp(coef(f0)[2]); p0 <- summary(f0)$coefficients[2, 4]
  auc0 <- auc_fun(dx[[xg]], dx$any_dr)
  f1 <- glm(any_dr ~ dx[[xg]] + age10 + male, data = dx, family = binomial())
  or1 <- exp(coef(f1)[2]); p1 <- summary(f1)$coefficients[2, 4]
  auc1 <- auc_fun(predict(f1, type = "response"), dx$any_dr)
  fb <- glm(any_dr ~ age10 + male, data = dx, family = binomial())
  aucb <- auc_fun(predict(fb, type = "response"), dx$any_dr)
  res_c <- rbind(res_c, data.frame(cohort = "Chinese", metric = nm, n = nrow(dx),
    events = sum(dx$any_dr), OR_perSD_crude = or0, p_crude = p0,
    OR_perSD_adj = or1, p_adj = p1, AUC_single = auc0,
    AUC_base_age_sex = aucb, AUC_base_plus_metric = auc1,
    deltaAUC = auc1 - aucb))
  cat(sprintf("%-26s n=%4d ev=%3d  ORadj=%.2f (p=%.3g)  AUC=%.3f  dAUC=%.3f\n",
      nm, nrow(dx), sum(dx$any_dr), or1, p1, auc0, auc1 - aucb))
}
write.csv(res_c, file.path(OUT, "phase4c_chinese_rc.csv"), row.names = FALSE)

cat("\n==== ALL DONE ====\n")
