## AUTO-GENERATED v2 from phase4b_temporal.R (corrected outcome)

# Phase 4b: temporal external validation within NHANES (2005-06 <-> 2007-08)
# Purpose: compensate for inability to validate the full model in the Chinese cohort
#   (Chinese cohort lacks HbA1c/duration/medication -> linkage attempt yielded only n=10).
# Compares robustness of FULL model (HbA1c+duration+insulin+...) vs LIPID-ONLY model.
suppressMessages({ library(rms); library(survey) })
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out"
NH  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
brier_fun <- function(p, y) { ok <- !is.na(p) & !is.na(y); mean((p[ok] - y[ok])^2) }
calib_slope_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  lp <- qlogis(pmin(pmax(p[ok], 1e-6), 1 - 1e-6))
  f <- glm(y[ok] ~ lp, family = binomial())
  c(intercept = as.numeric(coef(f)[1]), slope = as.numeric(coef(f)[2]))
}

## ---- load + unit conversion ----
nh <- read.csv(NH)
nh$tc_mmol  <- nh$tc  / 38.67
nh$hdl_mmol <- nh$hdl / 38.67
nh$tg_mmol  <- nh$tg  / 88.57
nh$ldl_mmol <- nh$ldl / 38.67
nh$apob_gl  <- nh$apob / 100
nh$tg_log   <- log2(nh$tg_mmol)
nh$age10    <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]
cyc <- sort(unique(dm[["X_cycle"]]))
cat("cycles:", cyc, "\n")
for (c in cyc) {
  g <- dm[dm[["X_cycle"]] == c, ]
  cat(sprintf("  %s: n=%d any_dr=%d (%.1f%%)\n", c, nrow(g), sum(g$any_dr), 100*mean(g$any_dr)))
}

fmla_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
             smoker_100 + insulin_use + oral_med + sys_bp
fmla_lip  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol

## ---- temporal validation (train one cycle -> test the other) ----
tv <- function(fmla, label) {
  res <- data.frame()
  for (a in cyc) for (b in cyc) {
    if (a == b) next
    v <- all.vars(fmla)
    tra <- dm[dm[["X_cycle"]] == a & complete.cases(dm[, v]), ]
    te  <- dm[dm[["X_cycle"]] == b & complete.cases(dm[, v]), ]
    dsn <- svydesign(ids = ~1, weights = ~WTMEC2YR, data = tra)   # 2-year weight per cycle
    f <- svyglm(fmla, design = dsn, family = binomial())
    p <- plogis(as.numeric(model.matrix(fmla, data = te) %*% coef(f)))
    cs <- calib_slope_fun(p, te$any_dr)
    res <- rbind(res, data.frame(model = label, train = a, test = b,
      n_train = nrow(tra), n_test = nrow(te), events_test = sum(te$any_dr),
      auc = auc_fun(p, te$any_dr), brier = brier_fun(p, te$any_dr),
      calib_slope = cs["slope"], calib_intercept = cs["intercept"]))
  }
  res
}

cat("\n==== temporal validation: FULL model ====\n")
r_full <- tv(fmla_full, "full")
print(r_full)
cat("\n==== temporal validation: LIPID-ONLY model (fitted within NHANES) ====\n")
r_lip <- tv(fmla_lip, "lipid_only")
print(r_lip)
cat(sprintf("\nmean temporal AUC: full=%.3f | lipid-only=%.3f\n",
    mean(r_full$auc), mean(r_lip$auc)))

## ---- Chinese coefficients applied to each cycle (cross-period stability) ----
cat("\n==== Chinese 'core' coefficients applied to each NHANES cycle ====\n")
dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
cf <- coef(m_core)
r_cn <- data.frame()
for (c in cyc) {
  g <- dm[dm[["X_cycle"]] == c, ]
  lp <- cf[1] + cf["age10"]*g$age10 + cf["male"]*g$male +
        cf["tc"]*g$tc_mmol + cf["hdl_c"]*g$hdl_mmol
  p <- plogis(lp)
  cs <- calib_slope_fun(p, g$any_dr)
  r_cn <- rbind(r_cn, data.frame(model = "chinese_core", cycle = c,
    n = sum(!is.na(p)), events = sum(g$any_dr[!is.na(p)]),
    auc = auc_fun(p, g$any_dr), brier = brier_fun(p, g$any_dr),
    calib_slope = cs["slope"], calib_intercept = cs["intercept"]))
}
print(r_cn)

all_res <- rbind(r_full, r_lip)
write.csv(all_res, file.path(OUT, "v2_phase4b_temporal.csv"), row.names = FALSE)
write.csv(r_cn, file.path(OUT, "v2_phase4b_chinese_by_cycle.csv"), row.names = FALSE)

cat("\n==== ALL DONE ====\n")
