## AUTO-GENERATED v2 from phase4d_mi_export.R (corrected outcome, option A)

# Phase 4d: sensitivity — export MI-pooled coefficients (m=50, maxit=50) to NHANES external validation
# Mirrors the CC-fitted core/std/LASSO external validation in phase4_external.R, but with
# coefficients pooled from multiple imputation of the FULL gradeable cohort (n=1,971).
# Unweighted + MEC-weighted AUC / calibration slope, to compare against CC export.
suppressMessages({ library(mice); library(survey); library(rms) })
options(survey.lonely.psu = "adjust")
set.seed(20260906)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------------- helpers (replicate phase7h) ----------------
auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aucw <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y) & !is.na(w); p <- p[ok]; y <- y[ok]; w <- w[ok]
  i1 <- which(y == 1); i0 <- which(y == 0)
  W1 <- sum(w[i1]); W0 <- sum(w[i0])
  if (W1 == 0 || W0 == 0) return(NA_real_)
  num <- 0
  for (i in i1) num <- num + w[i] * sum(w[i0] * (p[i] > p[i0]) + 0.5 * w[i0] * (p[i] == p[i0]))
  num / (W1 * W0)
}
calib_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y)
  lp <- qlogis(pmin(pmax(p[ok], 1e-6), 1 - 1e-6))
  f <- glm(y[ok] ~ lp, family = binomial())
  c(intercept = as.numeric(coef(f)[1]), slope = as.numeric(coef(f)[2]))
}
calibw <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  lp <- qlogis(pmin(pmax(p[ok], 1e-6), 1 - 1e-6))
  dd <- data.frame(y = y[ok], lp = lp, w = w[ok])
  f <- svyglm(y ~ lp, design = svydesign(ids = ~1, weights = ~w, data = dd), family = binomial())
  c(intercept = as.numeric(coef(f)[1]), slope = as.numeric(coef(f)[2]))
}

## ---------------- MI on full gradeable Chinese cohort ----------------
dev <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv")
g <- dev[dev$gradeable == 1, ]
d <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
d$tg_log <- log2(d$tg); d$tg <- NULL; d$age10 <- d$age / 10
cat(sprintf("MI dataset: n=%d (gradeable), events=%d, complete-cases=%d\n",
            nrow(d), sum(d$any_dr), sum(complete.cases(d))))

imp <- mice(d, m = 50, method = "pmm", maxit = 50, printFlag = FALSE, seed = 20260906)
cat("imputation done (m=50, maxit=50)\n")

# fit the three exported versions on each imputed dataset, pool coefficients (Rubin)
f_core <- with(imp, glm(any_dr ~ age10 + male + tc + hdl_c, family = binomial()))
f_std  <- with(imp, glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, family = binomial()))
f_full <- with(imp, glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_b, family = binomial()))

cf_core <- setNames(summary(pool(f_core))$estimate, summary(pool(f_core))$term)
cf_std  <- setNames(summary(pool(f_std))$estimate,  summary(pool(f_std))$term)
cf_full <- setNames(summary(pool(f_full))$estimate, summary(pool(f_full))$term)
cat("\n== pooled MI coefficients (log-odds) ==\n")
print(cf_core); print(cf_std); print(cf_full)

## ---------------- load NHANES + unit conversion (replicate phase4) ----------------
nh <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv")
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$apob_gl <- nh$apob / 100
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]
dm$wt4 <- dm$WTMEC2YR / 2
cat(sprintf("NHANES diabetes n=%d events=%d (%.1f%%)\n", nrow(dm), sum(dm$any_dr), 100*mean(dm$any_dr)))

lp_core <- cf_core["(Intercept)"] + cf_core["age10"]*dm$age10 + cf_core["male"]*dm$male +
           cf_core["tc"]*dm$tc_mmol + cf_core["hdl_c"]*dm$hdl_mmol
lp_std  <- cf_std["(Intercept)"] + cf_std["age10"]*dm$age10 + cf_std["male"]*dm$male +
           cf_std["tc"]*dm$tc_mmol + cf_std["hdl_c"]*dm$hdl_mmol +
           cf_std["tg_log"]*dm$tg_log + cf_std["ldl_c"]*dm$ldl_mmol
lp_full <- cf_full["(Intercept)"] + cf_full["age10"]*dm$age10 + cf_full["male"]*dm$male +
           cf_full["tg_log"]*dm$tg_log + cf_full["tc"]*dm$tc_mmol + cf_full["hdl_c"]*dm$hdl_mmol +
           cf_full["ldl_c"]*dm$ldl_mmol + cf_full["apo_b"]*dm$apob_gl
mods <- list(core_mi = plogis(lp_core), std_mi = plogis(lp_std), full_mi = plogis(lp_full))

cat("\n==== MI-pooled export to NHANES (unweighted + weighted) ====\n")
res <- data.frame()
for (nm in names(mods)) {
  p <- mods[[nm]]; y <- dm$any_dr; w <- dm$wt4
  a_u <- auc_fun(p, y); a_w <- aucw(p, y, w)
  cs_u <- calib_fun(p, y); cs_w <- calibw(p, y, w)
  n_ok <- sum(!is.na(p) & !is.na(y)); ev <- sum(y[!is.na(p) & !is.na(y)])
  res <- rbind(res, data.frame(model = nm, n = n_ok, events = ev,
    auc_unw = a_u, auc_w = a_w, slope_unw = cs_u["slope"], slope_w = cs_w["slope"],
    int_unw = cs_u["intercept"], int_w = cs_w["intercept"]))
  cat(sprintf("%-8s n=%4d ev=%3d AUC u=%.3f w=%.3f | slope u=%.3f w=%.3f\n",
      nm, n_ok, ev, a_u, a_w, cs_u["slope"], cs_w["slope"]))
}
print(res)
write.csv(res, file.path(OUT, "phase4d_mi_export.csv"), row.names = FALSE)

cat("\n==== CC export (from phase4) for comparison ====\n")
cat("core   n=1041 ev=302 AUC u=0.513 w=0.537 | slope u=0.309 w=0.791\n")
cat("std    n= 494 ev=148 AUC u=0.530 w=0.572 | slope u=0.378 w=0.935\n")
cat("lasso  n= 190 ev= 64 AUC u=0.477 w=0.468 | slope u=-0.116 w=-0.167\n")

cat("\n==== ALL DONE (4d) ====\n")
