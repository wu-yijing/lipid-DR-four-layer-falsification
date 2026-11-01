## =====================================================================
## revision_20261031: provenance and correction of Supplementary Table S42
##
## What this script establishes
##  (1) S42 was computed on the original outcome definition (1,119 diabetes
##      subgroup, ungradable photographs counted as non-cases): the complete-case
##      n and event counts of Wu 6/8 (910/268), Wu 4-variable (927/272) and
##      Guo (462/146) reproduce exactly on that frame.
##  (2) The Wu 6/8 AUC printed in S42 (0.531) is a unit error: HDL-C was entered
##      in mg/dL although the published coefficient (-0.635) is per mmol/L.
##      Passing mg/dL reproduces 0.531 to three decimals; passing mmol/L gives
##      0.707 on the same 910 participants. The 4-variable version, which has no
##      HDL term, reproduces the printed 0.719 exactly - so the pipeline is
##      otherwise sound.
##  (3) Guo 6/7 reproduces (0.708 vs 0.709 printed; 0.731 vs 0.727 with the
##      triglyceride component zeroed).
##  (4) Gao 2025 cannot be re-evaluated: only hazard ratios are published and
##      they are not openly accessible, and the reconstructable subset (BUN,
##      HbA1c, duration) is 1,012 participants against the 998 printed.
##
## Primary-frame values (884 with ascertainable retinal status) are reported for
## the two models that can be recomputed.
##
## Run: Rscript 4_external/revision_20261031_tableS42_primary_frame.R
## =====================================================================

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- if (is.na(this)) "." else normalizePath(file.path(dirname(this), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "output_tables"))) REPO <- "."
OUT <- file.path(REPO, "output_tables")

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

nh0 <- read.csv(file.path(SRC, "phase2_ext/NHANES_B.csv"))
nh0$hdl_mmol <- nh0$hdl / 38.67
nh0$tg_mmol  <- nh0$tg / 88.57

## Wu 2026, Table 3; intercept -12.45 is irrelevant to the AUC
wu_lp <- function(d, hdl, egfr = TRUE, age_dur_sbp_a1c_only = FALSE) {
  lp <- -12.450 + 0.045 * d$age + 0.174 * d$dm_duration + 0.031 * d$sys_bp +
        0.445 * d$hba1c_pct
  if (age_dur_sbp_a1c_only) return(lp)
  lp - 0.635 * hdl - 0.022 * d$egfr_ckdepi2021
}

## Guo 2024, Table 2 point score (albuminuria unavailable -> 0)
guo_score <- function(d, tg = TRUE) {
  s <- ifelse(d$male == 1, 2.5, 0) +
       ifelse(!is.na(d$smoker_100) & d$smoker_100 == 1, 4.5, 0)
  dur <- d$dm_duration
  s <- s + ifelse(dur < 5, 0, ifelse(dur < 10, 2.5, ifelse(dur < 15, 5,
           ifelse(dur < 20, 7.5, 10))))
  h <- d$hba1c_pct
  s <- s + ifelse(h < 7, 0, ifelse(h < 8, 2, ifelse(h < 9, 4, 6)))
  b <- d$sys_bp
  s <- s + ifelse(b < 130, 0, ifelse(b < 140, 1.5, ifelse(b < 150, 3, 4.5)))
  if (tg) s <- s + ifelse(d$tg_mmol >= 1.70, 4, 0)
  s
}

frame <- function(primary) {
  d <- nh0[nh0$diabetes == 1, ]
  if (primary) d <- d[!is.na(d$dr_level4), ]
  d
}

tab <- data.frame()
for (prim in c(FALSE, TRUE)) {
  d <- frame(prim)
  fr <- if (prim) "primary (884)" else "original (1,119)"

  k6 <- complete.cases(d[, c("age", "dm_duration", "sys_bp", "hba1c_pct", "hdl_mmol", "egfr_ckdepi2021")])
  tab <- rbind(tab, data.frame(
    frame = fr, model = "Wu 6 of 8 (HDL in mmol/L, corrected)",
    n = sum(k6), events = sum(d$any_dr[k6]),
    auc = auc_fun(wu_lp(d[k6, ], d$hdl_mmol[k6]), d$any_dr[k6])))

  tab <- rbind(tab, data.frame(
    frame = fr, model = "Wu 6 of 8 (HDL in mg/dL, as printed in S42)",
    n = sum(k6), events = sum(d$any_dr[k6]),
    auc = auc_fun(wu_lp(d[k6, ], d$hdl[k6]), d$any_dr[k6])))

  k4 <- complete.cases(d[, c("age", "dm_duration", "sys_bp", "hba1c_pct")])
  tab <- rbind(tab, data.frame(
    frame = fr, model = "Wu 4 variables (age, duration, SBP, HbA1c)",
    n = sum(k4), events = sum(d$any_dr[k4]),
    auc = auc_fun(wu_lp(d[k4, ], NULL, age_dur_sbp_a1c_only = TRUE), d$any_dr[k4])))

  kg <- complete.cases(d[, c("male", "smoker_100", "dm_duration", "hba1c_pct", "sys_bp", "tg_mmol")])
  tab <- rbind(tab, data.frame(
    frame = fr, model = "Guo 6 of 7",
    n = sum(kg), events = sum(d$any_dr[kg]),
    auc = auc_fun(guo_score(d[kg, ]), d$any_dr[kg])))

  tab <- rbind(tab, data.frame(
    frame = fr, model = "Guo 6 of 7, TG component set to 0",
    n = sum(kg), events = sum(d$any_dr[kg]),
    auc = auc_fun(guo_score(d[kg, ], tg = FALSE), d$any_dr[kg])))
}

tab$auc <- round(tab$auc, 3)
print(tab, row.names = FALSE)

cat("\nprinted in S42      : Wu6 0.531 | Wu4 0.719 | Guo 0.709 | GuoTG0 0.727\n")
cat("reproduced (orig.)  : Wu6 0.531 (mg/dL) / 0.707 (mmol/L) | Wu4 0.719 | Guo 0.708 | GuoTG0 0.731\n")

write.csv(tab, file.path(OUT, "tableS42_primary_frame.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT, "tableS42_primary_frame.csv"), "\n")
