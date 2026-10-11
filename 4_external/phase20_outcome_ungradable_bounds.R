## =====================================================================
## phase20: NHANES outcome-misclassification bounds for layers 3 and 4
##
## Reviewer item M2. The primary analyses classify the 1,093 participants
## whose fundus photographs could not be graded as non-cases (OPDUARMA /
## dr_level4 missing). 235 of them are in the diabetes subgroup. No
## sensitivity analysis for that decision was reported, although the
## Chinese cohort received three bounding analyses (Table S13).
##
## This script reruns the layer-3 (transported score) and layer-4
## (full reference model) discrimination analyses under three outcome
## definitions:
##   A  primary        : ungradable photographs counted as non-cases
##   B  upper bound    : ungradable photographs counted as cases
##   C  complete case  : ungradable photographs excluded
##
## Variant A reproduces the published values and is asserted against them,
## so the variants are only read once the baseline has been reproduced.
##
## Method is otherwise identical to phase4e_local_v2_dedup.R (same seeds,
## same folds for A and B, same estimators). Pure ASCII.
## Run: Rscript 4_external/phase20_outcome_ungradable_bounds.R
## =====================================================================
suppressMessages({ library(survey) })

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- if (is.na(this)) "." else normalizePath(file.path(dirname(this), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "output_tables"))) REPO <- "."
OUT <- file.path(REPO, "output_tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DEV <- file.path(SRC, "phase2_dev/dev_set_v2.csv")
NH  <- file.path(SRC, "phase2_ext/NHANES_B.csv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

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

## ---------------- Chinese coefficients (unaffected by the variant) ----
dev0 <- read.csv(DEV)
d <- dev0[dev0$gradeable == 1 &
          complete.cases(dev0[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
m_std  <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
cf <- coef(m_core); cs <- coef(m_std)
cat(sprintf("Chinese models fitted on n=%d (events=%d)\n", nrow(d), sum(d$any_dr)))

nh0 <- read.csv(NH)
ungradable <- is.na(nh0$dr_level4)

fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_full_lip <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
                smoker_100 + insulin_use + oral_med + sys_bp + tc_mmol + hdl_mmol
fml_lip_core <- any_dr ~ age10 + male + tc_mmol + hdl_mmol

run_variant <- function(variant) {
  nh <- nh0
  if (variant == "B") nh$any_dr[ungradable] <- 1
  if (variant == "C") nh <- nh[!ungradable, ]
  nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
  nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
  nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
  dm <- nh[nh$diabetes == 1, ]

  ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
                "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
  dref <- dm[complete.cases(dm[, ref_vars]), ]
  dref$wt4 <- dref$WTMEC2YR / 2

  set.seed(20260912)
  fold <- sample(rep(1:10, length.out = nrow(dref)))
  oof <- function(fmla, dat, wcol = NULL) {
    p <- rep(NA_real_, nrow(dat))
    for (k in 1:10) {
      tr <- dat[fold != k, ]; te <- dat[fold == k, ]
      f <- try(if (is.null(wcol)) glm(fmla, data = tr, family = binomial())
               else survey::svyglm(fmla, design = survey::svydesign(ids = ~1,
                      weights = as.formula(paste0("~", wcol)), data = tr), family = binomial()),
               silent = TRUE)
      if (inherits(f, "try-error")) next
      okk <- complete.cases(model.frame(fmla, data = te, na.action = na.pass))
      if (!any(okk)) next
      pr <- as.numeric(predict(f, newdata = te[okk, , drop = FALSE], type = "response"))
      stopifnot(length(pr) == sum(okk))
      p[which(fold == k)[okk]] <- pr
    }
    p
  }
  p_full <- oof(fml_full, dref, "wt4")
  p_full_lip <- oof(fml_full_lip, dref, "wt4")
  p_loc_core <- oof(fml_lip_core, dref, "wt4")
  p_tr_core <- plogis(cf[1] + cf["age10"] * dref$age10 + cf["male"] * dref$male +
                      cf["tc"] * dref$tc_mmol + cf["hdl_c"] * dref$hdl_mmol)

  ok <- !is.na(p_full) & !is.na(p_loc_core) & !is.na(p_tr_core)
  okf <- ok & !is.na(p_full_lip)
  y <- dref$any_dr
  a_full <- auc_fun(p_full[ok], y[ok])
  a_loc  <- auc_fun(p_loc_core[ok], y[ok])
  a_tr   <- auc_fun(p_tr_core[ok], y[ok])
  w_full <- auc_w_fun(p_full[ok], y[ok], dref$wt4[ok])
  w_loc  <- auc_w_fun(p_loc_core[ok], y[ok], dref$wt4[ok])
  w_tr   <- auc_w_fun(p_tr_core[ok], y[ok], dref$wt4[ok])
  a_fl   <- auc_fun(p_full[okf], y[okf]); a_fll <- auc_fun(p_full_lip[okf], y[okf])
  w_fl   <- auc_w_fun(p_full[okf], y[okf], dref$wt4[okf])
  w_fll  <- auc_w_fun(p_full_lip[okf], y[okf], dref$wt4[okf])

  ## transported core / std on the whole diabetes subgroup (as published)
  lp_core <- cf[1] + cf["age10"]*dm$age10 + cf["male"]*dm$male +
             cf["tc"]*dm$tc_mmol + cf["hdl_c"]*dm$hdl_mmol
  lp_std  <- cs[1] + cs["age10"]*dm$age10 + cs["male"]*dm$male + cs["tc"]*dm$tc_mmol +
             cs["hdl_c"]*dm$hdl_mmol + cs["tg_log"]*dm$tg_log + cs["ldl_c"]*dm$ldl_mmol
  pc <- plogis(lp_core); ps <- plogis(lp_std)
  okc <- !is.na(pc) & !is.na(dm$any_dr); oks <- !is.na(ps) & !is.na(dm$any_dr)
  a_core_all <- auc_fun(pc[okc], dm$any_dr[okc])
  a_std_all  <- auc_fun(ps[oks], dm$any_dr[oks])
  set.seed(20260912)
  ci_core <- quantile(replicate(2000, { i <- sample(sum(okc), sum(okc), TRUE)
                       auc_fun(pc[okc][i], dm$any_dr[okc][i]) }), c(0.025, 0.975), na.rm = TRUE)
  set.seed(20260912)
  ci_std <- quantile(replicate(2000, { i <- sample(sum(oks), sum(oks), TRUE)
                      auc_fun(ps[oks][i], dm$any_dr[oks][i]) }), c(0.025, 0.975), na.rm = TRUE)
  n_core <- sum(okc); n_std <- sum(oks)

  data.frame(
    variant = variant,
    dm_n = nrow(dm), dm_events = sum(dm$any_dr),
    dm_prev = mean(dm$any_dr),
    ref_n = nrow(dref), ref_events = sum(dref$any_dr),
    auc_transported_core_all = a_core_all,
    core_all_n = n_core, core_all_lo = ci_core[1], core_all_hi = ci_core[2],
    auc_transported_std_all  = a_std_all,
    std_all_n = n_std, std_all_lo = ci_std[1], std_all_hi = ci_std[2],
    auc_transported_core_ref = a_tr, auc_transported_core_ref_w = w_tr,
    auc_local_lipid_core = a_loc,   auc_local_lipid_core_w = w_loc,
    auc_full = a_full,              auc_full_w = w_full,
    auc_full_plus_lipids = a_fll,   auc_full_plus_lipids_w = w_fll,
    dAUC_full_minus_local = a_full - a_loc,
    dAUC_local_minus_transported = a_loc - a_tr,
    dAUC_full_minus_transported = a_full - a_tr,
    dAUC_lipids_over_full = a_fll - a_fl,
    dAUC_lipids_over_full_w = w_fll - w_fl,
    stringsAsFactors = FALSE)
}

cat("\n==== outcome-definition variants ====\n")
res <- do.call(rbind, lapply(c("A", "B", "C"), run_variant))
print(res, row.names = FALSE, digits = 4)

## ---------------- assert the baseline reproduces the published values ----
A <- res[res$variant == "A", ]
chk <- c(
  core_all = A$auc_transported_core_all,   core_pub = 0.49792,
  std_all  = A$auc_transported_std_all,    std_pub  = 0.49281,
  trans888 = A$auc_transported_core_ref,   trans888_pub = 0.48599,
  local888 = A$auc_local_lipid_core,       local888_pub = 0.51503,
  full888  = A$auc_full,                   full888_pub  = 0.73306,
  full_w   = A$auc_full_w,                 full_w_pub   = 0.74828)
cat("\n==== baseline reproduction check (variant A vs published) ====\n")
keys <- names(chk); k <- 1
ok_all <- TRUE
while (k <= length(keys)) {
  nm <- keys[k]; got <- chk[[k]]; want <- chk[[k + 1]]
  d1 <- abs(got - want)
  flag <- if (d1 < 0.002) "PASS" else "FAIL"
  if (flag == "FAIL") ok_all <- FALSE
  cat(sprintf("  %-10s estimated %.4f | published %.4f | %s\n", nm, got, want, flag))
  k <- k + 2
}
cat(sprintf("\nBASELINE %s\n", if (ok_all) "REPRODUCED" else "NOT REPRODUCED - read variants with caution"))

write.csv(res, file.path(OUT, "phase20_outcome_ungradable_bounds.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "phase20_outcome_ungradable_bounds.csv"), "\n")
cat("ALL DONE\n")
