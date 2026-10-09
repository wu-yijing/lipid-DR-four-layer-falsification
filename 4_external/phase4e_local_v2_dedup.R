## AUTO-GENERATED from phase4e_local_v2.R -- sensitivity analysis: one record per patient
## (most recent screening year), input dev_set_v2_dedup.csv.

## AUTO-GENERATED from phase4e_local.R by make_v2_scripts.py
## Reads dev_set_v2_dedup.csv (corrected DR outcome, option A) and writes to *_out_v2/.
## Method unchanged from the original script.

## =====================================================================
## P1: fair three-row comparison in NHANES (v6 items M3, M9)
##   (i)   transported Chinese lipid score  (previous comparator)
##   (ii)  locally fitted lipid-only model  (NEW; same folds, same subset)
##   (iii) full reference model             (benchmark)
##   plus: recomputed external-validation AUC CIs (percentile + BCa)
## Pure ASCII.  Run: Rscript phase4e_local.R
## =====================================================================
suppressMessages({ library(survey); library(boot) })
set.seed(20260912)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2d"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2_dedup.csv"
NH  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv"
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

## ---------------- load + convert ----------------
nh <- read.csv(NH)
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]
cat(sprintf("NHANES diabetes subgroup n=%d events=%d (%.1f%%)\n",
            nrow(dm), sum(dm$any_dr), 100 * mean(dm$any_dr)))

## Chinese coefficients for the transported score (core version)
dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
cf <- coef(m_core)
cat(sprintf("Chinese core model fitted on n=%d (events=%d)\n", nrow(d), sum(d$any_dr)))

## ---------------- reference subset ----------------
ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR / 2
cat(sprintf("reference subset n=%d events=%d (%.1f%%); lipids available: TC+HDL=%d, +TG+LDL=%d\n",
            nrow(dref), sum(dref$any_dr), 100 * mean(dref$any_dr),
            sum(complete.cases(dref[, c("tc_mmol","hdl_mmol")])),
            sum(complete.cases(dref[, c("tc_mmol","hdl_mmol","tg_log","ldl_mmol")]))))

## ---------------- common 10-fold split ----------------
fold <- sample(rep(1:10, length.out = nrow(dref)))
## 2026-10-09 FIX (alignment defect, see tools/audit_repo_consistency.py check F2).
## The previous body was
##     p[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% as.numeric(coef(f))))
## which is wrong whenever a model's predictors contain missing values: model.matrix()
## drops the incomplete rows, the replacement vector is then shorter than p[fold == k],
## and R recycles it with a warning instead of failing. Only fml_lip_std is affected
## here (tg_log and ldl_mmol are missing in 473 of the 888 reference records), and the
## consequence was that the committed local_lipid_only_std row reported n = 888 and
## AUC 0.4861 for a quantity that exists on 415 records and is 0.4988. Predictions are
## now written back to their own row positions and the length is asserted, so the
## failure mode cannot recur silently. Re-running this script therefore CHANGES the
## local_lipid_only_std row of phase4e_three_row.csv; the other rows are unchanged.
oof <- function(fmla, dat, wcol = NULL) {
  p <- rep(NA_real_, nrow(dat))
  for (k in 1:10) {
    tr <- dat[fold != k, ]; te <- dat[fold == k, ]
    f <- try(if (is.null(wcol)) glm(fmla, data = tr, family = binomial())
             else survey::svyglm(fmla, design = survey::svydesign(ids = ~1, weights = as.formula(paste0("~", wcol)), data = tr), family = binomial()), silent = TRUE)
    if (inherits(f, "try-error")) next
    okk <- complete.cases(model.frame(fmla, data = te, na.action = na.pass))
    if (!any(okk)) next
    pr <- as.numeric(predict(f, newdata = te[okk, , drop = FALSE], type = "response"))
    stopifnot(length(pr) == sum(okk))
    p[which(fold == k)[okk]] <- pr
  }
  p
}
fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_lip_core <- any_dr ~ age10 + male + tc_mmol + hdl_mmol
fml_lip_std  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol + tg_log + ldl_mmol

p_full <- oof(fml_full, dref, "wt4")
p_loc_core <- oof(fml_lip_core, dref, "wt4")
p_loc_std  <- oof(fml_lip_std, dref, "wt4")
p_tr_core <- plogis(cf[1] + cf["age10"] * dref$age10 + cf["male"] * dref$male +
                      cf["tc"] * dref$tc_mmol + cf["hdl_c"] * dref$hdl_mmol)

ok_core <- !is.na(p_full) & !is.na(p_loc_core) & !is.na(p_tr_core)
cat(sprintf("\ncommon evaluable subset (full + local lipid + transported): n=%d events=%d\n",
            sum(ok_core), sum(dref$any_dr[ok_core])))
a_full <- auc_fun(p_full[ok_core], dref$any_dr[ok_core])
a_loc  <- auc_fun(p_loc_core[ok_core], dref$any_dr[ok_core])
a_tr   <- auc_fun(p_tr_core[ok_core], dref$any_dr[ok_core])
w_full <- auc_w_fun(p_full[ok_core], dref$any_dr[ok_core], dref$wt4[ok_core])
w_loc  <- auc_w_fun(p_loc_core[ok_core], dref$any_dr[ok_core], dref$wt4[ok_core])
w_tr   <- auc_w_fun(p_tr_core[ok_core], dref$any_dr[ok_core], dref$wt4[ok_core])
cat("\n==== THREE-ROW COMPARISON (unweighted | MEC-weighted) ====\n")
cat(sprintf("  (i)   transported Chinese lipid score : AUC %.3f | %.3f\n", a_tr, w_tr))
cat(sprintf("  (ii)  locally fitted lipid-only model : AUC %.3f | %.3f\n", a_loc, w_loc))
cat(sprintf("  (iii) full reference model            : AUC %.3f | %.3f\n", a_full, w_full))
cat(sprintf("  dAUC (iii)-(ii) = %.3f   dAUC (ii)-(i) = %.3f   dAUC (iii)-(i) = %.3f\n",
            a_full - a_loc, a_loc - a_tr, a_full - a_tr))

ok_std <- !is.na(p_loc_std)
a_loc_std <- auc_fun(p_loc_std[ok_std], dref$any_dr[ok_std])
w_loc_std <- auc_w_fun(p_loc_std[ok_std], dref$any_dr[ok_std], dref$wt4[ok_std])
cat(sprintf("  (ii-b) locally fitted lipid-only (std set, n=%d): AUC %.3f | %.3f\n",
            sum(ok_std), a_loc_std, w_loc_std))

write.csv(data.frame(row = c("transported_Chinese_lipid", "local_lipid_only_core",
                             "local_lipid_only_std", "full_reference_model"),
                     n = c(sum(ok_core), sum(ok_core), sum(ok_std), sum(ok_core)),
                     AUC_unweighted = c(a_tr, a_loc, a_loc_std, a_full),
                     AUC_MEC_weighted = c(w_tr, w_loc, w_loc_std, w_full)),
          file.path(OUT, "phase4e_three_row.csv"), row.names = FALSE)

## ---------------- bootstrap CIs for the transported versions ----------------
cat("\n==== recomputed AUC CIs for the transported versions (B=2000) ====\n")
set.seed(20260912)
res <- data.frame()
for (nm in c("core", "std")) {
  if (nm == "core") { lp <- cf[1] + cf["age10"]*dm$age10 + cf["male"]*dm$male +
                              cf["tc"]*dm$tc_mmol + cf["hdl_c"]*dm$hdl_mmol }
  else { m_std <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
         cs <- coef(m_std)
         lp <- cs[1] + cs["age10"]*dm$age10 + cs["male"]*dm$male + cs["tc"]*dm$tc_mmol +
               cs["hdl_c"]*dm$hdl_mmol + cs["tg_log"]*dm$tg_log + cs["ldl_c"]*dm$ldl_mmol }
  p <- plogis(lp); y <- dm$any_dr
  ok <- !is.na(p) & !is.na(y)
  pp <- p[ok]; yy <- y[ok]; n <- length(yy)
  a <- auc_fun(pp, yy)
  B <- 2000
  bs <- replicate(B, { i <- sample(n, n, TRUE); auc_fun(pp[i], yy[i]) })
  pct <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  bca <- try(boot.ci(boot(data = data.frame(p = pp, y = yy), statistic = function(dd, i) auc_fun(dd$p[i], dd$y[i]), R = B, stype = "i"), type = "bca"), silent = TRUE)
  bc <- if (inherits(bca, "try-error")) c(NA, NA) else bca$bca[4:5]
  cat(sprintf("  %-5s n=%4d ev=%3d AUC=%.4f | percentile %.4f-%.4f | BCa %s\n",
              nm, n, sum(yy), a, pct[1], pct[2], if (is.na(bc[1])) "n/a" else sprintf("%.4f-%.4f", bc[1], bc[2])))
  res <- rbind(res, data.frame(model = nm, n = n, events = sum(yy), auc = a,
                               pct_lo = pct[1], pct_hi = pct[2], bca_lo = bc[1], bca_hi = bc[2]))
}
write.csv(res, file.path(OUT, "phase4e_auc_ci.csv"), row.names = FALSE)

cat("\n==== ALL DONE ====\n")
