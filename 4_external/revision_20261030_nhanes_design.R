# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Reproduces the analyses added in response to Reviewer 1.
# Inputs are the extracts written by tools/revision_20261030_extract_union.sh; the
# GLGC 2013 summary statistics are fetched by tools/revision_20261030_download_glgc.py.
## =====================================================================
## R1 major points 3 and 4 -- recomputation (v2)
##   #3 ungradable fundus photographs: missing outcome, not a non-case
##   #4 NHANES survey design (strata + PSU) and fasting-subsample weights
## Outputs: E:/workbuddy/_r1_review_20261029/out/
## =====================================================================
suppressMessages({ library(survey); library(foreign) })
setwd("E:/workbuddy/_r1_review_20261029")
dir.create("out", showWarnings = FALSE)
options(stringsAsFactors = FALSE)

B   <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv")
DEV <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv")
B$cycle <- B[[grep("cycle", names(B), ignore.case = TRUE, value = TRUE)[1]]]

## ---------- design variables ----------
dem_D <- read.xport("nhanes/DEMO_D.XPT"); dem_E <- read.xport("nhanes/DEMO_E.XPT")
tri_D <- read.xport("nhanes/TRIGLY_D.XPT"); tri_E <- read.xport("nhanes/TRIGLY_E.XPT")
dem <- rbind(data.frame(seqn = dem_D$SEQN, SDMVSTRA = dem_D$SDMVSTRA, SDMVPSU = dem_D$SDMVPSU,
                        WTMEC2YR_d = dem_D$WTMEC2YR),
             data.frame(seqn = dem_E$SEQN, SDMVSTRA = dem_E$SDMVSTRA, SDMVPSU = dem_E$SDMVPSU,
                        WTMEC2YR_d = dem_E$WTMEC2YR))
tri <- rbind(data.frame(seqn = tri_D$SEQN, WTSAF2YR = tri_D$WTSAF2YR, LBXAPB_x = tri_D$LBXAPB),
             data.frame(seqn = tri_E$SEQN, WTSAF2YR = tri_E$WTSAF2YR, LBXAPB_x = NA_real_))
B <- merge(B, dem, by = "seqn", all.x = TRUE)
B <- merge(B, tri, by = "seqn", all.x = TRUE)
cat("weight agreement rows:", sum(!is.na(B$WTMEC2YR_d)), " max diff:",
    max(abs(B$WTMEC2YR - B$WTMEC2YR_d), na.rm = TRUE), "\n")
B$WTMEC2YR <- B$WTMEC2YR_d
B$WTMEC4YR <- B$WTMEC2YR/2
B$WTSAF4YR <- B$WTSAF2YR/2
B$gradable <- !is.na(B$dr_level4)
cat("ApoB availability by cycle:\n"); print(table(B$cycle, !is.na(B$LBXAPB_x)))

## ---------- analysis frame ----------
nh <- B
nh$tc_mmol <- nh$tc/38.67; nh$hdl_mmol <- nh$hdl/38.67
nh$age10 <- nh$age/10
dm <- nh[nh$diabetes == 1, ]
ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
cat("\nsubgroup n =", nrow(dm), " frame n =", nrow(dref),
    " gradable", sum(dref$gradable), " ungradable", sum(!dref$gradable), "\n")

## ---------- #3 ungradable vs gradable ----------
dref$grade <- factor(ifelse(dref$gradable, "gradable", "ungradable"),
                     levels = c("gradable", "ungradable"))
des0 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = dref, nest = TRUE)
vars <- c("age","male","black","bmi","hba1c_pct","dm_duration","insulin_use","oral_med",
          "sys_bp","egfr_ckdepi2021","smoker_100","tc","hdl","tg")
grab <- function(m, lev, col) if (lev %in% as.character(m$grade)) m[m$grade == lev, col] else NA_real_
char <- do.call(rbind, lapply(vars, function(v) {
  m <- svyby(as.formula(paste0("~", v)), ~grade, des0, svymean, na.rm = TRUE)
  sg <- try(sqrt(svyvar(as.formula(paste0("~", v)), subset(des0, grade == "gradable"), na.rm = TRUE))[1], silent = TRUE)
  su <- try(sqrt(svyvar(as.formula(paste0("~", v)), subset(des0, grade == "ungradable"), na.rm = TRUE))[1], silent = TRUE)
  data.frame(variable = v, gradable = grab(m,"gradable",v), ungradable = grab(m,"ungradable",v),
             sd_g = if (inherits(sg,"try-error")) NA else sg,
             sd_u = if (inherits(su,"try-error")) NA else su, row.names = NULL)
}))
char$diff <- char$ungradable - char$gradable
char$smd  <- char$diff/sqrt((char$sd_g^2 + char$sd_u^2)/2)
write.csv(char, "out/t3_ungradable_characteristics.csv", row.names = FALSE)
cat("\n=== ungradable vs gradable, weighted ===\n")
print(char[, c("variable","gradable","ungradable","smd")], digits = 3)

## ---------- predictions ----------
set.seed(20260912)
dref$wt4 <- dref$WTMEC4YR
fold_all <- sample(rep(1:10, length.out = nrow(dref)))
oof_pred <- function(fmla, dat, subset, fold) {
  p <- rep(NA_real_, nrow(dat)); idx <- which(subset); f2 <- fold[idx]
  for (k in 1:10) {
    tr <- idx[f2 != k]; te <- idx[f2 == k]
    f <- try(svyglm(fmla, design = svydesign(ids = ~1, weights = ~wt4, data = dat[tr, ]),
                    family = binomial()), silent = TRUE)
    if (inherits(f, "try-error")) next
    p[te] <- plogis(as.numeric(model.matrix(fmla, data = dat[te, ]) %*% as.numeric(coef(f))))
  }
  p
}
fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_lip  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol
g <- dref$gradable
## (A) reproduced submitted pipeline: all 888, ungradable labelled non-case
pA_full <- oof_pred(fml_full, dref, rep(TRUE, nrow(dref)), fold_all)
pA_loc  <- oof_pred(fml_lip,  dref, rep(TRUE, nrow(dref)), fold_all)
## (B) primary: gradable only
pB_full <- oof_pred(fml_full, dref, g, fold_all)
pB_loc  <- oof_pred(fml_lip,  dref, g, fold_all)
## transported Chinese score (fixed; independent of fitting)
d <- DEV[DEV$gradeable == 1 & complete.cases(DEV[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$age10 <- d$age/10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial()); cf <- coef(m_core)
p_tr <- plogis(cf[1] + cf["age10"]*dref$age10 + cf["male"]*dref$male +
               cf["tc"]*dref$tc_mmol + cf["hdl_c"]*dref$hdl_mmol)

auc_w <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]; w <- w[ok]
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  n1 <- sum(w*y); n0 <- sum(w*(1-y))
  num <- sum(w*y*(cumsum(w*(1-y)) - w*(1-y)))
  if (anyDuplicated(p)) for (ix in split(seq_along(p), p)) if (length(ix) > 1)
    num <- num + 0.5*(sum(w[ix]*y[ix])*sum(w[ix]*(1-y[ix])))
  num/(n1*n0)
}

## ---------- #3 scenario table ----------
lb_none <- dref$any_dr; lb_none[!g] <- 0L
lb_case <- dref$any_dr; lb_case[!g] <- 1L
tab3 <- rbind(
  data.frame(scenario = "A. as submitted (ungradable = non-case, model fitted under same convention)",
             n = nrow(dref), events = sum(lb_none),
             AUC_full_w = auc_w(pA_full, lb_none, dref$wt4), AUC_lipid_w = auc_w(pA_loc, lb_none, dref$wt4),
             AUC_tr_w = auc_w(p_tr, lb_none, dref$wt4),
             AUC_full_u = auc_w(pA_full, lb_none, rep(1,nrow(dref))), AUC_lipid_u = auc_w(pA_loc, lb_none, rep(1,nrow(dref))),
             AUC_tr_u = auc_w(p_tr, lb_none, rep(1,nrow(dref)))),
  data.frame(scenario = "B. primary: ascertainable retinal status only (model refitted on this subset)",
             n = sum(g), events = sum(lb_none[g]),
             AUC_full_w = auc_w(pB_full[g], lb_none[g], dref$wt4[g]), AUC_lipid_w = auc_w(pB_loc[g], lb_none[g], dref$wt4[g]),
             AUC_tr_w = auc_w(p_tr[g], lb_none[g], dref$wt4[g]),
             AUC_full_u = auc_w(pB_full[g], lb_none[g], rep(1,sum(g))), AUC_lipid_u = auc_w(pB_loc[g], lb_none[g], rep(1,sum(g))),
             AUC_tr_u = auc_w(p_tr[g], lb_none[g], rep(1,sum(g)))),
  data.frame(scenario = "C. worst case: all ungradable counted as cases (predictions from A)",
             n = nrow(dref), events = sum(lb_case),
             AUC_full_w = auc_w(pA_full, lb_case, dref$wt4), AUC_lipid_w = auc_w(pA_loc, lb_case, dref$wt4),
             AUC_tr_w = auc_w(p_tr, lb_case, dref$wt4),
             AUC_full_u = auc_w(pA_full, lb_case, rep(1,nrow(dref))), AUC_lipid_u = auc_w(pA_loc, lb_case, rep(1,nrow(dref))),
             AUC_tr_u = auc_w(p_tr, lb_case, rep(1,nrow(dref))))
)
tab3$dAUC_full_minus_lipid_w <- tab3$AUC_full_w - tab3$AUC_lipid_w
write.csv(tab3, "out/t3_scenarios.csv", row.names = FALSE)
cat("\n=== #3 label-convention scenarios ===\n"); print(tab3, digits = 4)

## ---------- #4 design-based AUC CI: survey replicate weights ----------
options(survey.lonely.psu = "adjust")
dref$p_full <- pA_full; dref$p_loc <- pA_loc; dref$p_tr <- p_tr
desA <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = dref, nest = TRUE)
repA <- as.svrepdesign(desA, type = "subbootstrap", replicates = 500)
rwA  <- weights(repA, type = "analysis")
cat("\n=== #4 design-consistent AUC (strata + PSU, subbootstrap replicate weights) ===\n")
tab4 <- do.call(rbind, lapply(c("p_full","p_loc","p_tr"), function(v) {
  ests <- apply(rwA, 2, function(w) auc_w(dref[[v]], dref$any_dr, w))
  data.frame(model = v, auc = auc_w(dref[[v]], dref$any_dr, dref$WTMEC4YR),
             se_design = sd(ests), lo_design = quantile(ests, .025, names = FALSE),
             hi_design = quantile(ests, .975, names = FALSE), row.names = NULL)
}))
tab4$width_design <- tab4$hi_design - tab4$lo_design
boot_srs <- function(p, y, w, Bn = 500) {
  out <- numeric(Bn); n <- length(y)
  for (b in 1:Bn) { i <- sample(n, n, TRUE); out[b] <- auc_w(p[i], y[i], w[i]) }
  c(lo = quantile(out, 0.025, names = FALSE), hi = quantile(out, 0.975, names = FALSE))
}
tab4$lo_srs <- NA; tab4$hi_srs <- NA
for (i in seq_len(nrow(tab4))) {
  s <- boot_srs(dref[[tab4$model[i]]], dref$any_dr, dref$WTMEC4YR)
  tab4$lo_srs[i] <- s[1]; tab4$hi_srs[i] <- s[2]
}
tab4$width_srs <- tab4$hi_srs - tab4$lo_srs
tab4$width_ratio <- tab4$width_design/tab4$width_srs
print(tab4, digits = 4); write.csv(tab4, "out/t4_design_auc.csv", row.names = FALSE)

## ---------- #3 primary analysis: design-based CIs and delta-AUC ----------
cat("\n=== #3 primary (gradable only): AUC and delta-AUC with design-based CIs ===\n")
dB <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = dref[g, ], nest = TRUE)
repB <- as.svrepdesign(dB, type = "subbootstrap", replicates = 500)
rwB <- weights(repB, type = "analysis"); yB <- dref$any_dr[g]; wB <- dref$WTMEC4YR[g]
yb <- lb_none[g]
b3 <- do.call(rbind, lapply(list(list("full", as.numeric(pB_full[g])),
                                 list("lipid", as.numeric(pB_loc[g])),
                                 list("transported", as.numeric(p_tr[g]))), function(z) {
  p <- as.numeric(z[[2]])
  ests <- apply(rwB, 2, function(w) auc_w(p, yb, w))
  data.frame(model = z[[1]], auc = auc_w(p, yb, wB), lo = quantile(ests, .025, names = FALSE),
             hi = quantile(ests, .975, names = FALSE), row.names = NULL)
}))
d_est <- auc_w(as.numeric(pB_full[g]), yb, wB) - auc_w(as.numeric(pB_loc[g]), yb, wB)
d_boot <- apply(rwB, 2, function(w) auc_w(as.numeric(pB_full[g]), yb, w) - auc_w(as.numeric(pB_loc[g]), yb, w))
cat("delta-AUC (full - lipid) =", round(d_est, 4), " 95% CI",
    paste(round(quantile(d_boot, c(.025,.975)), 4), collapse = " to "), "\n")
print(b3, digits = 4); write.csv(b3, "out/t3_primary_auc_ci.csv", row.names = FALSE)
write.csv(data.frame(delta = d_est, lo = quantile(d_boot, .025, names = FALSE),
                     hi = quantile(d_boot, .975, names = FALSE), n = sum(g), events = sum(yb)),
          "out/t3_primary_delta.csv", row.names = FALSE)

## ---------- #4 clinical model ----------
cat("\n=== #4 clinical reference model: ids=~1 (submitted) vs design-consistent ===\n")
fml_clin <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
old  <- svyglm(fml_clin, design = svydesign(ids = ~1, weights = ~wt4, data = dref), family = quasibinomial())
newo <- svyglm(fml_clin, design = des0, family = quasibinomial())
ci_o <- confint(old); ci_n <- confint(newo)
cmp <- data.frame(term = names(coef(old)),
                  OR_old = exp(coef(old)), lo_old = exp(ci_o[,1]), hi_old = exp(ci_o[,2]),
                  OR_des = exp(coef(newo)), lo_des = exp(ci_n[,1]), hi_des = exp(ci_n[,2]), row.names = NULL)
cmp$width_ratio <- (log(cmp$hi_des) - log(cmp$lo_des))/(log(cmp$hi_old) - log(cmp$lo_old))
print(cmp, digits = 3); write.csv(cmp, "out/t4_clinical_model_design.csv", row.names = FALSE)

## ---------- #4 fasting-subsample weights ----------
cat("\n=== #4 MEC vs fasting-subsample weights, same subset ===\n")
fast_ok <- !is.na(dref$WTSAF4YR) & !is.na(dref$tg)
cat("fasting subsample rows in frame:", sum(fast_ok), "of", nrow(dref), "\n")
cmpw <- data.frame(
  model = c("lipid_core","transported"),
  auc_MEC  = c(auc_w(pA_loc[fast_ok], lb_none[fast_ok], dref$WTMEC4YR[fast_ok]),
               auc_w(p_tr[fast_ok],  lb_none[fast_ok], dref$WTMEC4YR[fast_ok])),
  auc_FAST = c(auc_w(pA_loc[fast_ok], lb_none[fast_ok], dref$WTSAF4YR[fast_ok]),
               auc_w(p_tr[fast_ok],  lb_none[fast_ok], dref$WTSAF4YR[fast_ok])),
  auc_unw  = c(auc_w(pA_loc[fast_ok], lb_none[fast_ok], rep(1, sum(fast_ok))),
               auc_w(p_tr[fast_ok],  lb_none[fast_ok], rep(1, sum(fast_ok)))), row.names = NULL)
print(cmpw, digits = 4); write.csv(cmpw, "out/t4_fasting_weight_sensitivity.csv", row.names = FALSE)
cat("\n==== ALL DONE ====\n")
