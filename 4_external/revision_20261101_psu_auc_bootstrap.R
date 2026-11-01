# -*- coding: utf-8 -*-
## Design-based intervals for the three layer-4 AUCs (Reviewer 1, point 4, second pass)
##
## Adds what "record-level percentile bootstrap" cannot give: an interval that reflects
## the NHANES design (31 strata, 62 PSUs, four-year MEC weights).
##   (a) weighted AUC  -> PSU subbootstrap replicate weights (survey::as.svrepdesign)
##   (b) unweighted AUC -> PSU cluster bootstrap (PSUs resampled within strata)
##   (c) delta-AUC     -> same resamples, paired
## Point estimates reproduce layer4_three_row_AUC_primary.csv exactly
## (0.7856 / 0.5397 / 0.4674 unweighted; 0.8039 weighted for the full model).
## Writes layer4_auc_design_ci.csv and layer4_delta_auc_design_ci.csv to LDF_OUT
## (default: a session temp directory, so a run leaves the working tree clean).
suppressMessages({ library(survey) })
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
DEM <- Sys.getenv("LDF_NHANES_XPT", "E:/workbuddy/_r1_review_20261029/nhanes")
OUT <- Sys.getenv("LDF_OUT", file.path(tempdir(), "psu_auc"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B   <- 1000
set.seed(20260912)

nh <- read.csv(file.path(SRC, "phase2_ext/NHANES_B.csv"))
nh$tc_mmol <- nh$tc/38.67; nh$hdl_mmol <- nh$hdl/38.67; nh$age10 <- nh$age/10
dm <- nh[nh$diabetes == 1, ]
dm$any_dr[is.na(dm$dr_level4)] <- NA          # ungradable = missing outcome (Methods 2.2)
ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
dref$wt4 <- dref$WTMEC2YR/2

dem_D <- foreign::read.xport(file.path(DEM, "DEMO_D.XPT"))
dem_E <- foreign::read.xport(file.path(DEM, "DEMO_E.XPT"))
dem <- rbind(data.frame(seqn = dem_D$SEQN, SDMVSTRA = dem_D$SDMVSTRA, SDMVPSU = dem_D$SDMVPSU),
             data.frame(seqn = dem_E$SEQN, SDMVSTRA = dem_E$SDMVSTRA, SDMVPSU = dem_E$SDMVPSU))
dref <- merge(dref, dem, by = "seqn", all.x = TRUE)

fold <- sample(rep(1:10, length.out = nrow(dref)))
oof <- function(fmla, dat) {
  p <- rep(NA_real_, nrow(dat))
  for (k in 1:10) {
    tr <- dat[fold != k, ]; te <- dat[fold == k, ]
    f <- try(svyglm(fmla, design = svydesign(ids = ~1, weights = ~wt4, data = tr),
                    family = binomial()), silent = TRUE)
    if (inherits(f, "try-error")) next
    p[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% as.numeric(coef(f))))
  }
  p
}
fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_lip  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol
p_full <- oof(fml_full, dref); p_loc <- oof(fml_lip, dref)

dev <- read.csv(file.path(SRC, "phase2_dev/dev_set_v2.csv"))
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
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

y <- dref$any_dr; w <- dref$WTMEC2YR/2
P <- list(full = p_full, lipid = p_loc, transported = p_tr)

## (a) weighted AUC, design-based
dB <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt4, data = dref, nest = TRUE)
set.seed(20261101)
rw <- weights(as.svrepdesign(dB, type = "subbootstrap", replicates = B), type = "analysis")
W <- lapply(P, function(p) apply(rw, 2, function(ww) auc_w(p, y, ww)))

## (b) unweighted AUC, PSU cluster bootstrap
idx <- split(seq_len(nrow(dref)), interaction(dref$SDMVSTRA, dref$SDMVPSU, drop = TRUE))
strata <- split(seq_along(idx), sapply(idx, function(ix) dref$SDMVSTRA[ix[1]]))
set.seed(20261101)
U <- lapply(P, function(p) numeric(B))
for (b in 1:B) {
  pick <- unlist(lapply(strata, function(ss) sample(ss, length(ss), replace = TRUE)), use.names = FALSE)
  rows <- unlist(idx[pick], use.names = FALSE)
  for (m in names(P)) U[[m]][b] <- auc_w(P[[m]][rows], y[rows], rep(1, length(rows)))
}

q <- function(x) quantile(x, c(.025, .975), names = FALSE)
out <- do.call(rbind, lapply(names(P), function(m) data.frame(
  model = m, n = length(y), events = sum(y),
  auc_unweighted = auc_w(P[[m]], y, rep(1, length(y))),
  lo_psu = q(U[[m]])[1], hi_psu = q(U[[m]])[2],
  auc_weighted = auc_w(P[[m]], y, w),
  lo_design = q(W[[m]])[1], hi_design = q(W[[m]])[2])))
d1 <- U$full - U$lipid; d2 <- U$full - U$transported; d3 <- U$lipid - U$transported
delta <- data.frame(
  contrast = c("full-lipid","full-transported","lipid-transported"),
  delta = c(auc_w(P$full,y,rep(1,length(y))) - auc_w(P$lipid,y,rep(1,length(y))),
            auc_w(P$full,y,rep(1,length(y))) - auc_w(P$transported,y,rep(1,length(y))),
            auc_w(P$lipid,y,rep(1,length(y))) - auc_w(P$transported,y,rep(1,length(y)))),
  lo_psu = c(q(d1)[1], q(d2)[1], q(d3)[1]),
  hi_psu = c(q(d1)[2], q(d2)[2], q(d3)[2]))
write.csv(out, file.path(OUT, "layer4_auc_design_ci.csv"), row.names = FALSE)
write.csv(delta, file.path(OUT, "layer4_delta_auc_design_ci.csv"), row.names = FALSE)
print(out, digits = 4); print(delta, digits = 4)
