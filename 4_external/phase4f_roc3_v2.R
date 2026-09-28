## AUTO-GENERATED from phase4f_roc3.R by make_v2_scripts.py
## Reads dev_set_v2.csv (corrected DR outcome, option A) and writes to *_out_v2/.
## Method unchanged from the original script.

## =====================================================================
## Figure 4 redraw input: ROC points for the three-row comparison
##   full reference model | locally fitted lipid-only | transported Chinese score
## Uses the identical fold split and subset as phase4e_local.R
## =====================================================================
suppressMessages({ library(survey) })
set.seed(20260912)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv"
NH  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv"

nh <- read.csv(NH)
nh$tc_mmol <- nh$tc/38.67; nh$hdl_mmol <- nh$hdl/38.67; nh$tg_mmol <- nh$tg/88.57
nh$ldl_mmol <- nh$ldl/38.67; nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age/10
dm <- nh[nh$diabetes == 1, ]
dev <- read.csv(DEV)
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age/10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial()); cf <- coef(m_core)

ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]; dref$wt4 <- dref$WTMEC2YR/2
fold <- sample(rep(1:10, length.out = nrow(dref)))
oof <- function(fmla, dat) {
  p <- rep(NA_real_, nrow(dat))
  for (k in 1:10) {
    tr <- dat[fold != k, ]; te <- dat[fold == k, ]
    f <- try(survey::svyglm(fmla, design = survey::svydesign(ids=~1, weights=~wt4, data=tr), family=binomial()), silent=TRUE)
    if (inherits(f, "try-error")) next
    p[fold == k] <- plogis(as.numeric(model.matrix(fmla, data=te) %*% as.numeric(coef(f))))
  }
  p
}
fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_lip  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol
p_full <- oof(fml_full, dref); p_loc <- oof(fml_lip, dref)
p_tr <- plogis(cf[1] + cf["age10"]*dref$age10 + cf["male"]*dref$male + cf["tc"]*dref$tc_mmol + cf["hdl_c"]*dref$hdl_mmol)
ok <- !is.na(p_full) & !is.na(p_loc) & !is.na(p_tr)
y <- dref$any_dr[ok]
roc <- function(p, y) {
  o <- order(-p); yy <- y[o]
  tpr <- cumsum(yy)/sum(yy); fpr <- cumsum(1-yy)/sum(1-yy)
  data.frame(fpr = c(0, fpr), tpr = c(0, tpr))
}
r1 <- roc(p_full[ok], y); r2 <- roc(p_loc[ok], y); r3 <- roc(p_tr[ok], y)
n <- max(nrow(r1), nrow(r2), nrow(r3))
idx <- function(r) round(seq(1, nrow(r), length.out = n))
out <- data.frame(fpr = r1$fpr[idx(r1)], tpr_full = r1$tpr[idx(r1)],
                  tpr_local = r2$tpr[idx(r2)], tpr_chinese = r3$tpr[idx(r3)])
write.csv(out, file.path(OUT, "phase4f_roc3.csv"), row.names = FALSE)
cat(sprintf("ROC rows=%d  (n=%d, events=%d)\n", nrow(out), sum(ok), sum(y)))
cat("==== ALL DONE ====\n")
