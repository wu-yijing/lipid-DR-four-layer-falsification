# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Reproduces the analyses added in response to Reviewer 1.
# Inputs are the extracts written by tools/revision_20261030_extract_union.sh; the
# GLGC 2013 summary statistics are fetched by tools/revision_20261030_download_glgc.py.
## Redraw Figure 4 from CURRENT predictions; tie-collapsed step ROC; unweighted (matches reported AUCs)
suppressMessages({ library(survey) })
set.seed(20260912)
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv"
NH  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv"
OUT <- "E:/workbuddy/Front  in  Endocrinology  投稿/投稿资料/审稿回应_InteractiveReview_20261029"

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
ok <- !is.na(p_full) & !is.na(p_loc) & !is.na(p_tr); y <- dref$any_dr[ok]

## step ROC with tie collapsing
roc_step <- function(p, y, npts = 2000) {
  u <- sort(unique(p), decreasing = TRUE)
  fpr <- tpr <- numeric(length(u)+1); fpr[1] <- 0; tpr[1] <- 0
  for (i in seq_along(u)) {
    pred <- p >= u[i]
    tpr[i+1] <- sum(y[pred] == 1)/sum(y == 1)
    fpr[i+1] <- sum(y[pred] == 0)/sum(y == 0)
  }
  data.frame(fpr = fpr, tpr = tpr)
}
rk <- function(p, y) {
  r <- rank(p, ties.method = "average"); n1 <- sum(y==1); n0 <- sum(y==0)
  (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)
}
R1 <- roc_step(p_full[ok], y); R2 <- roc_step(p_loc[ok], y); R3 <- roc_step(p_tr[ok], y)
cat(sprintf("full %.4f | local %.4f | transported %.4f\n", rk(p_full[ok],y), rk(p_loc[ok],y), rk(p_tr[ok],y)))
cat(sprintf("max fpr overlap check: npred=%d, points full/local/trans = %d/%d/%d\n", sum(ok), nrow(R1), nrow(R2), nrow(R3)))

png(file.path(OUT, "Fig4_check_corrected.png"), width = 1100, height = 1050, res = 160)
par(mar = c(4.6, 4.8, 1.2, 1.2), mgp = c(2.9, 0.75, 0), cex.axis = 1.05, cex.lab = 1.15)
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "1 \u2212 specificity", ylab = "Sensitivity", axes = FALSE)
box(col = "grey35"); axis(1, at = seq(0,1,0.2), col = "grey35", col.axis = "black")
axis(2, at = seq(0,1,0.2), las = 1, col = "grey35", col.axis = "black")
abline(0, 1, lty = 2, col = "grey55", lwd = 1.3)
lines(R1$fpr, R1$tpr, type = "s", col = "#B4639B", lwd = 2.6)
lines(R2$fpr, R2$tpr, type = "s", col = "#E1962B", lwd = 2.6)
lines(R3$fpr, R3$tpr, type = "s", col = "#1F6FB2", lwd = 2.6, lty = 1)
legend("bottomright", bty = "n", cex = 0.92, lwd = 2.6,
       col = c("#B4639B", "#E1962B", "#1F6FB2"),
       legend = c(sprintf("Full reference model (AUC %.3f)", rk(p_full[ok],y)),
                  sprintf("Local lipid-only model (AUC %.3f)", rk(p_loc[ok],y)),
                  sprintf("Transported Chinese score (AUC %.3f)", rk(p_tr[ok],y))))
dev.off()
write.csv(data.frame(fpr = R1$fpr, tpr_full = R1$tpr)[1:min(nrow(R1),nrow(R2),nrow(R3)), ],
          file.path(OUT, "fig4_roc_recomputed_head.csv"), row.names = FALSE)
cat("WROTE png + csv\n")
