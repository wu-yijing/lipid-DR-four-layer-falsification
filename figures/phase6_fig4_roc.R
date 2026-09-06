# Figure 4: ROC comparison on the comparable NHANES subset (n=888)
# Lipid-only (Chinese coefficients) vs Full model (10-fold out-of-fold)
suppressMessages({ library(rms); library(survey) })
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase6_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "fig4_roc.log"), split = TRUE)

nh <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_ext/NHANES_B.csv")
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$age10 <- nh$age / 10
dm <- nh[nh$diabetes == 1, ]

ref_vars <- c("age10","male","hba1c_pct","dm_duration","egfr_ckdepi2021","bmi",
              "smoker_100","insulin_use","oral_med","sys_bp","any_dr","WTMEC2YR")
dref <- dm[complete.cases(dm[, ref_vars]), ]
fmla <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
        smoker_100 + insulin_use + oral_med + sys_bp
fit1 <- function(dat) survey::svyglm(fmla, design = survey::svydesign(ids = ~1,
        weights = ~WTMEC2YR, data = dat), family = binomial())

# lipid-only predictions (Chinese core coefficients)
dev <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv")
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
cf <- coef(glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial()))
p_lip <- plogis(cf[1] + cf["age10"]*dref$age10 + cf["male"]*dref$male +
                cf["tc"]*dref$tc_mmol + cf["hdl_c"]*dref$hdl_mmol)

# full model out-of-fold predictions (same 10-fold scheme as phase 4)
fold <- sample(rep(1:10, length.out = nrow(dref)))
p_ref <- rep(NA_real_, nrow(dref))
for (k in 1:10) {
  tr <- dref[fold != k, ]; te <- dref[fold == k, ]
  p_ref[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% coef(fit1(tr))))
}

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1*(n1+1)/2) / (n1*n0)
}
roc_pts <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  o <- order(p, decreasing = TRUE); y <- y[o]
  data.frame(fpr = c(0, cumsum(1 - y) / sum(1 - y)),
             tpr = c(0, cumsum(y) / sum(y)))
}
auc_lip <- auc_fun(p_lip, dref$any_dr)
auc_ref <- auc_fun(p_ref, dref$any_dr)
cat(sprintf("n=%d events=%d | lipid-only AUC=%.3f | full-model AUC=%.3f\n",
            nrow(dref), sum(dref$any_dr), auc_lip, auc_ref))

BLUE <- "#2C7FB8"; RED <- "#D7301F"; GREY <- "#636363"
png(file.path(OUT, "fig4_roc_comparison.png"), width = 1900, height = 1900, res = 300)
par(mar = c(4.6, 4.6, 3, 1.5), family = "sans")
r_ref <- roc_pts(p_ref, dref$any_dr)
r_lip <- roc_pts(p_lip, dref$any_dr)
plot(tpr ~ fpr, data = r_ref, type = "l", col = BLUE, lwd = 2.4,
     xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity (FPR)",
     ylab = "Sensitivity (TPR)", main = "NHANES diabetes subgroup (n = 888)",
     cex.main = 1.05, cex.lab = 0.95)
lines(tpr ~ fpr, data = r_lip, col = RED, lwd = 2.4)
abline(0, 1, lty = 2, col = GREY, lwd = 1.2)
legend("bottomright", bty = "n", cex = 0.92, lwd = 2.4, col = c(BLUE, RED),
       legend = c(sprintf("Full model: HbA1c + duration + insulin + ... (AUC = %.3f)", auc_ref),
                  sprintf("Lipid-only bedside score (AUC = %.3f)", auc_lip)))
dev.off()
cat("fig4_roc_comparison.png written\n")
write.csv(data.frame(fpr = r_ref$fpr, tpr_full = r_ref$tpr, tpr_lipid = r_lip$tpr),
          file.path(OUT, "fig4_roc_curves.csv"), row.names = FALSE)
sink()
