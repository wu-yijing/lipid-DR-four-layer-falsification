## AUTO-GENERATED from phase2_model_v2.R -- sensitivity analysis: one record per patient
## (most recent screening year), input dev_set_v2_dedup.csv.

suppressMessages({ library(rms); library(glmnet) })
## AUTO-GENERATED from phase2_model.R by make_v2_scripts.py
## Reads dev_set_v2_dedup.csv (corrected DR outcome, option A) and writes to *_out_v2/.
## Method unchanged from the original script.

# Phase 2 modeling: bedside lipid DR risk score (development cohort)
# Target: any_dr (gradeable, worst-eye grade>=1)
# R-only ASCII; reads E:/workbuddy/.../phase2_dev/dev_set_v2_dedup.csv

set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out_v2d"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("==== load & subset ====\n")
d <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2_dedup.csv")
vars <- c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b")
sub <- d[d$gradeable == 1 & complete.cases(d[, c(vars, "any_dr")]), ]
cat(sprintf("modeling n=%d  events(any_dr)=%d  (%.1f%%)\n",
            nrow(sub), sum(sub$any_dr), 100*mean(sub$any_dr)))

sub$tg_log <- log2(sub$tg)
sub$age10 <- sub$age / 10
MV <- c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")
dd <- datadist(sub[, MV]); options(datadist = "dd")

cat("\n==== 2.2a full logistic (lrm) ====\n")
f_full <- lrm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
              data = sub, x = TRUE, y = TRUE)
print(f_full)

# store coefficient table IMMEDIATELY (before validate side-effects)
sd_map <- sapply(vars, function(v) sd(sub[[v]], na.rm = TRUE))
fc <- coef(f_full)
se <- sqrt(diag(vcov(f_full)))
pval <- 2 * pnorm(-abs(fc / se))
full_coef <- data.frame(
  variable = names(fc),
  coef = as.numeric(fc),
  OR = exp(as.numeric(fc)),
  OR_perSD = NA_real_,
  p = as.numeric(pval),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(full_coef))) {
  vn <- full_coef$variable[i]
  if (vn %in% names(sd_map)) full_coef$OR_perSD[i] <- exp(fc[i] * sd_map[vn])
}

cat("\n==== 2.2b LASSO (cv.glmnet, alpha=1) ====\n")
xmat <- model.matrix(~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b, sub)[, -1]
y <- sub$any_dr
cv <- cv.glmnet(xmat, y, family = "binomial", alpha = 1, nfolds = 10,
                type.measure = "deviance", standardize = TRUE)
b <- coef(cv, s = "lambda.min")
lasso_coef <- data.frame(var = rownames(b), coef = as.numeric(b))
lasso_coef <- lasso_coef[lasso_coef$coef != 0 | rownames(b) == "(Intercept)", ]
print(lasso_coef)

cat("\n==== 2.3 discrimination & calibration ====\n")
# store predictions BEFORE validate (validate may alter lrm object internals)
pred <- predict(f_full, type = "fitted")
lp <- predict(f_full)
val <- validate(f_full, method = "boot", B = 1000)
print(val)
auc_app <- (val["Dxy", "index.orig"] + 1) / 2
auc_corr <- (val["Dxy", "index.corrected"] + 1) / 2
cat(sprintf("Apparent AUC=%.3f  Optimism-corrected AUC=%.3f\n", auc_app, auc_corr))

brier <- mean((pred - sub$any_dr)^2, na.rm = TRUE)
cat(sprintf("Brier score=%.4f\n", brier))

g10 <- cut2(pred, g = 10)
cal_dec <- data.frame(
  pred_mean = as.numeric(tapply(pred, g10, mean)),
  obs_rate = as.numeric(tapply(sub$any_dr, g10, mean)),
  n = as.numeric(tapply(sub$any_dr, g10, length)),
  stringsAsFactors = FALSE)
cal_fit <- lrm(sub$any_dr ~ lp)
cal_slope <- as.numeric(coef(cal_fit)[2])
cal_int <- as.numeric(coef(cal_fit)[1])
cat(sprintf("calibration slope=%.3f  intercept=%.3f (ideal 1 / 0)\n", cal_slope, cal_int))
print(cal_dec)

cat("\n==== 2.4 NRI / IDI (self-implemented) ====\n")
fit0 <- glm(any_dr ~ age10 + male, data = sub, family = binomial)
fit1 <- glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
            data = sub, family = binomial)
p0 <- predict(fit0, type = "response")
p1 <- predict(fit1, type = "response")
ev <- sub$any_dr == 1
idi <- (mean(p1[ev]) - mean(p0[ev])) - (mean(p1[!ev]) - mean(p0[!ev]))
nri_cat <- function(cut) {
  up_e <- mean(p1[ev] >= cut & p0[ev] < cut); down_e <- mean(p1[ev] < cut & p0[ev] >= cut)
  up_n <- mean(p1[!ev] >= cut & p0[!ev] < cut); down_n <- mean(p1[!ev] < cut & p0[!ev] >= cut)
  c(up_e - down_e + down_n - up_n, up_e, down_e, down_n, up_n)
}
cuts <- c(median(p1), 0.10, 0.20)
nri_res <- t(sapply(cuts, function(c) c(cut = c, nri_cat(c))))
print(nri_res)
cat(sprintf("IDI=%.4f\n", idi))

cat("\n==== 2.5 nomogram & score formula ====\n")
# re-fit clean model (isolate from validate side-effects)
f_nom <- lrm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
             data = sub, x = TRUE, y = TRUE)
nom <- nomogram(f_nom, fun = plogis, funlabel = "Predicted DR risk",
                lp.at = seq(-2, 2, 0.5))
png(file.path(OUT, "phase2_nomogram.png"), width = 900, height = 640, res = 120)
plot(nom)
dev.off()
cat("nomogram written\n")

cf <- coef(f_nom)
lp_formula <- paste0("lp = ", round(cf[1], 4))
for (i in 2:length(cf)) lp_formula <- paste0(lp_formula,
    sprintf(" + %.4f*%s", cf[i], names(cf)[i]))
risk_formula <- "risk = 1 / (1 + exp(-lp))"
cat(lp_formula, "\n")
cat(risk_formula, "\n")

# ---- exports ----
write.csv(full_coef, file.path(OUT, "phase2_model_fit.csv"), row.names = FALSE)
write.csv(lasso_coef, file.path(OUT, "phase2_lasso.csv"), row.names = FALSE)
write.csv(cal_dec, file.path(OUT, "phase2_calibration.csv"), row.names = FALSE)
perf <- data.frame(
  metric = c("AUC_apparent","AUC_optimism_corrected","Brier",
             "calib_slope","calib_intercept","IDI"),
  value = c(auc_app, auc_corr, brier, cal_slope, cal_int, idi))
write.csv(perf, file.path(OUT, "phase2_performance.csv"), row.names = FALSE)
nri_out <- as.data.frame(nri_res); nri_out$IDI <- idi
write.csv(nri_out, file.path(OUT, "phase2_nri_idi.csv"), row.names = FALSE)
writeLines(c(lp_formula, risk_formula, "", "Note: age10 = age/10; tg_log = log2(TG).",
             "Bedside use: compute lp, then risk = 1/(1+exp(-lp))."),
          file.path(OUT, "phase2_score_formula.txt"))

cat("\n==== ALL DONE ====\n")
