## AUTO-GENERATED v2 (corrected outcome, option A)

# Phase 7b: 缺失诊断三件套（R2-1）—— 缺失率表 / MAR 依据 / 插补诊断
suppressMessages({ library(mice); library(lattice) })
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out_v2"

dev <- read.csv(file.path(P2, "dev_set_v2.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]
d <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
d$tg_log <- log2(d$tg); d$tg <- NULL
d$age10 <- d$age / 10; d$age <- NULL
d <- d[, c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
cat(sprintf("n=%d (gradeable), events=%d\n", nrow(d), sum(d$any_dr, na.rm = TRUE)))

vars_model <- c("tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
vars_all <- c("age10", "male", vars_model)

## ---- 1) per-variable missing rate table ----
rates <- do.call(rbind, lapply(vars_all, function(v) {
  miss <- is.na(d[[v]])
  data.frame(var = v,
             n_missing = sum(miss),
             pct_missing = round(100 * mean(miss), 2),
             missing_in_cases = sum(miss & d$any_dr == 1, na.rm = TRUE),
             missing_in_controls = sum(miss & d$any_dr == 0, na.rm = TRUE))
}))
cat("\n==== per-variable missing rates ====\n"); print(rates)
write.csv(rates, file.path(OUT, "missing_rates.csv"), row.names = FALSE)

md <- mice::md.pattern(d, plot = FALSE)
cat("\n==== missing pattern (md.pattern) ====\n"); print(md)
write.csv(md, file.path(OUT, "missing_pattern.csv"))

## ---- 2) MAR assessment: missingness (>=1 lipid) vs observed characteristics ----
d$miss_any_lipid <- rowSums(is.na(d[, vars_model])) > 0
cat(sprintf("\nmiss >=1 lipid: n=%d (%.1f%%)\n",
            sum(d$miss_any_lipid), 100 * mean(d$miss_any_lipid)))
sub <- d[complete.cases(d[, c("age10", "male", "any_dr")]), ]
m1 <- glm(miss_any_lipid ~ age10 + male + any_dr, data = sub, family = binomial())
s1 <- summary(m1)$coefficients
or_tab <- data.frame(term = rownames(s1), OR = exp(s1[, 1]),
                     lci95 = exp(s1[, 1] - 1.96 * s1[, 2]),
                     uci95 = exp(s1[, 1] + 1.96 * s1[, 2]),
                     p = s1[, 4])
cat("\n==== logistic: P(missing>=1 lipid) ~ age10 + male + any_dr ====\n"); print(or_tab)
write.csv(or_tab, file.path(OUT, "missing_assoc_or.csv"), row.names = FALSE)

## univariate tests (observed-characteristic comparison)
tt <- t.test(age10 ~ miss_any_lipid, data = sub)
xt_m <- table(sub$male, sub$miss_any_lipid)
xt_d <- table(sub$any_dr, sub$miss_any_lipid)
uni <- data.frame(
  comparison = c("age10 (t-test)", "male (chi2)", "any_dr (chi2)"),
  statistic = c(tt$statistic, suppressWarnings(chisq.test(xt_m))$statistic,
                suppressWarnings(chisq.test(xt_d))$statistic),
  p = c(tt$p.value, suppressWarnings(chisq.test(xt_m))$p.value,
        suppressWarnings(chisq.test(xt_d))$p.value))
cat("\n==== univariate missingness comparisons ====\n"); print(uni)
write.csv(uni, file.path(OUT, "missing_assoc_univariate.csv"), row.names = FALSE)

## ---- 3) imputation diagnostics (same mice call as phase2b_mi.R) ----
imp <- mice(d, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
cat("\nimputation done\n")

png(file.path(OUT, "imp_densityplot.png"), width = 1800, height = 1400, res = 150)
print(densityplot(imp, ~ tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
                  layout = c(3, 2), xlab = ""))
dev.off()
cat("density plot saved\n")

diag_tab <- do.call(rbind, lapply(vars_model, function(v) {
  obs <- d[[v]][!is.na(d[[v]])]
  imp_means <- sapply(1:20, function(i) mean(complete(imp, i)[[v]][is.na(d[[v]])]))
  imp_sds   <- sapply(1:20, function(i) sd(complete(imp, i)[[v]][is.na(d[[v]])]))
  data.frame(var = v,
             n_obs = length(obs), mean_obs = mean(obs), sd_obs = sd(obs),
             mean_imputed = mean(imp_means), sd_imputed = mean(imp_sds),
             shift_pct = round(100 * (mean(imp_means) - mean(obs)) / sd(obs), 2))
}))
cat("\n==== observed vs imputed distribution comparison ====\n"); print(diag_tab)
write.csv(diag_tab, file.path(OUT, "imp_diag_stats.csv"), row.names = FALSE)

## fraction of missing information (FMI) from pooled fits
fits <- with(imp, glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c +
                        apo_a1 + apo_b, family = binomial()))
pooled <- pool(fits)
fmi_tab <- data.frame(term = names(coef(fits$analyses[[1]])),
                      fmi = pooled$pooled$fmi,
                      lambda = pooled$pooled$lambda)
cat("\n==== FMI / lambda (pooled logistic) ====\n"); print(fmi_tab)
write.csv(fmi_tab, file.path(OUT, "imp_fmi.csv"), row.names = FALSE)

cat("\n==== ALL DONE (7b) ====\n")
