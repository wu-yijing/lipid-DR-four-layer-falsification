# 7f: 行级机制验证——A/B 规格下"插补值×any_dr"配对结构差异（解释 AUC 0.549 vs 0.623）
suppressMessages(library(mice))
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out"
sink(file.path(OUT, "phase7f_mechanism.log"), split = TRUE)

dev <- read.csv(file.path(P2, "dev_set.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b

dA <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
dA$tg_log <- log2(dA$tg); dA$tg <- NULL
dA$age10 <- dA$age / 10
dB <- dA[, c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]

impA <- mice(dA, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
impB <- mice(dB, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)

vars <- c("tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
mis <- sapply(dA[vars], is.na)

tab <- do.call(rbind, lapply(vars, function(v) {
  # 平均(跨 m)插补值，仅缺失行
  imp_matA <- sapply(1:20, function(i) complete(impA, i)[[v]][mis[, v]])
  imp_matB <- sapply(1:20, function(i) complete(impB, i)[[v]][mis[, v]])
  mA <- rowMeans(imp_matA); mB <- rowMeans(imp_matB)
  y <- dA$any_dr[mis[, v]]
  ttA <- try(t.test(mA ~ y), silent = TRUE); ttB <- try(t.test(mB ~ y), silent = TRUE)
  data.frame(var = v,
             n_imp = sum(mis[, v]),
             n_imp_case = sum(y),
             mean_case_A = mean(mA[y == 1]), mean_ctrl_A = mean(mA[y == 0]),
             diff_A = mean(mA[y == 1]) - mean(mA[y == 0]),
             p_A = if (inherits(ttA, "try-error")) NA else ttA$p.value,
             mean_case_B = mean(mB[y == 1]), mean_ctrl_B = mean(mB[y == 0]),
             diff_B = mean(mB[y == 1]) - mean(mB[y == 0]),
             p_B = if (inherits(ttB, "try-error")) NA else ttB$p.value,
             rowcorr_AB = suppressWarnings(cor(mA, mB)))
}))
cat("\n==== imputed rows: lipid by outcome, A vs B spec ====\n"); print(tab, digits = 3)
write.csv(tab, file.path(OUT, "mech_imp_by_outcome.csv"), row.names = FALSE)

# 拟合系数对比（pooled）
poolA <- pool(with(impA, glm(fml, family = binomial())))
poolB <- pool(with(impB, glm(fml, family = binomial())))
coef_tab <- data.frame(
  term = rownames(summary(poolA)),
  OR_A = exp(summary(poolA)$estimate), p_A = summary(poolA)$p.value,
  OR_B = exp(summary(poolB)$estimate), p_B = summary(poolB)$p.value)
cat("\n==== pooled coefficients (OR), A vs B ====\n"); print(coef_tab, digits = 3)
write.csv(coef_tab, file.path(OUT, "mech_pooled_coef.csv"), row.names = FALSE)

cat("\n==== ALL DONE (7f) ====\n")
sink()
