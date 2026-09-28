# 7d: 精确复现 phase2b_mi.R 的 MI AUC（同列结构/顺序/seed），并对比 v3 变体
suppressMessages(library(mice))
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
dev <- read.csv(file.path(P2, "dev_set.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]

auc_fun <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## --- A: phase2b 原样（保留 age，age10 追加在末尾；列序与原脚本一致） ---
dA <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
dA$tg_log <- log2(dA$tg); dA$tg <- NULL
dA$age10 <- dA$age / 10
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b
impA <- mice(dA, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
aucsA <- sapply(1:20, function(i) {
  di <- complete(impA, i)
  fi <- glm(fml, data = di, family = binomial())
  auc_fun(predict(fi, type = "response"), di$any_dr)
})
cat(sprintf("A (phase2b exact, age retained): mean=%.4f sd=%.4f range=[%.3f, %.3f]\n",
            mean(aucsA), sd(aucsA), min(aucsA), max(aucsA)))
print(round(aucsA, 3))

## --- B: v3 变体（删 age 只留 age10） ---
dB <- dA[, c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
impB <- mice(dB, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
aucsB <- sapply(1:20, function(i) {
  di <- complete(impB, i)
  fi <- glm(fml, data = di, family = binomial())
  auc_fun(predict(fi, type = "response"), di$any_dr)
})
cat(sprintf("B (v3 variant, age dropped):      mean=%.4f sd=%.4f range=[%.3f, %.3f]\n",
            mean(aucsB), sd(aucsB), min(aucsB), max(aucsB)))
print(round(aucsB, 3))
