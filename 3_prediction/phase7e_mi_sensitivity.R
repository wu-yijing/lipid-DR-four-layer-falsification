# 7e: MI 规格敏感性诊断——为什么 A-spec(保留age)=0.549 而 B-spec(删age)=0.623？
suppressMessages(library(mice))
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out"
sink(file.path(OUT, "phase7e_mi_sensitivity.log"), split = TRUE)

dev <- read.csv(file.path(P2, "dev_set.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]
auc_fun <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b

mkA <- function() {
  d <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
  d$tg_log <- log2(d$tg); d$tg <- NULL
  d$age10 <- d$age / 10
  d
}
mkB <- function() {
  d <- mkA()
  d[, c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
}

## ---- 1) predictorMatrix & logged events (seed 20260905, m=2 for speed) ----
iA <- mice(mkA(), m = 2, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
iB <- mice(mkB(), m = 2, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
cat("==== A loggedEvents ====\n"); print(iA$loggedEvents)
cat("==== B loggedEvents ====\n"); print(iB$loggedEvents)
cat("==== A predictors of tg_log ====\n"); print(iA$predictorMatrix["tg_log", ])
cat("==== B predictors of tg_log ====\n"); print(iB$predictorMatrix["tg_log", ])

## ---- 2) imputed-value distributions A vs B (from the m=2 runs) ----
cmp <- do.call(rbind, lapply(c("tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b"), function(v) {
  mis <- is.na(mkA()[[v]])
  obs <- mkA()[[v]][!mis]
  gA <- unlist(lapply(1:2, function(i) complete(iA, i)[[v]][mis]))
  gB <- unlist(lapply(1:2, function(i) complete(iB, i)[[v]][mis]))
  data.frame(var = v, mean_obs = mean(obs), mean_impA = mean(gA), mean_impB = mean(gB),
             sd_obs = sd(obs), sd_impA = sd(gA), sd_impB = sd(gB))
}))
cat("\n==== imputed vs observed, A vs B ====\n"); print(cmp)

## ---- 3) seed sensitivity: 3 seeds x m=5 per spec ----
run_spec <- function(d, seed, m = 5) {
  imp <- mice(d, m = m, method = "pmm", printFlag = FALSE, seed = seed, maxit = 10)
  mean(sapply(1:m, function(i) {
    di <- complete(imp, i)
    fi <- glm(fml, data = di, family = binomial())
    auc_fun(predict(fi, type = "response"), di$any_dr)
  }))
}
seeds <- c(11, 22, 33)
resA <- sapply(seeds, function(s) run_spec(mkA(), s))
resB <- sapply(seeds, function(s) run_spec(mkB(), s))
cat("\n==== seed sensitivity (mean MI AUC) ====\n")
cat(sprintf("A-spec (age retained): seeds 11/22/33 -> %s | mean=%.4f\n",
            paste(sprintf("%.4f", resA), collapse = ", "), mean(resA)))
cat(sprintf("B-spec (age dropped):  seeds 11/22/33 -> %s | mean=%.4f\n",
            paste(sprintf("%.4f", resB), collapse = ", "), mean(resB)))

## ---- 4) maxit sensitivity on the primary seed (convergence check) ----
for (mx in c(20, 40)) {
  aA <- run_spec(mkA(), 20260905)  # m=5 default in run_spec
  cat(sprintf("maxit=%d A-spec seed20260905 m=5: %.4f\n", mx, aA))
}
cat("\n==== ALL DONE (7e) ====\n")
sink()
