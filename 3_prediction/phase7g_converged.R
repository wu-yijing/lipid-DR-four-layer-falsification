# 7g: 收敛级 MI——双规格 × 3 种子，m=50, maxit=50（确定 MI AUC 的"收敛值域"）
suppressMessages(library(mice))
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out"
sink(file.path(OUT, "phase7g_converged.log"), split = TRUE)

dev <- read.csv(file.path(P2, "dev_set.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b
auc_fun <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

mkA <- function() {
  d <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
  d$tg_log <- log2(d$tg); d$tg <- NULL
  d$age10 <- d$age / 10
  d
}
mkB <- function() mkA()[, c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c",
                            "apo_a1", "apo_b", "any_dr")]

res <- data.frame()
t0 <- Sys.time()
for (spec in c("A", "B")) {
  d <- if (spec == "A") mkA() else mkB()
  for (sd_ in c(101, 202, 303)) {
    imp <- mice(d, m = 50, method = "pmm", printFlag = FALSE, seed = sd_, maxit = 50)
    aucs <- sapply(1:50, function(i) {
      di <- complete(imp, i)
      fi <- glm(fml, data = di, family = binomial())
      auc_fun(predict(fi, type = "response"), di$any_dr)
    })
    res <- rbind(res, data.frame(spec = spec, seed = sd_, m = 50, maxit = 50,
                                 auc_mean = mean(aucs), auc_sd = sd(aucs)))
    cat(sprintf("[%s seed=%d] mean=%.4f sd=%.4f  elapsed=%.1f min\n",
                spec, sd_, mean(aucs), sd(aucs),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
    write.csv(res, file.path(OUT, "phase7g_converged_partial.csv"), row.names = FALSE)
  }
}
cat("\n==== converged MI AUC summary ====\n")
print(res)
agg <- aggregate(auc_mean ~ spec, res, function(x) c(mean = mean(x), min = min(x), max = max(x)))
print(agg)
write.csv(res, file.path(OUT, "phase7g_converged.csv"), row.names = FALSE)
cat("\n==== ALL DONE (7g) ====\n")
sink()
