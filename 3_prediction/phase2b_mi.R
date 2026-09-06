# P2b: Multiple imputation sensitivity analysis for the Chinese cohort model
# Compares listwise-deletion primary model (n=1,207) with MI (m=20, gradeable n=1,966)
suppressMessages(library(mice))
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
sink(file.path(OUT, "phase2b_mi.log"), split = TRUE)

dev <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv")
g <- dev[dev$gradeable == 1, ]
d <- g[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
d$tg_log <- log2(d$tg); d$tg <- NULL
d$age10 <- d$age / 10
cat(sprintf("MI dataset: n=%d (gradeable), events=%d, complete-cases=%d\n",
    nrow(d), sum(d$any_dr), sum(complete.cases(d))))
md <- mice::md.pattern(d, plot = FALSE)
print(md)

imp <- mice(d, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905,
            maxit = 10)
cat("imputation done\n")

fits <- with(imp, glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c +
                        apo_a1 + apo_b, family = binomial()))
pooled <- pool(fits)
s <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)
cat("\n==== pooled MI results (OR scale) ====\n")
print(s)
write.csv(s, file.path(OUT, "phase2b_mi_pooled.csv"), row.names = FALSE)

# pooled AUC across imputations
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b
auc_fun <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aucs <- sapply(1:20, function(i) {
  di <- complete(imp, i)
  fi <- glm(fml, data = di, family = binomial())
  auc_fun(predict(fi, type = "response"), di$any_dr)
})
cat(sprintf("\npooled MI AUC: mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
    mean(aucs), sd(aucs), min(aucs), max(aucs)))
write.csv(data.frame(imputation = 1:20, auc = aucs),
          file.path(OUT, "phase2b_mi_auc.csv"), row.names = FALSE)

cat("\n==== comparison: listwise vs MI ====\n")
cat(sprintf("listwise: n=1,207, events=121, corrected AUC=0.669\n"))
cat(sprintf("MI:       n=1,966, events=%d, mean AUC=%.3f\n", sum(d$any_dr), mean(aucs)))
cat("\n==== ALL DONE ====\n")
sink()
