# 7c2: ΔAUC bootstrap 95% CI —— A 规格（=phase2b 预注册规格，保留 age）
# 点估计 ΔAUC = 0.669(CC, optimism-corrected) - 0.549(MI) = 0.120
# bootstrap: 对全可分级人群(n=1,966)重采样；每个 replicate 内按同一流程
#   重算 AUC_cc*(apparent) 与 AUC_mi*(m=20 平均)，CC 侧减固定 optimism 0.023
suppressMessages({ library(mice) })
P2 <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase7_out"
sink(file.path(OUT, "phase7c2_deltaauc.log"), split = TRUE)

OPT_CC <- 0.023
B <- 300

dev <- read.csv(file.path(P2, "dev_set.csv"), stringsAsFactors = FALSE)
g <- dev[dev$gradeable == 1, ]

mkA <- function(dat) {
  d <- dat[, c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr")]
  d$tg_log <- log2(d$tg); d$tg <- NULL
  d$age10 <- d$age / 10
  d
}
d0 <- mkA(g)

auc_fun <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fml <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b

## point estimates (must reproduce 0.692 / 0.549)
cc0 <- d0[complete.cases(d0), ]
f0 <- glm(fml, data = cc0, family = binomial())
auc_cc0 <- auc_fun(predict(f0, type = "response"), cc0$any_dr)
imp0 <- mice(d0, m = 20, method = "pmm", printFlag = FALSE, seed = 20260905, maxit = 10)
auc_mi0 <- mean(sapply(1:20, function(i) {
  di <- complete(imp0, i)
  fi <- glm(fml, data = di, family = binomial())
  auc_fun(predict(fi, type = "response"), di$any_dr)
}))
delta0 <- (auc_cc0 - OPT_CC) - auc_mi0
cat(sprintf("POINT: AUC_cc_app=%.4f (corr %.4f) AUC_mi=%.4f delta=%.4f\n",
            auc_cc0, auc_cc0 - OPT_CC, auc_mi0, delta0))

set.seed(20260905)
res <- data.frame()
t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample(nrow(d0), replace = TRUE)
  db <- d0[idx, ]
  cc <- db[complete.cases(db), ]
  ok <- FALSE
  if (nrow(cc) >= 50 && sum(cc$any_dr) >= 10 && sum(cc$any_dr) < nrow(cc)) {
    f_cc <- try(glm(fml, data = cc, family = binomial()), silent = TRUE)
    imp <- try(mice(db, m = 20, method = "pmm", printFlag = FALSE,
                    seed = 20260905 + b, maxit = 10), silent = TRUE)
    if (!inherits(f_cc, "try-error") && !inherits(imp, "try-error")) {
      auc_cc <- auc_fun(predict(f_cc, type = "response"), cc$any_dr)
      aucs <- sapply(1:20, function(i) {
        di <- complete(imp, i)
        fi <- try(glm(fml, data = di, family = binomial()), silent = TRUE)
        if (inherits(fi, "try-error")) NA else auc_fun(predict(fi, type = "response"), di$any_dr)
      })
      if (sum(!is.na(aucs)) >= 15) {
        res <- rbind(res, data.frame(rep = b, auc_cc_apparent = auc_cc,
                                     auc_mi = mean(aucs, na.rm = TRUE),
                                     delta = (auc_cc - OPT_CC) - mean(aucs, na.rm = TRUE)))
        ok <- TRUE
      }
    }
  }
  if (b %% 10 == 0 || b == B) {
    cat(sprintf("[%d/%d] ok=%d elapsed=%.1f min\n", b, B, nrow(res),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
    write.csv(res, file.path(OUT, "phase7c2_boot_partial.csv"), row.names = FALSE)
  }
  if (!ok) cat(sprintf("  rep %d skipped\n", b))
}

cat("\n==== bootstrap summary (A-spec) ====\n")
cat(sprintf("reps: %d/%d | point delta=%.4f | boot median=%.4f | boot SE=%.4f\n",
            nrow(res), B, delta0, median(res$delta), sd(res$delta)))
ci <- quantile(res$delta, probs = c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("95%% percentile CI: [%.4f, %.4f]\n", ci[1], ci[2]))
sum_tab <- data.frame(n_rep = nrow(res), point_delta = delta0,
                      boot_median = median(res$delta), boot_se = sd(res$delta),
                      ci_lower = ci[1], ci_upper = ci[2],
                      pct_delta_gt_0 = mean(res$delta > 0))
print(sum_tab)
write.csv(sum_tab, file.path(OUT, "phase7c2_deltaauc_summary.csv"), row.names = FALSE)
write.csv(res, file.path(OUT, "phase7c2_boot_full.csv"), row.names = FALSE)
cat("\n==== ALL DONE (7c2) ====\n")
sink()
