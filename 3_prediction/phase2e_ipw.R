## =====================================================================
## P1: IPW analysis of the complete-case selection + sample-size sensitivity
##     (completes phase2e_redo.R, which stopped at a weights argument bug)
## =====================================================================
suppressMessages({ library(pmsampsize) })
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv"
sink(file.path(OUT, "phase2e_ipw.log"), split = TRUE)

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_w_fun <- function(p, y, w) {
  ok <- !is.na(p) & !is.na(y) & !is.na(w)
  p <- p[ok]; y <- y[ok]; w <- w[ok]
  pos <- y == 1; neg <- y == 0
  if (!any(pos) || !any(neg)) return(NA_real_)
  M <- outer(p[pos], p[neg], function(a, b) (a > b) + 0.5 * (a == b))
  sum(w[pos] * as.vector(M %*% w[neg])) / (sum(w[pos]) * sum(w[neg]))
}

dev <- read.csv(DEV); g <- dev[dev$gradeable == 1, ]
d <- data.frame(age10 = g$age/10, male = g$male, tg_log = log2(g$tg), tc = g$tc,
                hdl_c = g$hdl_c, ldl_c = g$ldl_c, apo_a1 = g$apo_a1, apo_b = g$apo_b,
                any_dr = g$any_dr)
VARS <- c("age10","male","tg_log","tc","hdl_c","ldl_c","apo_a1","apo_b")
FML  <- as.formula(paste("any_dr ~", paste(VARS, collapse = " + ")))
d$cc <- complete.cases(d[, VARS])

sel <- glm(cc ~ age10 + male + any_dr, data = d, family = binomial())
cat("==== selection model P(lipids complete) ====\n"); print(summary(sel)$coefficients)
d$w_raw <- 1 / predict(sel, newdata = d, type = "response")
q <- quantile(d$w_raw[d$cc], c(0.01, 0.99))
d$w <- pmin(pmax(d$w_raw, q[1]), q[2])

dcc <- d[d$cc, ]
fit_unw <- glm(FML, data = dcc, family = binomial())
fit_ipw <- glm(FML, data = dcc, family = binomial(), weights = dcc$w)
a_unw <- auc_fun(predict(fit_unw, type = "response"), dcc$any_dr)
a_ipw <- auc_w_fun(predict(fit_ipw, type = "response"), dcc$any_dr, dcc$w)
cat(sprintf("\nCC subset n=%d events=%d\n", nrow(dcc), sum(dcc$any_dr)))
cat(sprintf("observed CC prevalence      = %.4f\n", mean(dcc$any_dr)))
cat(sprintf("IPW-weighted prevalence     = %.4f\n", weighted.mean(dcc$any_dr, dcc$w)))
cat(sprintf("full gradeable prevalence   = %.4f\n", mean(d$any_dr)))
cat(sprintf("unweighted CC AUC           = %.4f\n", a_unw))
cat(sprintf("IPW-weighted CC AUC         = %.4f\n", a_ipw))
cat(sprintf("difference (IPW - unweighted) = %+.4f\n", a_ipw - a_unw))
cat("\nInterpretation: if complete-case selection were inflating discrimination, IPW\n")
cat("re-weighting the complete cases back towards the full cohort should pull the AUC\n")
cat("down towards the MI value; the magnitude of that pull is reported above.\n")

cat("\n==== Riley pmsampsize sensitivity ====\n")
for (cv in c(0.67, 0.60, 0.55)) {
  ps <- pmsampsize(type = "b", csrsquared = 0.1, parameters = 8, prevalence = 0.14,
                   cstatistic = cv, shrinkage = 0.9)
  cat(sprintf("expected C=%.2f -> n required=%.0f  events=%.0f\n", cv, ps$sample_size, ps$events))
}
cat(sprintf("\nobserved: MI population n=%d, events=%d, 8 parameters\n", nrow(d), sum(d$any_dr)))

write.csv(data.frame(quantity = c("CC_n","CC_events","prevalence_CC","prevalence_IPW_weighted",
                                  "prevalence_full_gradeable","AUC_unweighted_CC","AUC_IPW_weighted_CC",
                                  "IPW_minus_unweighted"),
                     value = c(nrow(dcc), sum(dcc$any_dr), mean(dcc$any_dr),
                               weighted.mean(dcc$any_dr, dcc$w), mean(d$any_dr),
                               a_unw, a_ipw, a_ipw - a_unw)),
          file.path(OUT, "phase2e_ipw.csv"), row.names = FALSE)
cat("\n==== ALL DONE ====\n")
sink()
