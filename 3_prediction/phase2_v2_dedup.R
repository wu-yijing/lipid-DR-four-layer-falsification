## AUTO-GENERATED from phase2_v2.R -- sensitivity analysis: one record per patient
## (most recent screening year), input dev_set_v2_dedup.csv.

## =====================================================================
## phase2_v2.R -- layer-2 discrimination analysis on the CORRECTED outcome.
## Identical method to phase2e_redo2.R / phase2e_ipw.R (m=50, maxit=50, 3 seeds,
## pmm, B=1000 CC optimism, B=200 per-imputation MI optimism), but reads
## dev_set_v2_dedup.csv (option A: DR = explicit severity text on either eye).
## Writes *_v2 outputs to phase2_out_v2d.
## Pure ASCII. Run: Rscript phase2_v2.R
## =====================================================================
suppressMessages({ library(mice); library(pmsampsize) })
set.seed(20260912)

OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out_v2d"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2_dedup.csv"

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
  if (sum(pos) == 0 || sum(neg) == 0) return(NA_real_)
  M <- outer(p[pos], p[neg], function(a, b) (a > b) + 0.5 * (a == b))
  sum(w[pos] * as.vector(M %*% w[neg])) / (sum(w[pos]) * sum(w[neg]))
}
vif_fun <- function(dat, vars) {
  X <- as.data.frame(dat[, vars, drop = FALSE])
  sapply(vars, function(v) {
    f <- lm(as.formula(paste0("`", v, "` ~ .")), data = X[, c(v, setdiff(vars, v)), drop = FALSE])
    1 / (1 - summary(f)$r.squared)
  })
}
cond_num <- function(dat, vars) {
  M <- scale(as.matrix(dat[, vars, drop = FALSE])); s <- svd(M)$d; max(s) / min(s)
}

dev <- read.csv(DEV)
g <- dev[dev$gradeable == 1, ]
d <- data.frame(age10 = g$age / 10, male = g$male, tg_log = log2(g$tg),
                tc = g$tc, hdl_c = g$hdl_c, ldl_c = g$ldl_c,
                apo_a1 = g$apo_a1, apo_b = g$apo_b, any_dr = g$any_dr)
VARS <- c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
FML  <- as.formula(paste("any_dr ~", paste(VARS, collapse = " + ")))
cc_idx <- complete.cases(d[, VARS])
cat(sprintf("[v2] gradeable n=%d events=%d (%.1f%%); complete cases n=%d events=%d (%.1f%%)\n",
            nrow(d), sum(d$any_dr), 100 * mean(d$any_dr),
            sum(cc_idx), sum(d$any_dr[cc_idx]), 100 * mean(d$any_dr[cc_idx])))
cat(sprintf("EPV complete cases: %.2f on %d predictors (%.2f incl. intercept)\n",
            sum(d$any_dr[cc_idx]) / length(VARS), length(VARS),
            sum(d$any_dr[cc_idx]) / (length(VARS) + 1)))
cat("VIF (CC subset):\n"); print(round(vif_fun(d[cc_idx, ], VARS), 2))
cat(sprintf("condition number (CC subset) = %.1f\n", cond_num(d[cc_idx, ], VARS)))

## ---- missingness (on the corrected gradeable cohort) ----
mrate <- do.call(rbind, lapply(c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b"),
  function(v) data.frame(var = v, n_missing = sum(is.na(d[[v]])),
                         pct_missing = round(100 * mean(is.na(d[[v]])), 2),
                         missing_in_cases = sum(is.na(d[[v]]) & d$any_dr == 1),
                         missing_in_controls = sum(is.na(d[[v]]) & d$any_dr == 0))))
print(mrate, row.names = FALSE)
write.csv(mrate, file.path(OUT, "phase_v2d_missing_rates.csv"), row.names = FALSE)

## ---- complete-case model + optimism ----
dcc <- d[cc_idx, ]
fit_cc <- glm(FML, data = dcc, family = binomial())
auc_cc_app <- auc_fun(predict(fit_cc, type = "response"), dcc$any_dr)
opt <- replicate(1000, {
  s <- sample(nrow(dcc), replace = TRUE)
  f <- glm(FML, data = dcc[s, ], family = binomial())
  auc_fun(predict(f, type = "response"), dcc[s, ]$any_dr) -
    auc_fun(predict(f, newdata = dcc, type = "response"), dcc$any_dr)
})
auc_cc_cor <- auc_cc_app - mean(opt, na.rm = TRUE)
cat(sprintf("\nCC model: apparent AUC=%.4f optimism=%.4f corrected=%.4f\n",
            auc_cc_app, mean(opt, na.rm = TRUE), auc_cc_cor))
write.csv(data.frame(variable = names(coef(fit_cc)), coef = as.numeric(coef(fit_cc)),
                     OR = exp(as.numeric(coef(fit_cc))), p = summary(fit_cc)$coefficients[, 4],
                     VIF = c(NA, round(vif_fun(dcc, VARS), 2))),
          file.path(OUT, "phase_v2d_cc_model.csv"), row.names = FALSE)

## ---- MI x 3 seeds ----
seeds <- c(20260912, 20260913, 20260914)
mi_auc <- c(); mi_cc2all <- c(); pool_tabs <- list(); vif_all <- list()
for (si in seq_along(seeds)) {
  cat(sprintf("---- MI seed %d ----\n", seeds[si]))
  imp <- mice(d, m = 50, method = "pmm", maxit = 50, seed = seeds[si], printFlag = FALSE)
  s <- summary(pool(with(imp, glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
                                  family = binomial()))), conf.int = TRUE)
  pool_tabs[[si]] <- data.frame(seed = seeds[si], term = s$term, estimate = s$estimate,
                                std.error = s$std.error, p.value = s$p.value, OR = exp(s$estimate),
                                lower = exp(s$estimate - 1.96 * s$std.error),
                                upper = exp(s$estimate + 1.96 * s$std.error))
  cf <- setNames(s$estimate, s$term)
  a <- numeric(50); a2 <- numeric(50); vv <- matrix(NA, 50, length(VARS))
  for (i in 1:50) {
    di <- complete(imp, i)
    fi <- glm(FML, data = di, family = binomial())
    a[i] <- auc_fun(predict(fi, type = "response"), di$any_dr)
    lp <- cf["(Intercept)"] + as.matrix(di[, VARS]) %*% cf[VARS]
    a2[i] <- auc_fun(plogis(as.numeric(lp)), di$any_dr)
    vv[i, ] <- vif_fun(di, VARS)
  }
  colnames(vv) <- VARS; vif_all[[si]] <- colMeans(vv, na.rm = TRUE)
  mi_auc <- c(mi_auc, a); mi_cc2all <- c(mi_cc2all, a2)
  cat(sprintf("seed %d: per-imputation AUC mean=%.4f sd=%.4f range=[%.4f,%.4f] | pooled-coef AUC mean=%.4f\n",
              seeds[si], mean(a), sd(a), min(a), max(a), mean(a2)))
}
cat(sprintf("\n=== [v2] MI pooled AUC over 3x50: mean=%.4f sd=%.4f range=[%.4f,%.4f]\n",
            mean(mi_auc), sd(mi_auc), min(mi_auc), max(mi_auc)))
cat(sprintf("=== [v2] MI pooled-coefficient AUC: mean=%.4f\n", mean(mi_cc2all)))
write.csv(do.call(rbind, pool_tabs), file.path(OUT, "phase_v2d_mi_coef.csv"), row.names = FALSE)
write.csv(data.frame(variable = VARS, VIF_mean = round(vif_all[[1]], 2)),
          file.path(OUT, "phase_v2d_mi_vif.csv"), row.names = FALSE)

## ---- MI optimism ----
imp1 <- mice(d, m = 50, method = "pmm", maxit = 50, seed = seeds[1], printFlag = FALSE)
o <- sapply(1:50, function(i) {
  di <- complete(imp1, i)
  mean(replicate(200, {
    s <- sample(nrow(di), replace = TRUE)
    f <- glm(FML, data = di[s, ], family = binomial())
    auc_fun(predict(f, type = "response"), di[s, ]$any_dr) -
      auc_fun(predict(f, newdata = di, type = "response"), di$any_dr)
  }), na.rm = TRUE)
})
auc_mi_cor <- mean(mi_auc) - mean(o, na.rm = TRUE)
cat(sprintf("=== [v2] MI optimism=%.4f -> corrected MI AUC=%.4f\n", mean(o, na.rm = TRUE), auc_mi_cor))

## ---- IPW ----
d$cc <- cc_idx
ok_age <- !is.na(d$age10)
sel <- glm(cc ~ age10 + male + any_dr, data = d[ok_age, ], family = binomial())
cat("\n---- selection model P(complete cases) ----\n"); print(summary(sel)$coefficients)
cat(sprintf("selection model fitted on %d rows (%d dropped for missing age10)\n",
            sum(ok_age), sum(!ok_age)))
w_raw <- rep(NA_real_, nrow(d)); w_raw[ok_age] <- 1 / fitted(sel)
q <- quantile(w_raw[which(d$cc & ok_age)], c(0.01, 0.99), na.rm = TRUE)
dcc$w <- pmin(pmax(w_raw[cc_idx], q[1]), q[2])
dcc$w[is.na(dcc$w)] <- 1
fit_ipw <- glm(FML, data = dcc, family = binomial(), weights = w)
auc_ipw <- auc_w_fun(predict(fit_ipw, type = "response"), dcc$any_dr, dcc$w)
cat(sprintf("IPW-weighted CC AUC=%.4f (unweighted %.4f) | IPW prevalence=%.4f vs CC %.4f vs full %.4f\n",
            auc_ipw, auc_cc_app, weighted.mean(dcc$any_dr, dcc$w), mean(dcc$any_dr), mean(d$any_dr)))
msel <- summary(sel)$coefficients
write.csv(data.frame(term = rownames(msel), estimate = msel[, 1], se = msel[, 2],
                     p = msel[, 4]),
          file.path(OUT, "phase_v2d_selection_model.csv"), row.names = FALSE)

## ---- decomposition ----
dec <- data.frame(
  row = c("1_CCcoef_to_CCsubset_apparent", "2_CCcoef_to_CCsubset_corrected",
          "3_CCcoef_to_MIcohort", "4_MIcoef_to_MIcohort_apparent",
          "5_MIcoef_to_MIcohort_corrected", "6_IPW_CC_subset"),
  AUC = c(auc_cc_app, auc_cc_cor, mean(mi_cc2all), mean(mi_auc), auc_mi_cor, auc_ipw))
print(dec, row.names = FALSE)
write.csv(dec, file.path(OUT, "phase_v2d_decomposition.csv"), row.names = FALSE)

## ---- delta-AUC bootstrap (CC vs MI, same procedure as phase7c2) ----
B <- 300
set.seed(7)
pa <- predict(fit_cc, type = "response")
pmi <- rowMeans(sapply(1:50, function(i) {
  di <- complete(imp1, i)
  predict(glm(FML, data = di, family = binomial()), newdata = dcc, type = "response")
}))
delta <- replicate(B, {
  s <- sample(nrow(dcc), replace = TRUE)
  auc_fun(pa[s], dcc$any_dr[s]) - auc_fun(pmi[s], dcc$any_dr[s])
})
cat(sprintf("\n[v2] dAUC (CC - MI) on the same CC subset: mean=%.3f 95%%CI %.3f to %.3f  P(>0)=%.3f\n",
            mean(delta), quantile(delta, .025), quantile(delta, .975), mean(delta > 0)))
write.csv(data.frame(metric = c("dAUC_cc_minus_mi", "lo", "hi", "p_gt0"),
                     value = c(mean(delta), quantile(delta, .025), quantile(delta, .975), mean(delta > 0))),
          file.path(OUT, "phase_v2d_deltaauc.csv"), row.names = FALSE)

## ---- pmsampsize ----
p <- mean(d$any_dr)
grid <- expand.grid(S = 8, C = c(0.67, 0.66, 0.65, 0.61, 0.60, 0.55))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  r <- tryCatch(pmsampsize(type = "b", cstatistic = grid$C[i], parameters = grid$S[i],
                           prevalence = p, shrinkage = 0.90), error = function(e) NULL)
  if (is.null(r)) return(data.frame(S = grid$S[i], C = grid$C[i], n_required = NA, events = NA))
  data.frame(S = grid$S[i], C = grid$C[i], n_required = r$sample_size,
             events = ceiling(r$sample_size * p))
}))
res$meets_n <- res$n_required <= nrow(d)
cat(sprintf("\n[v2] pmsampsize (prevalence=%.4f, S=8):\n", p)); print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "phase_v2d_pmsampsize.csv"), row.names = FALSE)

cat("\n==== ALL DONE (v2) ====\n")
