## =====================================================================
## P1 re-analysis of the Chinese cohort layer (v6 audit items M4-a/b/c/d/e/f)
##   - single age coding (age10 only; the age + age10 redundancy is removed)
##   - MI m=50, maxit=50, 3 seeds; apparent + optimism-corrected pooled AUC
##   - MI coefficient table with VIF; CC model with VIF
##   - three-row AUC decomposition (CC->CC, CC->MI, MI->MI)
##   - IPW analysis of the complete-case selection
##   - Riley pmsampsize sensitivity at C = 0.60
## Pure ASCII. Run:  Rscript phase2e_redo.R
## =====================================================================
suppressMessages({ library(mice); library(pmsampsize) })
set.seed(20260912)

OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "phase2e_redo.log"), split = TRUE)

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_w_fun <- function(p, y, w) {           # weighted (IPW) AUC via weighted Mann-Whitney
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
  M <- scale(as.matrix(dat[, vars, drop = FALSE]))
  s <- svd(M)$d
  max(s) / min(s)
}

## ---------------- data: single age coding ----------------
dev <- read.csv(DEV)
g <- dev[dev$gradeable == 1, ]
d <- data.frame(age10 = g$age / 10, male = g$male, tg_log = log2(g$tg),
                tc = g$tc, hdl_c = g$hdl_c, ldl_c = g$ldl_c,
                apo_a1 = g$apo_a1, apo_b = g$apo_b, any_dr = g$any_dr)
VARS <- c("age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
FML  <- as.formula(paste("any_dr ~", paste(VARS, collapse = " + ")))
cc_idx <- complete.cases(d[, VARS])
cat(sprintf("gradeable n=%d events=%d (%.1f%%); complete cases n=%d events=%d (%.1f%%)\n",
            nrow(d), sum(d$any_dr), 100 * mean(d$any_dr),
            sum(cc_idx), sum(d$any_dr[cc_idx]), 100 * mean(d$any_dr[cc_idx])))
cat("collinearity of the design matrix (variance inflation factors):\n")
print(round(vif_fun(d[cc_idx, ], VARS), 2))
cat(sprintf("condition number (CC subset) = %.1f\n", cond_num(d[cc_idx, ], VARS)))

## ---------------- complete-case model + optimism correction ----------------
dcc <- d[cc_idx, ]
fit_cc <- glm(FML, data = dcc, family = binomial())
auc_cc_app <- auc_fun(predict(fit_cc, type = "response"), dcc$any_dr)
B <- 1000
opt <- replicate(B, {
  s <- sample(nrow(dcc), replace = TRUE)
  f <- glm(FML, data = dcc[s, ], family = binomial())
  a1 <- auc_fun(predict(f, type = "response"), dcc[s, ]$any_dr)
  a2 <- auc_fun(predict(f, newdata = dcc, type = "response"), dcc$any_dr)
  a1 - a2
})
auc_cc_cor <- auc_cc_app - mean(opt, na.rm = TRUE)
cat(sprintf("\nCC model: apparent AUC=%.4f  optimism=%.4f  corrected=%.4f\n",
            auc_cc_app, mean(opt, na.rm = TRUE), auc_cc_cor))
write.csv(data.frame(variable = names(coef(fit_cc)), coef = as.numeric(coef(fit_cc)),
                     OR = exp(as.numeric(coef(fit_cc))),
                     p = summary(fit_cc)$coefficients[, 4],
                     VIF = c(NA, round(vif_fun(dcc, VARS), 2))),
          file.path(OUT, "phase2e_cc_model.csv"), row.names = FALSE)

## ---------------- MI (m = 50, maxit = 50) across 3 seeds ----------------
seeds <- c(20260912, 20260913, 20260914)
mi_all_auc <- list(); mi_pool <- list(); mi_vif <- list(); mi_cc_to_all <- list()
for (si in seq_along(seeds)) {
  cat(sprintf("\n---- MI seed %d ----\n", seeds[si]))
  imp <- mice(d, m = 50, method = "pmm", maxit = 50, seed = seeds[si], printFlag = FALSE)
  fits <- with(imp, glm(any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b,
                        family = binomial()))
  pooled <- pool(fits)
  s <- summary(pooled, conf.int = TRUE)
  if (si == 1) {
    mi_pool[[1]] <- data.frame(term = s$term, estimate = s$estimate, std.error = s$std.error,
                               p.value = s$p.value, OR = exp(s$estimate),
                               lower = exp(s$estimate - 1.96 * s$std.error),
                               upper = exp(s$estimate + 1.96 * s$std.error))
  }
  cf_pool <- setNames(s$estimate, s$term)
  aucs <- numeric(50); aucs_pool <- numeric(50); aucs_cc2all <- numeric(50); vifs <- matrix(NA, 50, 8)
  for (i in 1:50) {
    di <- complete(imp, i)
    fi <- glm(FML, data = di, family = binomial())
    aucs[i] <- auc_fun(predict(fi, type = "response"), di$any_dr)
    lp_pool <- cf_pool["(Intercept)"] + as.matrix(di[, VARS]) %*% cf_pool[VARS]
    aucs_pool[i] <- auc_fun(plogis(as.numeric(lp_pool)), di$any_dr)
    lp_cc <- as.numeric(predict(fit_cc, newdata = di, type = "response"))
    aucs_cc2all[i] <- auc_fun(lp_cc, di$any_dr)
    vifs[i, ] <- vif_fun(di, VARS)
  }
  mi_all_auc[[si]] <- aucs; mi_cc_to_all[[si]] <- aucs_cc2all
  mi_vif[[si]] <- colMeans(vifs, na.rm = TRUE)
  cat(sprintf("seed %d: per-imputation AUC mean=%.4f sd=%.4f range=[%.4f, %.4f]\n",
              seeds[si], mean(aucs), sd(aucs), min(aucs), max(aucs)))
  cat(sprintf("         pooled-coefficient AUC mean=%.4f | CC-coefficients applied to MI cohort mean=%.4f\n",
              mean(aucs_pool), mean(aucs_cc2all)))
}

mi_auc <- unlist(mi_all_auc)
cat(sprintf("\n=== MI apparent pooled AUC across 3 seeds x 50 imputations: mean=%.4f sd=%.4f range=[%.4f,%.4f]\n",
            mean(mi_auc), sd(mi_auc), min(mi_auc), max(mi_auc)))
cat(sprintf("=== CC coefficients applied to the MI cohort: mean=%.4f\n", mean(unlist(mi_cc_to_all))))

## MI-side optimism correction (bootstrap on the first seed, all 50 imputations, B=200)
imp1 <- mice(d, m = 50, method = "pmm", maxit = 50, seed = seeds[1], printFlag = FALSE)
Bmi <- 200
opt_mi <- sapply(1:50, function(i) {
  di <- complete(imp1, i)
  mean(replicate(Bmi, {
    s <- sample(nrow(di), replace = TRUE)
    f <- glm(FML, data = di[s, ], family = binomial())
    auc_fun(predict(f, type = "response"), di[s, ]$any_dr) -
      auc_fun(predict(f, newdata = di, type = "response"), di$any_dr)
  }), na.rm = TRUE)
})
auc_mi_cor <- mean(mi_auc) - mean(opt_mi, na.rm = TRUE)
cat(sprintf("=== MI optimism (mean over 50 imputations, B=%d per imputation)=%.4f -> corrected MI AUC=%.4f\n",
            Bmi, mean(opt_mi, na.rm = TRUE), auc_mi_cor))

## ---------------- IPW: model the complete-case selection ----------------
d$cc <- cc_idx
sel <- glm(cc ~ age10 + male + any_dr, data = d, family = binomial())
cat("\n---- selection model P(complete cases) ----\n"); print(summary(sel)$coefficients)
d$w_raw <- 1 / fitted(sel)
q <- quantile(d$w_raw[d$cc], c(0.01, 0.99))
dcc$w <- pmin(pmax(d$w_raw[cc_idx], q[1]), q[2])
fit_ipw <- glm(FML, data = dcc, family = binomial(), weights = w)
auc_ipw_app <- auc_w_fun(predict(fit_ipw, type = "response"), dcc$any_dr, dcc$w)
cat(sprintf("IPW-weighted CC AUC (apparent)=%.4f  (unweighted apparent CC AUC=%.4f)\n",
            auc_ipw_app, auc_cc_app))
cat(sprintf("IPW weighted prevalence=%.4f vs observed CC prevalence=%.4f vs full-cohort prevalence=%.4f\n",
            weighted.mean(dcc$any_dr, dcc$w), mean(dcc$any_dr), mean(d$any_dr)))

## ---------------- sample size sensitivity (C = 0.60) ----------------
cat("\n---- Riley pmsampsize sensitivity ----\n")
for (cv in c(0.67, 0.60)) {
  ps <- pmsampsize(type = "b", csrsquared = 0.1, parameters = 8, prevalence = 0.14,
                   cstatistic = cv, shrinkage = 0.9)
  cat(sprintf("expected C=%.2f -> n required=%s  events=%s\n", cv,
              format(ps$sample_size), format(ps$events)))
}

## ---------------- outputs ----------------
write.csv(data.frame(seed = rep(seeds, each = 50), imputation = rep(1:50, 3), auc = mi_auc),
          file.path(OUT, "phase2e_mi_auc.csv"), row.names = FALSE)
write.csv(mi_pool[[1]], file.path(OUT, "phase2e_mi_coef.csv"), row.names = FALSE)
write.csv(data.frame(variable = names(mi_vif[[1]]), VIF_mean_over_50_imputations = round(mi_vif[[1]], 2)),
          file.path(OUT, "phase2e_mi_vif.csv"), row.names = FALSE)
write.csv(data.frame(
  row = c("1_CCcoef_to_CCsubset_apparent", "2_CCcoef_to_CCsubset_corrected",
          "3_CCcoef_to_MIcohort", "4_MIcoef_to_MIcohort_apparent",
          "5_MIcoef_to_MIcohort_corrected", "6_IPW_CC_subset"),
  AUC = c(auc_cc_app, auc_cc_cor, mean(unlist(mi_cc_to_all)), mean(mi_auc), auc_mi_cor, auc_ipw_app)),
  file.path(OUT, "phase2e_decomposition.csv"), row.names = FALSE)
cat("\n==== ALL DONE ====\n")
sink()
