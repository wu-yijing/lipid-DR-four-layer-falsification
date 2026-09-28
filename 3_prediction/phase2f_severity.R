## =====================================================================
## P1 follow-up analyses (v6 items M5 outcome coding, M6 MR strength)
##   1) full coding table of the raw severity field
##   2) per-target cis instrument F statistics + positive control
##   3) MVMR conditional instrument strength (residualised approximation,
##      Sanderson et al. 2019 / Broadbent et al. 2021 formulation)
##   4) MVMR Q_A heterogeneity
## Pure ASCII. Run: Rscript phase2f_severity.R
## =====================================================================
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
DEV <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set.csv"
P1  <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase2f_severity.log"), split = TRUE)

## ---------------- 1) severity coding table ----------------
dev <- read.csv(DEV)
cat("==== raw severity field: value counts ====\n")
tb <- table(dev$severity, useNA = "ifany"); print(tb)
tab <- data.frame(severity = names(tb), n = as.vector(tb),
                  gradeable = as.vector(tapply(dev$gradeable, dev$severity, function(x) sum(x, na.rm = TRUE))),
                  any_dr = as.vector(tapply(dev$any_dr, dev$severity, function(x) sum(x, na.rm = TRUE))),
                  dr_raw_pos = as.vector(tapply(dev$dr_raw, dev$severity, function(x) sum(x, na.rm = TRUE))),
                  pdr = as.vector(tapply(dev$pdr, dev$severity, function(x) sum(x, na.rm = TRUE))))
print(tab)
write.csv(tab, file.path(OUT, "phase2f_severity_table.csv"), row.names = FALSE)

cat("\n==== cross-tabulations ====\n")
cat("-- gradeable x any_dr --\n"); print(table(dev$gradeable, dev$any_dr))
cat("-- gradeable x dr_raw --\n"); print(table(dev$gradeable, dev$dr_raw))
cat("-- severity x gradeable --\n"); print(table(dev$severity, dev$gradeable))
cat("-- severity x any_dr --\n"); print(table(dev$severity, dev$any_dr))
cat("\nprevalence: full sample = %.4f | gradeable = %.4f\n" ,
    mean(dev$any_dr), mean(dev$any_dr[dev$gradeable == 1]))
cat("among gradeable: worst-eye grade distribution (severity) = \n")
print(table(dev$severity[dev$gradeable == 1], dev$any_dr[dev$gradeable == 1]))
cat(sprintf("gradeable n=%d  any_dr=%d  PDR among gradeable=%d\n",
            sum(dev$gradeable == 1), sum(dev$any_dr[dev$gradeable == 1]),
            sum(dev$pdr[dev$gradeable == 1], na.rm = TRUE)))
cat("excluded (not gradeable) n=%d ; of these, dr_raw positive=%d ; any_dr positive=%d\n",
    sum(dev$gradeable != 1), sum(dev$dr_raw[dev$gradeable != 1], na.rm = TRUE),
    sum(dev$any_dr[dev$gradeable != 1], na.rm = TRUE))
cat("compare age / sex / available lipids between excluded and included:\n")
inc <- dev$gradeable == 1
for (v in c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")) {
  cat(sprintf("  %-7s included mean=%.3f (missing %.1f%%) | excluded mean=%.3f (missing %.1f%%)\n",
              v, mean(dev[[v]][inc], na.rm = TRUE), 100 * mean(is.na(dev[[v]][inc])),
              mean(dev[[v]][!inc], na.rm = TRUE), 100 * mean(is.na(dev[[v]][!inc]))))
}

## ---------------- 2) cis instrument strength + positive control ----------------
cat("\n==== cis drug-target instruments ====\n")
iv <- read.csv(file.path(P1, "iv_druggable.csv"))
iv$F <- (iv$beta / iv$se)^2
agg <- do.call(rbind, lapply(split(iv, iv$target), function(x) data.frame(
  target = x$target[1], lipid_trait = x$lipid_trait[1], n_iv = nrow(x),
  mean_F = mean(x$F), min_F = min(x$F),
  positive_control = all(x$beta > 0),
  beta_range = sprintf("%.3f to %.3f", min(x$beta), max(x$beta)),
  p_range = sprintf("%.2g to %.2g", max(x$p), min(x$p)))))
print(agg)
write.csv(agg, file.path(OUT, "phase2f_cis_strength.csv"), row.names = FALSE)
cat(sprintf("\nall 25 cis variants: mean F=%.1f  min F=%.1f  max F=%.1f\n",
            mean(iv$F), min(iv$F), max(iv$F)))
cat(sprintf("variants with F<10: %d of %d\n", sum(iv$F < 10), nrow(iv)))

## ---------------- 3) MVMR conditional instrument strength ----------------
cat("\n==== MVMR core (TG + LDL + HDL): conditional instrument strength ====\n")
h <- read.csv(file.path(P1, "harmonised_data.csv"))
h <- h[h$mr_keep, ]
ex <- c("TG", "LDL", "HDL")
w <- list()
for (e in ex) {
  sub <- h[h$exposure == e, c("SNP", "beta.exposure", "se.exposure")]
  names(sub)[2:3] <- c(paste0("b_", e), paste0("se_", e))
  w[[e]] <- sub
}
m <- Reduce(function(a, b) merge(a, b, by = "SNP"), w)
cat(sprintf("shared instruments across TG/LDL/HDL: %d SNPs\n", nrow(m)))
B <- as.matrix(m[, paste0("b_", ex)])
S <- as.matrix(m[, paste0("se_", ex)])
pv <- ncol(B)
condF <- numeric(pv); names(condF) <- ex
for (j in 1:pv) {
  others <- setdiff(1:pv, j)
  W <- diag(1 / S[, j]^2)
  Xo <- B[, others, drop = FALSE]
  gam <- solve(t(Xo) %*% W %*% Xo) %*% t(Xo) %*% W %*% B[, j]
  r <- B[, j] - Xo %*% gam
  condF[j] <- sum((r / S[, j])^2) / (nrow(B) - length(others))
}
cat("conditional F (residualised approximation):\n"); print(round(condF, 1))
margF <- colSums((B / S)^2) / nrow(B)
cat("marginal F (same instrument set):\n"); print(round(margF, 1))
write.csv(data.frame(exposure = ex, marginal_F = round(margF, 1), conditional_F = round(condF, 1)),
          file.path(OUT, "phase2f_conditional_f.csv"), row.names = FALSE)

## ---------------- 4) MVMR Q_A ----------------
mm <- read.csv(file.path(P1, "mvmr_core_TG_LDL_HDL.csv"))
theta <- setNames(mm$b, c("HDL", "LDL", "TG")[match(mm$exposure, c("HDL", "LDL", "TG"))])
th <- c(TG = mm$b[mm$exposure == "TG"], LDL = mm$b[mm$exposure == "LDL"], HDL = mm$b[mm$exposure == "HDL"])
hh <- h[h$exposure == "TG", c("SNP", "beta.outcome", "se.outcome")]
for (e in c("LDL", "HDL")) {
  s2 <- h[h$exposure == e, c("SNP", "beta.outcome", "se.outcome")]
  names(s2)[2:3] <- paste0(c("by", "sey"), "_", e)
  hh <- merge(hh, s2, by = "SNP")
}
Bx <- as.matrix(m[, paste0("b_", ex)]); Sx <- as.matrix(m[, paste0("se_", ex)])
rownames(Bx) <- m$SNP; rownames(Sx) <- m$SNP
hbar <- hh[hh$SNP %in% m$SNP, ]
hbar <- hbar[match(m$SNP, hbar$SNP), ]
pred <- as.numeric(Bx %*% th[ex])
varQ <- hbar$se.outcome^2 + as.numeric((Sx^2) %*% (th[ex]^2))
QA <- sum((hbar$beta.outcome - pred)^2 / varQ)
dfQ <- nrow(Bx) - length(ex)
cat(sprintf("\nMVMR Q_A = %.1f (df=%d), P = %.4g\n", QA, dfQ, pchisq(QA, dfQ, lower.tail = FALSE)))
cat(sprintf("MVMR-Egger-type intercept check not computed (MVMR package unavailable for R 4.5.2)\n"))

cat("\n==== ALL DONE ====\n")
sink()
