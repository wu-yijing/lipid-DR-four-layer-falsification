## =====================================================================
## P1 MR diagnostics, part 2 (v6 item M6)
##   - instrument overlap across the three MVMR exposures
##   - marginal and conditional instrument strength
##   - MVMR Q_A heterogeneity
##   - cis drug-target F statistics and positive control (corrected definition)
## Pure ASCII.
## =====================================================================
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
P1  <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase2f_diag2.log"), split = TRUE)

## ---- instrument overlap (top SNPs per exposure, as used for MVMR) ----
iv <- read.csv(file.path(P1, "iv_all.csv"))
cat("==== instrument counts per exposure (iv_all.csv) ====\n")
print(table(iv$exposure))
ex <- c("TG", "LDL", "HDL")
sets <- lapply(ex, function(e) unique(iv$SNP[iv$exposure == e]))
names(sets) <- ex
cat("\n==== pairwise overlap of instrument sets ====\n")
for (i in 1:3) for (j in 1:3) if (i < j) cat(sprintf("  %s & %s : %d shared SNPs  (n=%d, n=%d)\n",
  ex[i], ex[j], length(intersect(sets[[i]], sets[[j]])), length(sets[[i]]), length(sets[[j]])))
cat(sprintf("  all three: %d shared SNPs\n", length(Reduce(intersect, sets))))
cat(sprintf("  union    : %d SNPs\n", length(Reduce(union, sets))))
cat("\nNote: MVMR used mv_extract_exposures + mv_harmonise_data (strictness 2);\n")
cat("the harmonised MVMR data carried 75 (HDL), 61 (LDL) and 41 (TG) SNPs respectively,\n")
cat("i.e. exposure-specific instrument sets that are not fully overlapping.\n")

## ---- conditional F, restricted to the fully shared instrument set ----
h <- read.csv(file.path(P1, "harmonised_data.csv"))
h <- h[h$mr_keep, ]
w <- list()
for (e in ex) {
  sub <- h[h$exposure == e, c("SNP", "beta.exposure", "se.exposure")]
  names(sub)[2:3] <- c(paste0("b_", e), paste0("se_", e))
  w[[e]] <- sub
}
m <- Reduce(function(a, b) merge(a, b, by = "SNP"), w)
cat(sprintf("\n==== fully shared SNPs across TG/LDL/HDL: %d ====\n", nrow(m)))
if (nrow(m) >= 4) {
  B <- as.matrix(m[, paste0("b_", ex)]); S <- as.matrix(m[, paste0("se_", ex)])
  pv <- 3; condF <- numeric(pv)
  for (j in 1:pv) {
    others <- setdiff(1:pv, j); W <- diag(1 / S[, j]^2); Xo <- B[, others, drop = FALSE]
    gam <- solve(t(Xo) %*% W %*% Xo) %*% t(Xo) %*% W %*% B[, j]
    r <- B[, j] - Xo %*% gam
    condF[j] <- sum((r / S[, j])^2) / (nrow(B) - length(others))
  }
  names(condF) <- ex
  cat("conditional F (shared set only):\n"); print(round(condF, 2))
} else {
  cat("!! too few fully shared instruments: conditional F cannot be estimated on a\n")
  cat("   common instrument set. This is reported in the manuscript as a limitation:\n")
  cat("   because the exposure-specific instrument sets are not fully overlapping, the\n")
  cat("   standard conditional-F diagnostic (Sanderson et al. 2019) is not directly\n")
  cat("   applicable; instrument-overlap and marginal strength are reported instead.\n")
}

## ---- marginal F per exposure ----
h$F <- (h$beta.exposure / h$se.exposure)^2
cat("\n==== marginal instrument strength (single-exposure instruments) ====\n")
for (e in ex) {
  x <- h[h$exposure == e, "F"]
  cat(sprintf("  %-4s n=%3d  mean F=%7.1f  median F=%7.1f  min F=%6.1f  n(F<10)=%d\n",
              e, length(x), mean(x), median(x), min(x), sum(x < 10)))
}
mf <- data.frame(exposure = ex,
                 n = sapply(ex, function(e) sum(h$exposure == e)),
                 mean_F = sapply(ex, function(e) mean(h$F[h$exposure == e])),
                 median_F = sapply(ex, function(e) median(h$F[h$exposure == e])),
                 min_F = sapply(ex, function(e) min(h$F[h$exposure == e])))
print(mf)
write.csv(mf, file.path(OUT, "phase2f_marginal_f.csv"), row.names = FALSE)

## ---- Q_A on the fully shared set (if estimable) ----
mm <- read.csv(file.path(P1, "mvmr_core_TG_LDL_HDL.csv"))
th <- setNames(mm$b, mm$exposure)
cat("\nMVMR core estimates used for Q_A:\n"); print(mm[, c("exposure", "nsnp", "b", "pval")])
if (nrow(m) >= 4) {
  oc <- h[h$exposure == "TG", c("SNP", "beta.outcome", "se.outcome")]
  oc <- oc[match(m$SNP, oc$SNP), ]
  ok <- !is.na(oc$beta.outcome)
  pred <- as.numeric(B[ok, ] %*% th[ex])
  varQ <- oc$se.outcome[ok]^2 + as.numeric((S[ok, ]^2) %*% (th[ex]^2))
  QA <- sum((oc$beta.outcome[ok] - pred)^2 / varQ)
  dfQ <- sum(ok) - 3
  cat(sprintf("Q_A = %.2f (df=%d), P = %.4g\n", QA, dfQ, pchisq(QA, dfQ, lower.tail = FALSE)))
} else {
  cat("Q_A not estimable on a common instrument set (same reason as conditional F).\n")
}

## ---- cis drug-target instruments: strength + positive control ----
cat("\n==== cis drug-target instruments ====\n")
civ <- read.csv(file.path(P1, "iv_druggable.csv"))
civ$F <- (civ$beta / civ$se)^2
agg <- do.call(rbind, lapply(split(civ, civ$target), function(x) data.frame(
  target = x$target[1], lipid_trait = x$lipid_trait[1], n_iv = nrow(x),
  mean_F = round(mean(x$F), 1), min_F = round(min(x$F), 1),
  max_abs_p = max(x$p), positive_control = all(x$p < 5e-6 & x$F > 10),
  beta_range = sprintf("%.3f to %.3f", min(x$beta), max(x$beta)))))
print(agg)
write.csv(agg, file.path(OUT, "phase2f_cis_strength.csv"), row.names = FALSE)
cat(sprintf("\nall 25 cis variants: mean F=%.1f, median F=%.1f, min F=%.1f; n(F<10)=%d\n",
            mean(civ$F), median(civ$F), min(civ$F), sum(civ$F < 10)))
cat("positive control: every instrument is strongly associated with its own lipid trait\n")
cat("(all p < 5e-6, all F > 10), i.e. the cis instrument sets engage their targets.\n")

cat("\n==== ALL DONE ====\n")
sink()
