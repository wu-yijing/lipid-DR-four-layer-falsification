# -*- coding: utf-8 -*-
# Revision round 2 (2026-11-02). Reviewer 1, point 6, second pass.
#
# The first pass supplied the jointly clumped instrument set and the conditional
# strength/heterogeneity diagnostics, but left the exposure-association
# covariance out of the model: MVMR::ivw_mvmr() accepts `gencov` and never reads
# it, and MendelianRandomization::mr_mvivw() uses `correl.x` only inside the
# conditional-F calculation. This script supplies the matrix and re-estimates.
#
# Input : 1_mr/union_joint_set.csv (shipped; 133 variants x 8 columns of public
#         GLGC 2013 and FinnGen R13 summary statistics, no individual-level
#         records). Written by 1_mr/revision_20261030_mvmr_joint_set.R from the
#         GLGC and FinnGen extracts; the copy in this repository is kept so that
#         the script runs from a clean clone without re-downloading the sources.
# Output: <outdir>/mvmr_exposure_correlation.csv
#         <outdir>/mvmr_covariance_adjusted.csv
#         <outdir>/mvmr_covariance_sensitivity.csv
#
# Run from the repository root:
#   Rscript 1_mr/revision_20261102_mvmr_covariance.R [joint_csv] [outdir]

suppressMessages(library(MVMR))

args      <- commandArgs(trailingOnly = TRUE)
joint_csv <- if (length(args) >= 1) args[1] else "1_mr/union_joint_set.csv"
outdir    <- if (length(args) >= 2) args[2] else "output_tables"
if (!file.exists(joint_csv)) {
  stop(sprintf(paste0("joint-set table not found: %s\n",
                      "Run 1_mr/revision_20261030_mvmr_joint_set.R first, or pass ",
                      "the path as the first argument."), joint_csv))
}

d   <- read.csv(joint_csv, stringsAsFactors = FALSE)
X   <- as.matrix(d[, c("bT", "bL", "bH")])      # TG, LDL-C, HDL-C
seX <- as.matrix(d[, c("sT", "sL", "sH")])
Y   <- d$bY
seY <- d$sY
n   <- nrow(d); p <- ncol(X)

## ---- point estimate: weighted least squares, identical to ivw_mvmr() ---------
w    <- 1 / seY^2
W    <- diag(w)
A    <- t(X) %*% W %*% X
th   <- solve(A) %*% t(X) %*% W %*% Y
res  <- as.numeric(Y - X %*% th)
s2   <- sum(w * res^2) / (n - p)          # equals summary(lm)$sigma^2
se_lm <- sqrt(diag(solve(A))) * sqrt(s2)

## ---- exposure correlation matrix and the covariance it implies ---------------
## Sigma_j[k, l] = R[k, l] * se_jk * se_jl, the sampling covariance of the three
## exposure estimates at variant j when the exposures come from the same GWAS.
R <- cor(X)
sigma_list <- function(Rk) lapply(seq_len(n), function(j) Rk * (seX[j, ] %o% seX[j, ]))

## ---- covariance-adjusted variance -------------------------------------------
## M-estimation: A^-1 [sum_j w_j^2 (theta' Sigma_j theta) Sigma_j] A^-1 is the
## exposure-uncertainty contribution; with Sigma_j = 0 the expression reduces to
## s2 * A^-1, i.e. the weighted-regression variance that ivw_mvmr() reports.
se_adj <- function(Rk) {
  S <- sigma_list(Rk)
  M <- matrix(0, p, p)
  for (j in seq_len(n)) {
    M <- M + w[j]^2 * as.numeric(t(th) %*% S[[j]] %*% th) * S[[j]]
  }
  sqrt(diag(s2 * solve(A) + solve(A) %*% M %*% solve(A)))
}

r_in <- format_mvmr(BXGs = X, BYG = Y, seBXGs = seX, seBYG = seY, RSID = d$rsid)

cat("variants:", n, "\n\n")
cat("=== exposure correlation matrix R ===\n"); print(round(R, 4))

cat("\n=== conditional F: gencov = 0 vs gencov = R ===\n")
f0 <- suppressWarnings(strength_mvmr(r_in, gencov = 0))
fR <- strength_mvmr(r_in, gencov = sigma_list(R))
cat(sprintf("%-10s TG %.2f  LDL-C %.2f  HDL-C %.2f\n", "gencov=0", f0[1], f0[2], f0[3]))
cat(sprintf("%-10s TG %.2f  LDL-C %.2f  HDL-C %.2f\n", "gencov=R",  fR[1], fR[2], fR[3]))

cat("\n=== conditional estimates ===\n")
out <- data.frame(exposure = c("TG", "LDL-C", "HDL-C"), b = as.numeric(th),
                  se_independence = se_lm, se_covariance = se_adj(R))
out$or_independence <- exp(out$b)
out$lci_independence <- exp(out$b - 1.96 * out$se_independence)
out$uci_independence <- exp(out$b + 1.96 * out$se_independence)
out$or_covariance <- exp(out$b)
out$lci_covariance <- exp(out$b - 1.96 * out$se_covariance)
out$uci_covariance <- exp(out$b + 1.96 * out$se_covariance)
out$p_independence <- 2 * pnorm(-abs(out$b / out$se_independence))
out$p_covariance <- 2 * pnorm(-abs(out$b / out$se_covariance))
print(out[, c("exposure", "or_independence", "lci_independence", "uci_independence",
              "p_independence", "or_covariance", "lci_covariance", "uci_covariance",
              "p_covariance")], digits = 4, row.names = FALSE)

cat("\n=== sensitivity of the adjusted SE to the assumed correlation ===\n")
mk <- function(off) { M <- matrix(off, 3, 3); diag(M) <- 1; M }
grid <- list("R = 0 (independence)" = matrix(0, 3, 3),
             "R as estimated" = R, "R x 0.5" = R * 0.5, "R x 2" = R * 2,
             "off-diagonals +0.60" = mk(0.60), "off-diagonals -0.60" = mk(-0.60))
sens <- do.call(rbind, lapply(names(grid), function(nm) {
  g <- if (all(grid[[nm]] == 0)) 0 else sigma_list(grid[[nm]])
  f <- suppressWarnings(strength_mvmr(r_in, gencov = g))
  s <- se_adj(grid[[nm]])
  data.frame(case = nm, condF_TG = f[1], condF_LDL = f[2], condF_HDL = f[3],
             se_TG = s[1], se_LDL = s[2], se_HDL = s[3])
}))
print(sens, digits = 4, row.names = FALSE)

cat("\n=== validation: Sigma = 0 must reproduce the weighted-regression SE ===\n")
cat(sprintf("se_lm    %.6f %.6f %.6f\n", se_lm[1], se_lm[2], se_lm[3]))
cat(sprintf("se_adj(0) %.6f %.6f %.6f\n", se_adj(matrix(0, 3, 3))[1],
            se_adj(matrix(0, 3, 3))[2], se_adj(matrix(0, 3, 3))[3]))
cat(sprintf("\nrelative change in SE from the covariance term: %.3e\n",
            max(abs(se_adj(R) / se_lm - 1))))

dir.create(outdir, showWarnings = FALSE)
write.csv(round(R, 4), file.path(outdir, "mvmr_exposure_correlation.csv"))
write.csv(out, file.path(outdir, "mvmr_covariance_adjusted.csv"), row.names = FALSE)
write.csv(sens, file.path(outdir, "mvmr_covariance_sensitivity.csv"), row.names = FALSE)
cat(sprintf("\nwritten to %s/: mvmr_exposure_correlation.csv, mvmr_covariance_adjusted.csv, mvmr_covariance_sensitivity.csv\n", outdir))
