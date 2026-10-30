# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Reproduces the analyses added in response to Reviewer 1.
# Inputs are the extracts written by tools/revision_20261030_extract_union.sh; the
# GLGC 2013 summary statistics are fetched by tools/revision_20261030_download_glgc.py.
suppressMessages(library(MVMR))
d <- read.csv("union_joint_set.csv")
cat("N joint =", nrow(d), "\n")
bx <- as.matrix(d[, c("bT", "bL", "bH")])
by <- d$bY
se <- as.matrix(d[, c("sT", "sL", "sH")])
nm <- c("TG", "LDL", "HDL")
for (j in 1:3) {
  w <- 1 / se[, j]^2
  b <- sum(bx[, j] * by * w) / sum(bx[, j]^2 * w)
  s <- sqrt(1 / sum(bx[, j]^2 * w))
  cat(sprintf("%s IVW on the joint set: OR %.3f (%.3f-%.3f) p=%.3g\n",
              nm[j], exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s), 2 * pnorm(-abs(b / s))))
}
r <- format_mvmr(BXGs = bx, BYG = by, seBXGs = se, seBYG = d$sY, RSID = d$rsid)
cat("\nMVMR-IVW:\n")
iv <- ivw_mvmr(r)
print(iv)
cat("\nQ_A (pleiotropy_mvmr):\n")
print(pleiotropy_mvmr(r))
cat("\nConditional F (gencov = 0):\n")
print(strength_mvmr(r, gencov = 0))
