# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Reproduces the analyses added in response to Reviewer 1.
# Inputs are the extracts written by tools/revision_20261030_extract_union.sh; the
# GLGC 2013 summary statistics are fetched by tools/revision_20261030_download_glgc.py.
suppressMessages(library(MVMR))
d <- read.csv("union_joint_set.csv")
h <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv", stringsAsFactors = FALSE)
tg_inst <- unique(h$SNP[h$exposure == "TG" & h$mr_keep])
cat("TG instruments:", length(tg_inst), "; in joint set:", sum(d$rsid %in% tg_inst), "\n")

## one-exposure MVMR on the TG instruments must reproduce the published IVW (0.86)
s <- d[d$rsid %in% tg_inst, ]
r1 <- format_mvmr(BXGs = as.matrix(s[, "bT", drop = FALSE]), BYG = s$bY,
                  seBXGs = as.matrix(s[, "sT", drop = FALSE]), seBYG = s$sY, RSID = s$rsid)
cat("\none-exposure MVMR on TG instruments (should be 0.86):\n")
print(exp(ivw_mvmr(r1)[1, 1:2]))

## three-exposure model, same SNPs as the published univariable analysis
r3 <- format_mvmr(BXGs = as.matrix(d[, c("bT", "bL", "bH")]), BYG = d$bY,
                  seBXGs = as.matrix(d[, c("sT", "sL", "sH")]), seBYG = d$sY, RSID = d$rsid)
iv <- ivw_mvmr(r3)
cat("\nMVMR on the joint set:\n")
print(cbind(OR = exp(iv[, 1]), lo = exp(iv[, 1] - 1.96 * iv[, 2]),
            hi = exp(iv[, 1] + 1.96 * iv[, 2]), p = iv[, 4]))

## also: univariable IVW on the 55 TG instruments using the joint-set columns
w <- 1 / s$sY^2
b <- sum(s$bT * s$bY * w) / sum(s$bT^2 * w)
se <- sqrt(1 / sum(s$bT^2 * w))
cat(sprintf("\nmanual IVW on the same %d SNPs: OR %.3f (%.3f-%.3f)\n", nrow(s), exp(b),
            exp(b - 1.96 * se), exp(b + 1.96 * se)))
