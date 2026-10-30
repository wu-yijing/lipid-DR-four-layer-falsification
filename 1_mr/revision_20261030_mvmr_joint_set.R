# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Reproduces the analyses added in response to Reviewer 1.
# Inputs are the extracts written by tools/revision_20261030_extract_union.sh; the
# GLGC 2013 summary statistics are fetched by tools/revision_20261030_download_glgc.py.
# Run from the repository root; writes the two MVMR tables under output_tables/.
## Multivariable MR on a jointly clumped union of instruments (R1 point 6)
## Exposure order: TG, LDL-C, HDL-C. Outcome: FinnGen R13 DM_RETINOPATHY.
suppressMessages({library(MVMR)})
setwd("E:/workbuddy/_r1_review_20261029/mvmr")

rd <- function(f) {
  x <- read.table(f, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  names(x)[1:5] <- c("rsid", "a1", "a2", "beta", "se")
  x
}
tg  <- rd("extract_TG.tsv")
ldl <- rd("extract_LDL.tsv")
hdl <- rd("extract_HDL.tsv")
fg  <- rd("extract_FG.tsv")
fg  <- fg[!duplicated(fg$rsid), ]
## FinnGen reports beta for the ALT allele; put the effect allele first so that
## the same alignment rule can be applied to every source.
fg  <- data.frame(rsid = fg$rsid, a1 = fg$a2, a2 = fg$a1, beta = fg$beta,
                  se = fg$se, stringsAsFactors = FALSE)

cat(sprintf("extracted: TG %d, LDL %d, HDL %d, FinnGen %d of 134\n",
            nrow(tg), nrow(ldl), nrow(hdl), nrow(fg)))

comp <- function(x) chartr("acgtACGT", "tgcaTGCA", x)
is_pal <- function(a1, a2) {
  (toupper(a1) == "A" & toupper(a2) == "T") | (toupper(a1) == "T" & toupper(a2) == "A") |
  (toupper(a1) == "C" & toupper(a2) == "G") | (toupper(a1) == "G" & toupper(a2) == "C")
}

## reference allele set: TG
ref <- tg[, c("rsid", "a1", "a2")]
m <- merge(ref, ldl, by = "rsid", suffixes = c("_ref", "_x"))

align <- function(beta, a1, a2, ra1, ra2) {
  same <- toupper(a1) == toupper(ra1) & toupper(a2) == toupper(ra2)
  swap <- toupper(a1) == toupper(ra2) & toupper(a2) == toupper(ra1)
  same_c <- toupper(comp(a1)) == toupper(ra1) & toupper(comp(a2)) == toupper(ra2)
  swap_c <- toupper(comp(a1)) == toupper(ra2) & toupper(comp(a2)) == toupper(ra1)
  out <- rep(NA_real_, length(beta))
  out[same | same_c] <- beta[same | same_c]
  out[swap | swap_c] <- -beta[swap | swap_c]
  out
}

build <- function(df, ref, tag) {
  x <- merge(ref, df, by = "rsid", suffixes = c("_r", "_x"))
  x$b <- align(x$beta, x$a1_x, x$a2_x, x$a1_r, x$a2_r)
  x$pal <- is_pal(x$a1_x, x$a2_x)
  cat(sprintf("%s: matched %d, aligned %d, palindromic %d\n", tag, nrow(x),
              sum(!is.na(x$b)), sum(x$pal, na.rm = TRUE)))
  x
}

T <- build(tg,  ref, "TG")
L <- build(ldl, ref, "LDL")
H <- build(hdl, ref, "HDL")
F <- build(fg,  ref, "FinnGen")

key <- Reduce(intersect, list(T$rsid, L$rsid, H$rsid, F$rsid))
d <- data.frame(rsid = key,
                bT = T$b[match(key, T$rsid)], sT = T$se[match(key, T$rsid)],
                bL = L$b[match(key, L$rsid)], sL = L$se[match(key, L$rsid)],
                bH = H$b[match(key, H$rsid)], sH = H$se[match(key, H$rsid)],
                bY = F$b[match(key, F$rsid)], sY = F$se[match(key, F$rsid)],
                pal = T$pal[match(key, T$rsid)])
d <- d[stats::complete.cases(d) & !d$pal, ]
cat(sprintf("complete, non-palindromic joint set: %d variants\n", nrow(d)))

bx  <- as.matrix(d[, c("bT", "bL", "bH")])
seb <- as.matrix(d[, c("sT", "sL", "sH")])
r   <- format_mvmr(BXGs = bx, BYG = d$bY, seBXGs = seb, seBYG = d$sY, RSID = d$rsid)

cat("\n--- conditional F (gencov fixed at 0) ---\n")
f0 <- strength_mvmr(r, gencov = 0)

## genetic correlation matrix, used both for gencov and for Q_A
pcor <- cor(bx)
cat("\n--- exposure correlation matrix (from per-SNP effects) ---\n"); print(round(pcor, 3))

glist <- lapply(seq_len(nrow(d)), function(i) {
  S <- diag(seb[i, ]); pcor %*% S
})
cat("\n--- conditional F (gencov from the estimated exposure correlation) ---\n")
f1 <- try(strength_mvmr(r, gencov = glist), silent = TRUE)
if (inherits(f1, "try-error")) cat("gencov version failed\n") else print(f1)

cat("\n--- MVMR-IVW (mutually adjusted, common instrument set) ---\n")
iv <- ivw_mvmr(r)
print(iv)

cat("\n--- Q_A heterogeneity ---\n")
q <- try(qhet_mvmr(r, pcor, CI = 0.95, iterations = 1000), silent = TRUE)
if (inherits(q, "try-error")) { cat("qhet failed:", as.character(q), "\n") } else print(q)

write.csv(d, "union_joint_set.csv", row.names = FALSE)
saveRDS(list(f0 = f0, f1 = f1, iv = iv, q = q, pcor = pcor), "mvmr_union_results.rds")
cat("\nsaved union_joint_set.csv and mvmr_union_results.rds\n")
