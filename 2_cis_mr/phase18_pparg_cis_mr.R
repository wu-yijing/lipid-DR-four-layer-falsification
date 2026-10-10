# PPARG cis-MR: instrument -> FinnGen R13 diabetic retinopathy, with backward
# comparison of instruments the manuscript already reports.
#
# Run:
#   Rscript 2_cis_mr/phase18_pparg_cis_mr.R
#
# Inputs
#   PPARG exposure instruments  : _ext_mr/pparg_cis/pparg_cis_ivs.tsv
#   PPARG outcome subset        : _ext_mr/pparg_cis/finn_dr_pparg_subset.tsv
#   six-target inputs (for the self-check): <P1>/iv_druggable.csv,
#                                           <P1>/finngen_drug_subset.tsv,
#                                           <P1>/druggable_mr_results.csv
#
# Outputs (output_tables/)
#   cis_drug_target_MR_PPARG.csv            estimates, same shape as the six-target table
#   layer1_cis_target_MDE_PPARG.csv         adds MDE at 80% / 90% power
#   cis_instrument_strength_PPARG.csv       F statistics
#
# Conventions (identical to phase5c / phase13 so the rows are directly comparable)
#   * harmonise by allele identity, action = 2
#   * 1 instrument -> Wald ratio; >=2 -> inverse-variance weighted
#   * MDE_OR(power) = exp((Z_alpha/2 + Z_power) * se), two-sided alpha = 0.05
suppressMessages({ library(TwoSampleMR) })
set.seed(20261010)

P1  <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
EXT <- "E:/workbuddy/_ext_mr/pparg_cis"
OUT <- "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

Z_A <- 1.959963984540054; Z80 <- 0.8416212335729143; Z90 <- 1.2815515655446004
NCASE <- 15353; NCONTROL <- 88091   # FinnGen R13 DM_RETINOPATHY

# ---- inverse-variance weighted, fixed and random effects, from ratio inputs ----
ivw <- function(b, se) {
  w  <- 1 / se^2
  b_i <- sum(w * b) / sum(w)
  se_fe <- 1 / sqrt(sum(w))
  Q  <- sum(w * (b - b_i)^2); df <- length(b) - 1
  se_re <- se_fe * sqrt(max(1, Q / df))
  list(b = b_i, se_fe = se_fe, se_re = se_re, Q = Q, Q_df = df,
       Q_p = pchisq(Q, df, lower.tail = FALSE),
       p_fe = 2 * pnorm(-abs(b_i / se_fe)),
       p_re = 2 * pnorm(-abs(b_i / se_re)))
}

mde <- function(se) c(exp((Z_A + Z80) * se), exp((Z_A + Z90) * se))

# Correlated-instrument IVW (GLS on the ratio estimates).
# Rho is the correlation between the instruments' allele dosages *oriented to the
# exposure effect alleles*, taken from the 1000G EUR reference panel.
# Covariance of the ratio estimates by the delta method, assuming the exposure and
# outcome GWAS samples are disjoint (GLGC 2013 vs FinnGen R13):
#   Cov(b_i, b_j) = rho_ij * (se_i se_j + b_i b_j (se_Xi/b_Xi)(se_Xj/b_Xj))
corr_ivw <- function(b, se, bx, sex, rho) {
  k <- length(b)
  S <- matrix(0, k, k)
  relx <- sex / bx
  for (i in 1:k) for (j in 1:k) {
    R <- if (i == j) 1 else rho
    S[i, j] <- R * (se[i] * se[j] + b[i] * b[j] * relx[i] * relx[j])
  }
  w <- tryCatch(solve(S, rep(1, k)), error = function(e) NULL)
  if (is.null(w)) return(NULL)
  bg <- sum(w * b) / sum(w); sg <- 1 / sqrt(sum(w))
  list(b = bg, se = sg, OR = exp(bg),
       lci = exp(bg - Z_A * sg), uci = exp(bg + Z_A * sg),
       p = 2 * pnorm(-abs(bg / sg)), MDE80 = exp((Z_A + Z80) * sg))
}

# ---- read PPARG inputs -------------------------------------------------------
iv <- read.delim(file.path(EXT, "pparg_cis_ivs.tsv"), stringsAsFactors = FALSE)
oc <- read.delim(file.path(EXT, "finn_dr_pparg_subset.tsv"), stringsAsFactors = FALSE)
names(oc)[names(oc) == "rsids"] <- "SNP"

cat("=== PPARG cis instruments ===\n"); print(iv[, c("lipid_trait","rsid","ea","nea","beta","se","p","F")])
cat("\n=== outcome subset (FinnGen R13 DR) ===\n")
print(oc[, c("SNP","ref","alt","beta","sebeta","pval","af_alt")])

res <- data.frame(); strength <- data.frame(); per_snp <- data.frame()
for (tr in unique(iv$lipid_trait)) {
  e <- iv[iv$lipid_trait == tr, ]
  e_dat <- data.frame(
    SNP = e$rsid, beta.exposure = as.numeric(e$beta), se.exposure = as.numeric(e$se),
    effect_allele.exposure = e$ea, other_allele.exposure = e$nea,
    eaf.exposure = as.numeric(e$eaf), pval.exposure = as.numeric(e$p),
    exposure = paste0("PPARG (", tr, ")"), id.exposure = paste0("PPARG_", tr))
  o <- oc[oc$SNP %in% e$rsid, ]
  o_dat <- data.frame(
    SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$sebeta,
    effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
    eaf.outcome = o$af_alt, pval.outcome = o$pval,
    outcome = "FinnGen_R13_DR", id.outcome = "FinnGen_R13_DR",
    ncase.outcome = NCASE, ncontrol.outcome = NCONTROL)
  dat <- harmonise_data(e_dat, o_dat, action = 2)
  keep <- dat[dat$mr_keep, ]
  cat(sprintf("\n--- PPARG / %s : %d instrument(s), %d retained after harmonisation ---\n",
              tr, nrow(dat), nrow(keep)))
  print(keep[, c("SNP","effect_allele.exposure","other_allele.exposure",
                 "beta.exposure","beta.outcome","se.outcome")])

  b_ratio <- keep$beta.outcome / keep$beta.exposure
  se_ratio <- keep$se.outcome / abs(keep$beta.exposure)
  per_snp <- rbind(per_snp, data.frame(
    target = "PPARG", lipid_trait = tr, rsid = keep$SNP,
    ea = keep$effect_allele.exposure, b = b_ratio, se = se_ratio,
    OR = exp(b_ratio), lci = exp(b_ratio - 1.96 * se_ratio),
    uci = exp(b_ratio + 1.96 * se_ratio),
    p = 2 * pnorm(-abs(b_ratio / se_ratio)),
    F = (keep$beta.exposure / keep$se.exposure)^2))
  if (nrow(keep) == 1) {
    method <- "Wald ratio"; b <- b_ratio; se_used <- se_ratio; p <- 2 * pnorm(-abs(b / se_used))
    Q <- NA; Qp <- NA; se_fe <- se_ratio; se_re <- NA
    cat(sprintf("  Wald ratio: b=%.4f se=%.4f OR=%.3f p=%.4g\n", b, se_used, exp(b), p))
  } else {
    v <- ivw(b_ratio, se_ratio)
    # Manuscript convention (established empirically by the self-check below):
    # fixed-effects IVW is the reported estimator, random effects are shown for
    # reference. Q is printed so the choice is auditable.
    method <- "Inverse variance weighted"; b <- v$b
    se_fe <- v$se_fe; se_re <- v$se_re; se_used <- v$se_fe; p <- v$p_fe
    Q <- v$Q; Qp <- v$Q_p
    cat(sprintf("  IVW  b=%.4f  se_FE=%.4f se_RE=%.4f  OR=%.3f (%.3f-%.3f) p_FE=%.4g p_RE=%.4g  Q=%.3f df=%d p=%.3f\n",
                b, v$se_fe, v$se_re, exp(b), exp(b - 1.96 * se_used), exp(b + 1.96 * se_used),
                v$p_fe, v$p_re, v$Q, v$Q_df, v$Q_p))
    # cross-check against the package estimator
    mi <- mr(keep, method_list = c("mr_ivw"))
    cat(sprintf("  [package mr_ivw] b=%.4f se=%.4f p=%.4g nsnp=%d\n",
                mi$b[1], mi$se[1], mi$pval[1], mi$nsnp[1]))
  }

  mm <- mde(se_used)
  res <- rbind(res, data.frame(
    target = "PPARG", lipid_trait = tr, n_iv = nrow(dat), n_harmonised = nrow(keep),
    method = method, b = b, se = se_used, se_FE = se_fe, se_RE = se_re,
    p = p, Q = Q, Q_p = Qp,
    OR = exp(b), or_lci95 = exp(b - 1.96 * se_used), or_uci95 = exp(b + 1.96 * se_used),
    MDE_OR_80 = mm[1], MDE_OR_90 = mm[2]))
  fs <- (as.numeric(e$beta) / as.numeric(e$se))^2
  strength <- rbind(strength, data.frame(
    target = "PPARG", lipid_trait = tr, n_iv = nrow(keep),
    mean_F = round(mean(fs), 1), min_F = round(min(fs), 1),
    max_abs_p = max(as.numeric(e$p)), beta_range = sprintf("%.3f to %.3f", min(as.numeric(e$beta)), max(as.numeric(e$beta)))))
}

mde_tab <- data.frame(
  target = res$target, lipid_trait = res$lipid_trait, n_iv = res$n_iv,
  n_harmonised = res$n_harmonised, method = res$method,
  mean_F = strength$mean_F[match(res$lipid_trait, strength$lipid_trait)],
  min_F  = strength$min_F[match(res$lipid_trait, strength$lipid_trait)],
  OR = res$OR, lci = res$or_lci95, uci = res$or_uci95, p = res$p,
  se_FE = res$se_FE, se_RE = res$se_RE, Q = res$Q, Q_p = res$Q_p,
  MDE_OR_80pct_power = res$MDE_OR_80, MDE_OR_90pct_power = res$MDE_OR_90)
mde_tab$excludes_OR_above <- mde_tab$MDE_OR_80pct_power

cat("\n=== PPARG cis-MR summary ===\n")
print(mde_tab[, c("target","lipid_trait","n_harmonised","method","mean_F","min_F","OR","lci","uci","p",
                  "se_FE","se_RE","Q","Q_p","MDE_OR_80pct_power","MDE_OR_90pct_power")], row.names = FALSE)

cat("\n=== per-instrument Wald ratios (for a 2-SNP set this is also leave-one-out) ===\n")
print(per_snp[, c("lipid_trait","rsid","ea","b","se","OR","lci","uci","p","F")], row.names = FALSE)

# ---- correlation-aware sensitivity for any multi-instrument arm --------------
# Empirical r between the PPARG/TC instruments, oriented to the exposure effect
# alleles: +0.173 (r2 = 0.030) in 1000G EUR. Reported for both signs so the reader
# can see the full sensitivity range.
RHO_TC <- 0.173
corr_sens <- data.frame()
cat("\n=== correlation-aware IVW (sensitivity) ===\n")
tc_keep <- NULL
for (tr in unique(iv$lipid_trait)) {
  e <- iv[iv$lipid_trait == tr, ]
  if (nrow(e) < 2) { cat(sprintf("%s: single instrument, not applicable\n", tr)); next }
  e_dat <- data.frame(SNP = e$rsid, beta.exposure = as.numeric(e$beta), se.exposure = as.numeric(e$se),
    effect_allele.exposure = e$ea, other_allele.exposure = e$nea, eaf.exposure = as.numeric(e$eaf),
    pval.exposure = as.numeric(e$p), exposure = tr, id.exposure = tr)
  o <- oc[oc$SNP %in% e$rsid, ]
  o_dat <- data.frame(SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$sebeta,
    effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
    eaf.outcome = o$af_alt, pval.outcome = o$pval, outcome = "DR", id.outcome = "DR",
    ncase.outcome = NCASE, ncontrol.outcome = NCONTROL)
  k <- harmonise_data(e_dat, o_dat, action = 2); k <- k[k$mr_keep, ]
  k <- k[match(e$rsid, k$SNP), ]
  b <- k$beta.outcome / k$beta.exposure; se <- k$se.outcome / abs(k$beta.exposure)
  for (rho in c(RHO_TC, -RHO_TC, 0)) {
    g <- corr_ivw(b, se, k$beta.exposure, k$se.exposure, rho)
    cat(sprintf("  %s rho=%+.3f -> OR=%.3f (%.3f-%.3f) p=%.4g MDE80=%.3f\n",
                tr, rho, g$OR, g$lci, g$uci, g$p, g$MDE80))
    corr_sens <- rbind(corr_sens, data.frame(
      target = "PPARG", lipid_trait = tr, n_iv = nrow(k), rho = rho, method = "Correlated IVW (GLS)",
      OR = g$OR, lci = g$lci, uci = g$uci, p = g$p, MDE_OR_80pct_power = g$MDE80))
  }
}

write.csv(corr_sens, file.path(OUT, "cis_pparg_correlated_ivw.csv"), row.names = FALSE)

write.csv(res,      file.path(OUT, "cis_drug_target_MR_PPARG.csv"), row.names = FALSE)
write.csv(mde_tab,  file.path(OUT, "layer1_cis_target_MDE_PPARG.csv"), row.names = FALSE)
write.csv(strength, file.path(OUT, "cis_instrument_strength_PPARG.csv"), row.names = FALSE)
write.csv(per_snp,  file.path(OUT, "cis_pparg_per_instrument_wald.csv"), row.names = FALSE)

# ---- self-check: reproduce the published six-target cis results -------------
# Also establishes empirically which SE convention the manuscript's
# layer1_cis_target_MDE table is built on (fixed or random effects), so the
# PPARG rows can be reported on the same footing.
cat("\n=== self-check: manuscript's six-target cis estimates ===\n")
iv6 <- read.csv(file.path(P1, "iv_druggable.csv"), stringsAsFactors = FALSE)
oc6 <- read.delim(file.path(P1, "finngen_drug_subset.tsv"), stringsAsFactors = FALSE)
names(oc6)[names(oc6) == "rsids"] <- "SNP"
pub <- read.csv(file.path(P1, "druggable_mr_results.csv"), stringsAsFactors = FALSE)
conv <- data.frame()
for (tg in unique(iv6$target)) {
  e <- iv6[iv6$target == tg, ]
  e_dat <- data.frame(SNP = e$rsid, beta.exposure = e$beta, se.exposure = e$se,
    effect_allele.exposure = toupper(e$ea), other_allele.exposure = toupper(e$nea),
    eaf.exposure = e$eaf, pval.exposure = e$p, exposure = tg, id.exposure = tg)
  o <- oc6[oc6$SNP %in% e$rsid, ]
  o_dat <- data.frame(SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$se,
    effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
    eaf.outcome = o$af_alt, pval.outcome = o$pval,
    outcome = "DR", id.outcome = "DR", ncase.outcome = NCASE, ncontrol.outcome = NCONTROL)
  k <- harmonise_data(e_dat, o_dat, action = 2); k <- k[k$mr_keep, ]
  b <- k$beta.outcome / k$beta.exposure; se_r <- k$se.outcome / abs(k$beta.exposure)
  if (nrow(k) == 1) { b_i <- b; se_fe <- se_r; se_re <- NA; Q <- NA_real_
  } else { v <- ivw(b, se_r); b_i <- v$b; se_fe <- v$se_fe; se_re <- v$se_re; Q <- v$Q }
  ref <- pub[pub$target == tg, ]
  which_se <- if (nrow(k) == 1) "ratio" else if (abs(se_fe - ref$se) < 1e-6) "FIXED" else if (abs(se_re - ref$se) < 1e-6) "RANDOM" else "?"
  cat(sprintf("%-9s n=%d  mine OR=%.4f se_FE=%.5f se_RE=%.5f | published OR=%.4f se=%.5f  -> %s\n",
              tg, nrow(k), exp(b_i), se_fe, se_re, ref$OR, ref$se, which_se))
  conv <- rbind(conv, data.frame(target = tg, n_iv = nrow(k), se_FE = se_fe, se_RE = se_re,
                                 published_se = ref$se, convention = which_se))
}

cat("\n==== ALL DONE ====\n")
