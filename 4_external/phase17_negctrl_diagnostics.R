# Negative-control diagnostics: MR-Egger intercept and Cochran's Q for the
# external negative-control outcome (appendicitis). Reuses the same harmonised
# exposure data and external-outcome subset as phase15.
suppressMessages(library(TwoSampleMR))

harm   <- Sys.getenv("LDF_HARM", "E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv")
extf   <- Sys.getenv("LDF_EXT")
lab    <- Sys.getenv("LDF_EXT_LABEL", "external")
outdir <- Sys.getenv("LDF_OUT", "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")

d <- read.csv(harm, stringsAsFactors = FALSE)
d <- d[d$mr_keep %in% c(TRUE, "TRUE"), ]
exp <- data.frame(SNP = d$SNP, effect_allele = d$effect_allele.exposure,
                  other_allele = d$other_allele.exposure, beta = d$beta.exposure,
                  se = d$se.exposure, eaf = d$eaf.exposure, pval = d$pval.exposure,
                  samplesize = d$samplesize.exposure, exposure = d$exposure,
                  stringsAsFactors = FALSE)
exp <- format_data(exp, type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize",
                   phenotype_col = "exposure")

o <- read.delim(extf, stringsAsFactors = FALSE, check.names = FALSE)
o <- o[!duplicated(o$SNP), ]
out <- format_data(o, type = "outcome", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
                   effect_allele_col = "EA", other_allele_col = "OA", eaf_col = "EAF",
                   pval_col = "P", samplesize_col = "N")

h <- harmonise_data(exp, out, action = 2)
h <- h[h$mr_keep, ]

ple <- mr_pleiotropy_test(h); ple$outcome_label <- lab
het <- mr_heterogeneity(h);   het$outcome_label <- lab

# Main estimates under the manuscript's reporting convention (random-effects IVW),
# because Cochran's Q is significant for every exposure.
est <- mr(h, method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median"))
z <- 1.959963984540054; zp <- 0.8416212335729143
est$OR <- exp(est$b); est$lo <- exp(est$b - z * est$se); est$hi <- exp(est$b + z * est$se)
est$MDE_OR_80pct_power <- exp((z + zp) * est$se)
est$outcome_label <- lab

# Fixed- and random-effects IVW computed analytically from the harmonised data,
# matching the convention used by 1_mr/phase12_power_mde.py
# (se_fixed = 1/sqrt(sum w); se_random = se_fixed * sqrt(max(1, Q/(k-1)))).
ivw <- do.call(rbind, lapply(split(h, h$exposure), function(g) {
  ratio <- g$beta.outcome / g$beta.exposure          # Wald ratio per variant
  w     <- g$beta.exposure^2 / g$se.outcome^2        # IVW weight (Burgess)
  b     <- sum(w * ratio) / sum(w)
  se_fe <- 1 / sqrt(sum(w))
  Q     <- sum(w * (ratio - b)^2)
  se_re <- se_fe * sqrt(max(1, Q / (nrow(g) - 1)))
  data.frame(exposure = g$exposure[1], nsnp = nrow(g), b = b,
             se_fe = se_fe, se_re = se_re, Q = Q,
             p_fe = 2 * pnorm(-abs(b / se_fe)), p_re = 2 * pnorm(-abs(b / se_re)),
             stringsAsFactors = FALSE)
}))
ivw$OR_fe <- exp(ivw$b); ivw$lo_fe <- exp(ivw$b - z * ivw$se_fe); ivw$hi_fe <- exp(ivw$b + z * ivw$se_fe)
ivw$OR_re <- exp(ivw$b); ivw$lo_re <- exp(ivw$b - z * ivw$se_re); ivw$hi_re <- exp(ivw$b + z * ivw$se_re)
ivw$MDE_fe <- exp((z + zp) * ivw$se_fe); ivw$MDE_re <- exp((z + zp) * ivw$se_re)
write.csv(ivw, file.path(outdir, paste0("layer1_negctrl_ivw_fe_re_", lab, ".csv")), row.names = FALSE)


write.csv(ple, file.path(outdir, paste0("layer1_negctrl_pleiotropy_", lab, ".csv")), row.names = FALSE)
write.csv(het, file.path(outdir, paste0("layer1_negctrl_heterogeneity_", lab, ".csv")), row.names = FALSE)
write.csv(est, file.path(outdir, paste0("layer1_negctrl_estimates_", lab, ".csv")), row.names = FALSE)

cat("== MR-Egger intercept (", lab, ") ==\n", sep = ""); print(ple)
cat("\n== Cochran Q, IVW (", lab, ") ==\n", sep = "")
print(het[het$method == "Inverse variance weighted", c("exposure", "Q", "Q_df", "Q_pval")])
cat("\n== IVW fixed vs random effects (", lab, ") ==\n", sep = "")
print(ivw[, c("exposure", "nsnp", "OR_fe", "lo_fe", "hi_fe", "p_fe", "MDE_fe",
              "OR_re", "lo_re", "hi_re", "p_re", "MDE_re")], row.names = FALSE)

