# Run the six layer-1 lipid instrument sets against an external outcome GWAS subset.
# Column names of the external file are supplied through environment variables so the
# same script serves the replication outcome and the negative-control outcome.
suppressMessages(library(TwoSampleMR))

harm  <- Sys.getenv("LDF_HARM", "E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv")
extf  <- Sys.getenv("LDF_EXT")                       # tab-separated subset of the external GWAS
lab   <- Sys.getenv("LDF_EXT_LABEL", "external")
outdir<- Sys.getenv("LDF_OUT", "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")
nsnp_col <- Sys.getenv("LDF_N", "")
c_snp <- Sys.getenv("LDF_C_SNP"); c_ea <- Sys.getenv("LDF_C_EA"); c_oa <- Sys.getenv("LDF_C_OA")
c_b   <- Sys.getenv("LDF_C_B");   c_se <- Sys.getenv("LDF_C_SE"); c_eaf<- Sys.getenv("LDF_C_EAF")
c_p   <- Sys.getenv("LDF_C_P")

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
names(o)[names(o) == c_snp] <- "SNP"
names(o)[names(o) == c_ea]  <- "effect_allele"
names(o)[names(o) == c_oa]  <- "other_allele"
names(o)[names(o) == c_b]   <- "beta"
names(o)[names(o) == c_se]  <- "se"
if (nzchar(c_eaf)) names(o)[names(o) == c_eaf] <- "eaf"
if (nzchar(c_p))   names(o)[names(o) == c_p]   <- "pval"
if (nzchar(nsnp_col)) names(o)[names(o) == nsnp_col] <- "samplesize"
o <- o[!duplicated(o$SNP), ]
out <- format_data(o, type = "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   eaf_col = if ("eaf" %in% names(o)) "eaf" else NULL,
                   pval_col = if ("pval" %in% names(o)) "pval" else NULL,
                   samplesize_col = if ("samplesize" %in% names(o)) "samplesize" else NULL)

h <- harmonise_data(exp, out, action = 2)
h <- h[h$mr_keep, ]
res <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res$outcome_label <- lab
z <- 1.959963984540054; zp <- 0.8416212335729143
res$MDE_OR_80pct_power <- exp((z + zp) * res$se)
res$OR <- exp(res$b); res$lo <- exp(res$b - z * res$se); res$hi <- exp(res$b + z * res$se)
keep <- c("outcome_label","exposure","method","nsnp","b","se","pval","OR","lo","hi","MDE_OR_80pct_power")
res <- res[, keep]
write.csv(res, file.path(outdir, paste0("layer1_external_", lab, ".csv")), row.names = FALSE)
cat("=== IVW results for", lab, "===\n")
print(res[res$method == "Inverse variance weighted", ], row.names = FALSE)
cat("\nharmonised SNPs per exposure:/n"); print(table(h$exposure))
