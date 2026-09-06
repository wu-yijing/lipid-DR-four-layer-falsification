OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_mvmr2.log"), split = TRUE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(data.table))
cat("==== Phase 1.7 Multivariable MR ====\n")
cat("Time:", format(Sys.time()), "\n\n")

N_CASE <- 15353; N_CTRL <- 88091

## ---- local outcome: FinnGen R13 DM_RETINOPATHY (reuse phase1 formatting) ----
sub <- fread(file.path(OUT, "finngen_dr_subset.tsv"), data.table = FALSE)
names(sub) <- sub("#", "", names(sub))
eaf_col <- if ("af_alt" %in% names(sub)) "af_alt" else if ("af" %in% names(sub)) "af" else NA
if (is.na(eaf_col)) stop("no eaf column found")
snp_vec <- strsplit(sub$rsids, ";", fixed = TRUE)
out_df <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
  data.frame(SNP = snp_vec[[i]],
             beta = sub$beta[i], se = sub$sebeta[i], pval = sub$pval[i],
             effect_allele = toupper(sub$alt[i]), other_allele = toupper(sub$ref[i]),
             eaf = sub[[eaf_col]][i],
             chr = sub[[1]][i], pos = sub$pos[i],
             samplesize = N_CASE + N_CTRL, ncase = N_CASE, ncontrol = N_CTRL,
             stringsAsFactors = FALSE)
}))
outcome_dat <- format_data(out_df, type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  pval_col = "pval", chr_col = "chr", pos_col = "pos",
  samplesize_col = "samplesize", ncase_col = "ncase", ncontrol_col = "ncontrol")
cat("formatted outcome SNPs:", nrow(outcome_dat), "\n\n")

run_mvmr <- function(id_exposures, model_name) {
  cat("==== Model:", model_name, "====\n")
  cat("exposures:", paste(id_exposures, collapse = ", "), "\n")
  exp_dat <- mv_extract_exposures(id_exposure = id_exposures)
  cat("extracted instruments:", nrow(exp_dat), "SNP-exposure pairs\n")
  mvdat <- mv_harmonise_data(exp_dat, outcome_dat, harmonise_strictness = 2)
  cat("harmonised MVMR SNPs:", nrow(mvdat), "\n")
  res <- mv_multiple(mvdat, intercept = FALSE)
  r <- res$result
  r$or <- exp(r$b)
  r$or_lci95 <- exp(r$b - 1.96 * r$se)
  r$or_uci95 <- exp(r$b + 1.96 * r$se)
  print(r[, c("exposure", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")])
  write.csv(r, file.path(OUT, paste0("mvmr_", model_name, ".csv")), row.names = FALSE)
  cat("\n")
  invisible(r)
}

# Core model: TG + LDL + HDL (all GLGC, per 1-SD, directly comparable)
run_mvmr(c("ieu-a-302", "ieu-a-300", "ieu-a-299"), "core_TG_LDL_HDL")

Sys.sleep(10)
# Full model: 6 lipids
run_mvmr(c("ieu-a-302", "ieu-a-300", "ieu-a-299", "ieu-a-301", "ukb-d-30630_raw", "ukb-d-30640_raw"), "full_6lipid")

cat("==== done ====\n")
sink()
