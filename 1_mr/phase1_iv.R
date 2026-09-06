# Phase 1.3a: extract IVs for 6 locked lipid exposures (API clump r2<0.001 kb=10000 p<5e-8)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_iv.log"), split = TRUE)
suppressMessages(library(TwoSampleMR))
cat("==== Phase 1.3a IV extraction ====\n")
cat("Time:", format(Sys.time()), "\n\n")

EXPOSURES <- list(
  LDL   = "ieu-a-300",        # GLGC 2013 LDL cholesterol
  HDL   = "ieu-a-299",        # GLGC 2013 HDL cholesterol
  TC    = "ieu-a-301",        # GLGC 2013 total cholesterol
  TG    = "ieu-a-302",        # GLGC 2013 triglycerides
  ApoA1 = "ukb-d-30630_raw",  # Neale lab UKB Apolipoprotein A (label-verified)
  ApoB  = "ukb-d-30640_raw"   # Neale lab UKB Apolipoprotein B (label-verified)
)
sleep_sec <- 10

inst_list <- list()
for (nm in names(EXPOSURES)) {
  id <- EXPOSURES[[nm]]
  cat("  -", nm, "(", id, "): ", flush = TRUE)
  d <- tryCatch({
    extract_instruments(id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
  }, error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(d) || nrow(d) == 0) { inst_list[[nm]] <- NULL; next }
  d$exposure <- nm
  d$exposure_gwas_id <- id
  inst_list[[nm]] <- d
  cat(nrow(d), "SNPs\n", flush = TRUE)
  write.csv(d, file.path(OUT, paste0("iv_", nm, ".csv")), row.names = FALSE)
  Sys.sleep(sleep_sec)
}

inst <- do.call(rbind, inst_list)
if (is.null(inst) || nrow(inst) == 0) {
  cat("NO INSTRUMENTS EXTRACTED - abort\n"); sink(); quit(status = 1)
}
write.csv(inst, file.path(OUT, "iv_all.csv"), row.names = FALSE)
# unique snp list for local outcome filtering
snps <- unique(inst$SNP)
writeLines(snps, file.path(OUT, "iv_snps.txt"))
cat("\ntotal unique IV SNPs:", length(snps), "\n")
for (nm in names(EXPOSURES)) {
  k <- sum(inst$exposure == nm)
  cat(sprintf("  %-6s %3d SNPs\n", nm, k))
}
# per-SNP F-stat
inst$F_stat <- (inst$beta.exposure / inst$se.exposure)^2
fstat <- aggregate(F_stat ~ exposure, data = inst, FUN = mean)
names(fstat)[2] <- "mean_F_stat"
cat("\nmean F-stat per exposure:\n"); print(fstat)
write.csv(fstat, file.path(OUT, "iv_mean_fstat.csv"), row.names = FALSE)

cat("\n==== IV extraction done ====\n")
sink()
