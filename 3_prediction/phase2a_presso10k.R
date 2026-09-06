# P2a: MR-PRESSO with NbDistribution = 10,000 (publication-grade robustness)
# MRPRESSO 1.0 interface: mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, ...)
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
suppressMessages(library(MRPRESSO))
sink(file.path(OUT, "presso10k.log"), split = TRUE)

h <- read.csv(file.path(OUT, "harmonised_data.csv"))
exps <- c("TG", "LDL", "HDL", "TC", "ApoA1", "ApoB")
summ <- data.frame()
for (nm in exps) {
  subd <- h[h$exposure == nm & h$mr_keep,
            c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
  if (nrow(subd) < 4) { cat(sprintf("[%s] skipped (n<4)\n", nm)); next }
  cat(sprintf("[%s] n=%d PRESSO Nb=10000 running ...", nm, nrow(subd)))
  flush.console()
  t0 <- Sys.time()
  pr <- tryCatch(
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              data = subd, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              SignifThreshold = 0.05, NbDistribution = 10000, seed = 20260905),
    error = function(e) { cat(" ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(pr)) next
  saveRDS(pr, file.path(OUT, paste0("presso10k_", nm, ".rds")))
  parse_p <- function(s) {
    s <- as.character(s)
    ifelse(grepl("<", s), 0.00005, suppressWarnings(as.numeric(s)))
  }
  mres <- pr[["MR-PRESSO results"]]
  gt <- mres[["Global Test"]]                      # list: $RSSobs, $Pvalue (string)
  pg_str <- as.character(gt$Pvalue)
  pg <- parse_p(pg_str)
  main <- pr[["Main MR results"]]                  # cols: MR Analysis / Causal Estimate / Sd / P-value
  cori <- which(main[["MR Analysis"]] == "Outlier-corrected")
  if (!length(cori)) cori <- 1
  bcorr <- as.numeric(main[["Causal Estimate"]][cori])
  pcorr <- as.numeric(main[["P-value"]][cori])
  ot <- mres[["Outlier Test"]]                     # cols: RSSobs, Pvalue (string)
  nout <- sum(parse_p(ot$Pvalue) < 0.05, na.rm = TRUE)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  summ <- rbind(summ, data.frame(exposure = nm, n_snp = nrow(subd),
      global_p = pg, global_p_raw = pg_str, corrected_or = exp(bcorr),
      corrected_p = pcorr, n_outliers = nout, minutes = mins))
  cat(sprintf(" done [%.1f min] global P=%s  corrected OR=%.3f (p=%s)  outliers=%s\n",
      mins, format(pg, digits = 3), exp(bcorr), format(pcorr, digits = 3), nout))
  flush.console()
  write.csv(summ, file.path(OUT, "presso10k_summary.csv"), row.names = FALSE)
}
cat("\n==== ALL DONE ====\n")
print(summ)
sink()
