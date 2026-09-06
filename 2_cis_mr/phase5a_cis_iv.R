# Improvement 3, step 1: extract cis-region instruments for 6 lipid drug targets from GLGC
# Targets (hg19 coords, +/-500kb): HMGCR/NPC1L1/PCSK9 -> LDL; APOC3/LPL/ANGPTL3 -> TG
suppressMessages(library(ieugwasr))
set.seed(20260905)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "druggable_cis.log"), split = TRUE)

targets <- data.frame(
  target = c("HMGCR", "NPC1L1", "PCSK9", "APOC3", "LPL", "ANGPTL3"),
  region = c("5:74124000-75197000", "7:44052000-45259000", "1:55005221-56028551",
             "11:116162381-117164780", "8:19341030-20369239", "1:62092402-63184930"),
  lipid  = c("LDL", "LDL", "LDL", "TG", "TG", "TG"),
  stringsAsFactors = FALSE)
lipid_id <- c(LDL = "ieu-a-300", TG = "ieu-a-302")

iv_all <- data.frame()
for (i in seq_len(nrow(targets))) {
  tg_ <- targets$target[i]; reg <- targets$region[i]
  lid <- lipid_id[[targets$lipid[i]]]
  cat(sprintf("[%s] region %s in %s ... ", tg_, reg, lid)); flush.console()
  a <- tryCatch(associations(reg, lid), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(a)) next
  a <- a[a$p < 1e-5, ]
  cat(sprintf("%d SNPs p<1e-5; ", nrow(a)))
  if (nrow(a) == 0) { cat("skip\n"); next }
  cd <- data.frame(rsid = a$rsid, pval = a$p, id = lid)
  cl <- tryCatch(ld_clump(cd, clump_kb = 500, clump_r2 = 0.01, pop = "EUR"),
                 error = function(e) { cat("clump ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(cl) && nrow(cl) > 0) keep <- a[a$rsid %in% cl$rsid, ]
  else keep <- a[which.min(a$p), ]
  keep$target <- tg_; keep$lipid_trait <- targets$lipid[i]
  iv_all <- rbind(iv_all, keep)
  cat(sprintf("kept %d independent cis-IVs\n", nrow(keep)))
  Sys.sleep(1)
}
if (nrow(iv_all)) {
  write.csv(iv_all[, c("target", "lipid_trait", "rsid", "chr", "position", "ea", "nea",
                       "beta", "se", "p", "eaf")],
            file.path(OUT, "iv_druggable.csv"), row.names = FALSE)
  writeLines(iv_all$rsid, file.path(OUT, "iv_druggable.txt"))
  cat("\nIVs per target:\n"); print(table(iv_all$target))
  cat("total cis-IVs:", nrow(iv_all), "\n")
} else cat("NO IVs extracted\n")
sink()
