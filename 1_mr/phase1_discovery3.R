# Verify UKB ApoA1/ApoB dataset ids
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_discovery3.log"), split = TRUE)
suppressMessages(library(ieugwasr))
cat("==== verify UKB apo ids ====\n")
ids <- c("ukb-d-30620_raw","ukb-d-30630_raw","ukb-d-30620","ukb-d-30630")
d <- tryCatch(gwasinfo(ids), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(d)) {
  cols <- intersect(c("id","trait","ncase","ncontrol","nsnps","population","consortium","year","author","build"), names(d))
  print(d[, cols])
  write.csv(d, file.path(OUT, "ukb_apo_info.csv"), row.names = FALSE)
}
cat("\n==== done ====\n")
sink()
