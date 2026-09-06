OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_discovery4.log"), split = TRUE)
suppressMessages(library(ieugwasr))
cat("==== verify apoB id ====\n")
ids <- c("ukb-d-30640_raw","ukb-d-30650_raw")
d <- tryCatch(gwasinfo(ids), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(d)) {
  cols <- intersect(c("id","trait","nsnps","population","year","author","build"), names(d))
  print(d[, cols])
  write.csv(d, file.path(OUT, "ukb_apoB_info.csv"), row.names = FALSE)
}
cat("\n==== done ====\n")
sink()
