# Phase 1.1 lean discovery: query small ID batches (full catalog times out on trial tier)
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "phase1_discovery2.log"), split = TRUE)

suppressMessages(library(ieugwasr))
cat("==== lean discovery ====\n")

show <- function(d) {
  cols <- intersect(c("id","trait","ncase","ncontrol","nsnps","population","consortium","year","author","build"), names(d))
  print(d[, cols])
}

cat("\n[batch A] GLGC lipids\n")
ids_a <- c("ieu-a-295","ieu-a-296","ieu-a-297","ieu-a-298","ieu-a-299","ieu-a-300",
           "ieu-a-301","ieu-a-302","ieu-a-303","ieu-a-304","ieu-a-305")
a <- tryCatch(gwasinfo(ids_a), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(a)) { show(a); write.csv(a, file.path(OUT, "glgc_info.csv"), row.names = FALSE) }

cat("\n[batch B] ApoA1/ApoB candidates\n")
ids_b <- c("met-d-APOB","met-d-APOA1","met-d-APOA2","met-a-APOB","met-a-APOA1",
           "ukb-d-30600_raw","ukb-d-30690_raw","ukb-d-30600","ukb-d-30690",
           "ebi-a-GCST90014021","ebi-a-GCST90014023","ieu-a-1185","ieu-a-1186")
b <- tryCatch(gwasinfo(ids_b), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(b)) { show(b); write.csv(b, file.path(OUT, "apo_info.csv"), row.names = FALSE) }

cat("\n[batch C] IEU finn-a DR candidates (in case hosted)\n")
ids_c <- c("finn-a-DM_RETINOPATHY","finn-a-E4_DMRETINO","finn-a-H7_BCKRNDRETINOPAT")
cc <- tryCatch(gwasinfo(ids_c), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(cc)) { show(cc); write.csv(cc, file.path(OUT, "finn_info.csv"), row.names = FALSE) }

cat("\n==== lean discovery done ====\n")
sink()
