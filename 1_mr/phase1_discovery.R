# Phase 1.1/1.2 exposure & outcome GWAS discovery
# ASCII-only script (Rscript on GBK locale corrupts non-ASCII literals)
# Outputs: log + candidate CSVs in E:/workbuddy/2026-09-04-18-19-16/phase1_out
OUTDIR <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUTDIR, "phase1_discovery.log")
sink(logf, split = TRUE)

cat("==== Phase 1.1/1.2 GWAS discovery ====\n")
cat("Time:", format(Sys.time()), "\n\n")

library(ieugwasr)

au <- tryCatch(user(), error = function(e) list(error = conditionMessage(e)))
cat("[auth] "); print(au); cat("\n")

cat("[1] gwasinfo() full catalog...\n")
info <- tryCatch(gwasinfo(), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (is.null(info)) {
  Sys.sleep(30)
  info <- gwasinfo()
}
cat("    catalog rows:", nrow(info), "\n\n")

tr <- tolower(ifelse(is.na(info$trait), "", info$trait))
sel_cols <- intersect(c("id","trait","ncase","ncontrol","population","consortium","year","author","build"), names(info))

cat("[2] GLGC lipid candidates (ieu-a-29x/30x):\n")
glgc_ids <- info[grepl("^ieu-a-(29[5-9]|30[0-5])$", info$id), ]
print(glgc_ids[, sel_cols])
write.csv(glgc_ids, file.path(OUTDIR, "discovery_glgc_lipids.csv"), row.names = FALSE)

cat("\n[2b] larger lipid GWAS (ieu-b / ebi-a):\n")
big <- info[grepl("^(ieu-b|ebi-a)-", info$id) &
            grepl("low density lipoprotein|high density lipoprotein|total cholesterol|triglycerides", tr), ]
print(head(big[, sel_cols], 25))
write.csv(big, file.path(OUTDIR, "discovery_lipid_large.csv"), row.names = FALSE)

cat("\n[3] ApoA1 / ApoB candidates:\n")
apo <- info[grepl("apolipoprotein", tr), ]
print(head(apo[, sel_cols], 40))
write.csv(apo, file.path(OUTDIR, "discovery_apo.csv"), row.names = FALSE)

cat("\n[4] outcome candidates (retinopathy):\n")
ret <- info[grepl("retinopathy", tr), ]
print(ret[, sel_cols])
write.csv(ret, file.path(OUTDIR, "discovery_retinopathy.csv"), row.names = FALSE)

cat("\n[5] FinnGen diabetes endpoints (backup):\n")
fg <- info[grepl("^finn-", info$id) & grepl("diabet", tr), ]
print(head(fg[, intersect(c("id","trait","ncase","ncontrol","population","year"), names(fg))], 20))
write.csv(fg, file.path(OUTDIR, "discovery_finngen_diabetes.csv"), row.names = FALSE)

cat("\n==== discovery done ====\n")
sink()
