# PPARG cis-MR: server-side LD clumping cross-check.
#
# Purpose
#   The PPARG instruments used in the main analysis were clumped locally with
#   PLINK 2 against the 1000 Genomes Phase 3 EUR panel, because the OpenGWAS
#   token was unavailable.  This script re-runs the same clumping through the
#   OpenGWAS LD service (ieugwasr::ld_clump without a local bfile), i.e. the code
#   path the manuscript's six-target arm used, and diffs the two instrument sets.
#
# Requires
#   OPENGWAS_JWT in ~/.Renviron, or OPENGWAS_JWT set in the environment.
#
# Run
#   Rscript 2_cis_mr/phase19_pparg_ld_clump_server.R
suppressMessages({ library(ieugwasr) })

EXT <- "E:/workbuddy/_ext_mr/pparg_cis"
OUT <- "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables"

cat("== OpenGWAS identity check ==\n")
me <- tryCatch(user(), error = function(e) e)
if (inherits(me, "error")) {
  cat("!! user() failed:", conditionMessage(me), "\n")
  if (grepl("401|token", conditionMessage(me), ignore.case = TRUE)) {
    cat("!! No usable token - aborting before any clumping call.\n")
    quit(status = 1)
  }
} else {
  print(utils::str(me))
}

# local instrument sets, for comparison
local_idx <- list(
  TG = "rs10440120",
  TC = c("rs7616006", "rs7638333"))

cmp <- data.frame()
for (trait in c("TG", "TC")) {
  inp <- file.path(EXT, sprintf("pparg_%s_clump_in.txt", trait))
  cd  <- read.delim(inp, stringsAsFactors = FALSE)   # columns SNP, P
  names(cd) <- c("rsid", "pval")
  cd$id <- trait                                     # ieugwasr expects an id column
  cat(sprintf("\n== %s: %d candidates -> server-side ld_clump(r2<0.01, 500 kb, EUR)\n",
              trait, nrow(cd)))
  cl <- tryCatch(
    ld_clump(cd, clump_kb = 500, clump_r2 = 0.01, pop = "EUR"),
    error = function(e) { cat("!! ld_clump failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(cl)) {
    cmp <- rbind(cmp, data.frame(trait = trait, rsid = NA, source = "server FAILED"))
    next
  }
  srv <- sort(unique(cl$rsid))
  loc <- sort(local_idx[[trait]])
  cat("  server :", paste(srv, collapse = ", "), "\n")
  cat("  local  :", paste(loc, collapse = ", "), "\n")
  cat("  identical set:", setequal(srv, loc), "\n")
  cmp <- rbind(cmp,
    data.frame(trait = trait, rsid = loc, source = "local PLINK2"),
    data.frame(trait = trait, rsid = srv, source = "OpenGWAS server"))
}

if (nrow(cmp)) write.csv(cmp, file.path(OUT, "cis_pparg_ld_clump_crosscheck.csv"), row.names = FALSE)
cat("\n==== ALL DONE ====\n")
