# six-target cis-MR: server-side vs local PLINK2 LD clumping cross-check.
#
# Purpose
#   The committed six-target instruments (2_cis_mr/phase5a_cis_iv.R) were obtained
#   with the OpenGWAS LD service (ieugwasr::ld_clump, no local bfile).  This script
#   re-derives the candidate lists for the same six regions and clumps them twice:
#   (i) through the same server-side path, and (ii) locally with PLINK 2 against the
#   1000 Genomes Phase 3 EUR panel, using identical thresholds (p < 1e-5,
#   r2 < 0.01 within 500 kb).  It then diffs the two sets and both against the
#   instrument list the manuscript reports (iv_druggable.csv).
#
# Requires
#   OPENGWAS_JWT in ~/.Renviron ; PLINK2_BIN and G1000_EUR below must exist.
#
# Run
#   Rscript 2_cis_mr/phase19b_six_target_ld_clump_crosscheck.R

suppressMessages(library(ieugwasr))

OUT  <- "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables"
COMMITTED <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out/iv_druggable.csv"
PLINK2_BIN <- "E:/R_document/plink2_win64_20260110/plink2.exe"
G1000_EUR  <- "E:/workbuddy/2026-09-13-06-39-47/ld/g1000_eur"
WORK <- "E:/workbuddy/2026-10-14-16-04-57/.workbuddy/tmp/sixtarget_clump"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

targets <- data.frame(
  target = c("HMGCR", "NPC1L1", "PCSK9", "APOC3", "LPL", "ANGPTL3"),
  region = c("5:74124000-75197000", "7:44052000-45259000", "1:55005221-56028551",
             "11:116162381-117164780", "8:19341030-20369239", "1:62092402-63184930"),
  lipid  = c("LDL", "LDL", "LDL", "TG", "TG", "TG"),
  stringsAsFactors = FALSE)
lipid_id <- c(LDL = "ieu-a-300", TG = "ieu-a-302")

committed <- NULL
if (file.exists(COMMITTED)) {
  cm <- read.csv(COMMITTED, stringsAsFactors = FALSE)
  committed <- split(cm$rsid, cm$target)
}

cmp <- data.frame()
for (i in seq_len(nrow(targets))) {
  tg <- targets$target[i]; reg <- targets$region[i]; lid <- lipid_id[[targets$lipid[i]]]
  cat(sprintf("\n== %s  %s  (%s, %s)\n", tg, reg, targets$lipid[i], lid))

  a <- tryCatch(associations(reg, lid), error = function(e) {
    cat("   fetch ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(a)) { next }
  a <- a[a$p < 1e-5, ]
  cat(sprintf("   candidates p<1e-5: %d\n", nrow(a)))
  if (nrow(a) == 0) { next }

  # ---- (i) server-side ----
  cd  <- data.frame(rsid = a$rsid, pval = a$p, id = lid)
  srv <- tryCatch(ld_clump(cd, clump_kb = 500, clump_r2 = 0.01, pop = "EUR"),
                  error = function(e) { cat("   server ERR:", conditionMessage(e), "\n"); NULL })
  srv_set <- if (is.null(srv)) character(0) else sort(unique(srv$rsid))

  # ---- (ii) local PLINK 2, identical thresholds ----
  cin  <- file.path(WORK, sprintf("%s_in.txt", tg))
  cpre <- file.path(WORK, sprintf("%s_clump", tg))
  write.table(data.frame(SNP = a$rsid, P = a$p), cin,
              row.names = FALSE, quote = FALSE, sep = "\t")
  cf   <- paste0(cpre, ".clumps")
  if (file.exists(cf)) file.remove(cf)
  pargs <- c("--bfile", G1000_EUR, "--clump", cin,
             "--clump-id-field", "SNP", "--clump-p-field", "P",
             "--clump-kb", "500", "--clump-p1", "1e-5",
             "--clump-r2", "0.01", "--out", cpre)
  cmd_out <- suppressWarnings(system2(PLINK2_BIN, pargs, stdout = TRUE, stderr = TRUE))
  loc_set <- character(0)
  if (file.exists(cf)) {
    # NB: PLINK2's .clumps header starts with '#', which read.table would treat as a
    # comment and silently drop, leaving the first data row as the header (0 rows).
    cl <- read.table(cf, header = TRUE, comment.char = "", stringsAsFactors = FALSE,
                     check.names = FALSE)
    loc_set <- sort(unique(cl$ID))
  } else {
    cat("   !! local PLINK2 produced no .clumps (status ",
        format(attr(cmd_out, "status")), ")\n", sep = "")
    cat(paste(tail(as.character(cmd_out), 14), collapse = "\n"), "\n\n")
  }

  cat("   server :", paste(srv_set, collapse = ", "), "\n")
  cat("   local  :", paste(loc_set, collapse = ", "), "\n")
  cat("   identical set:", setequal(srv_set, loc_set), "\n")
  if (!is.null(committed[[tg]])) {
    cat("   committed:", paste(sort(committed[[tg]]), collapse = ", "), "\n")
    cat("   server == committed:", setequal(srv_set, committed[[tg]]),
        "| local == committed:", setequal(loc_set, committed[[tg]]), "\n")
  }

  add <- function(set, src) {
    if (length(set) == 0) set <- NA_character_
    data.frame(target = tg, lipid = targets$lipid[i], n_candidates = nrow(a),
               rsid = set, source = src, stringsAsFactors = FALSE)
  }
  cmp <- rbind(cmp, add(srv_set, "OpenGWAS server"), add(loc_set, "local PLINK2"))
  Sys.sleep(1)
}

if (nrow(cmp)) {
  p <- file.path(OUT, "cis_six_target_ld_clump_crosscheck.csv")
  write.csv(cmp, p, row.names = FALSE)
  cat("\nwritten:", p, "\n")
}
cat("\n==== ALL DONE ====\n")
