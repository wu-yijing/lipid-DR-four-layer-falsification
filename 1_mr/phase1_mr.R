# Phase 1.3c-1.6: format local FinnGen outcome, harmonise, MR, sensitivity, gate export
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_mr.log"), split = TRUE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
cat("==== Phase 1.3c-1.6 MR pipeline ====\n")
cat("Time:", format(Sys.time()), "\n\n")

N_CASE <- 15353; N_CTRL <- 88091
sleep_sec <- 8

## ---- 1) load IVs ----
inst <- fread(file.path(OUT, "iv_all.csv"), data.table = FALSE)
cat("IVs loaded:", nrow(inst), "SNPs,", length(unique(inst$exposure)), "exposures\n\n")

## ---- 2) local outcome: FinnGen R13 DM_RETINOPATHY subset ----
sub <- fread(file.path(OUT, "finngen_dr_subset.tsv"), data.table = FALSE)
names(sub) <- sub("#", "", names(sub))
cat("outcome subset rows:", nrow(sub), "\n")
cat("columns:", paste(names(sub), collapse = ", "), "\n")
eaf_col <- if ("af_alt" %in% names(sub)) "af_alt" else if ("af" %in% names(sub)) "af" else NA
if (is.na(eaf_col)) stop("no eaf column found")

# expand multi-rsid rows: one outcome row per target SNP
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
cat("outcome rows after rsid expansion:", nrow(out_df), "\n")

outcome_dat <- format_data(out_df, type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  pval_col = "pval", chr_col = "chr", pos_col = "pos",
  samplesize_col = "samplesize", ncase_col = "ncase", ncontrol_col = "ncontrol")
cat("formatted outcome SNPs:", nrow(outcome_dat), "\n\n")

## ---- 3) harmonise ----
dat <- harmonise_data(exposure_dat = inst, outcome_dat = outcome_dat, action = 2)
write.csv(dat, file.path(OUT, "harmonised_data.csv"), row.names = FALSE)
cat("harmonised:", nrow(dat), "SNP-exposure pairs;",
    sum(dat$mr_keep), "kept\n")
dropped <- dat[!dat$mr_keep, ]
if (nrow(dropped)) {
  cat("dropped (not kept):", nrow(dropped), "\n")
  cat("  remove flag:", sum(dropped$remove, na.rm = TRUE),
      " palindromic:", sum(dropped$palindromic, na.rm = TRUE),
      " ambiguous:", sum(dropped$ambiguous, na.rm = TRUE), "\n")
}
cat("\n")

## ---- 4) MR main ----
cat("[1.4] MR main methods...\n")
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res <- generate_odds_ratios(res)
disp_cols <- c("exposure","outcome","method","nsnp","b","se","pval","or","or_lci95","or_uci95")
disp_cols <- intersect(disp_cols, names(res))
print(res[, disp_cols])
write.csv(res, file.path(OUT, "mr_main.csv"), row.names = FALSE)
cat("\n")

## ---- 5) sensitivity ----
cat("[1.5] sensitivity...\n")
het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
if (!is.null(het)) { print(het[, c("exposure","method","Q","Q_df","Q_pval")]); write.csv(het, file.path(OUT, "mr_heterogeneity.csv"), row.names = FALSE) }
ple <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
if (!is.null(ple)) { print(ple); write.csv(ple, file.path(OUT, "mr_egger_intercept.csv"), row.names = FALSE) }
loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
if (!is.null(loo)) write.csv(loo, file.path(OUT, "mr_leaveoneout.csv"), row.names = FALSE)

presso_all <- data.frame()
for (nm in unique(dat$exposure)) {
  subd <- dat[dat$exposure == nm & dat$mr_keep, ]
  if (nrow(subd) < 4) { cat("  PRESSO skip (nsnps<4):", nm, "\n"); next }
  cat("  PRESSO:", nm, "...", flush = TRUE)
  pr <- tryCatch(run_mr_presso(subd, NbDistribution = 1000),
                 error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(pr)) {
    g <- pr[["Main MR results"]]
    g$exposure <- nm
    pg <- tryCatch(pr[["Global Test"]]$`P-global`[1], error = function(e) NA)
    if (is.null(pg) || is.na(pg)) pg <- tryCatch(pr[["Global Test"]]$P.global[1], error = function(e) NA)
    g$global_p <- pg
    presso_all <- rbind(presso_all, g)
    cat("global P =", pg, "\n")
    write.csv(g, file.path(OUT, paste0("presso_", nm, ".csv")), row.names = FALSE)
  }
  Sys.sleep(sleep_sec)
}
if (nrow(presso_all)) write.csv(presso_all, file.path(OUT, "mr_presso_all.csv"), row.names = FALSE)
cat("\n")

## ---- 6) gate + prior export ----
cat("[1.6] gate logic...\n")
ivw <- res[res$method == "Inverse variance weighted", ]
eg  <- res[res$method == "MR Egger", ]
wm  <- res[res$method == "Weighted median", ]

gate <- data.frame()
for (nm in unique(ivw$exposure)) {
  i <- ivw[ivw$exposure == nm, ]; e <- eg[eg$exposure == nm, ]; w <- wm[wm$exposure == nm, ]
  signs <- c(sign(i$b), if (nrow(e)) sign(e$b), if (nrow(w)) sign(w$b))
  sign_consistent <- length(unique(signs)) == 1
  q_row <- if (!is.null(het)) het[het$exposure == nm & het$method == "Inverse variance weighted", ] else data.frame()
  eg_int <- if (!is.null(ple)) ple[ple$exposure == nm, ] else data.frame()
  p_presso <- NA
  if (nrow(presso_all)) {
    pr_row <- presso_all[presso_all$exposure == nm, ]
    if (nrow(pr_row) && "global_p" %in% names(pr_row)) p_presso <- pr_row$global_p[1]
  }
  pass_p     <- i$pval < 0.05
  pass_het   <- if (nrow(q_row)) q_row$Q_pval > 0.05 else NA
  pass_pleio <- if (nrow(eg_int)) eg_int$pval > 0.05 else NA
  pass_presso <- if (!is.na(p_presso)) p_presso > 0.05 else TRUE
  gate_pass <- isTRUE(pass_p) && isTRUE(sign_consistent) &&
               (is.na(pass_het) || pass_het) && (is.na(pass_pleio) || pass_pleio) && pass_presso
  gate <- rbind(gate, data.frame(
    exposure = nm, nsnps = i$nsnp,
    beta_ivw = i$b, se_ivw = i$se, p_ivw = i$pval, or_ivw = i$or,
    or_lci = i$or_lci95, or_uci = i$or_uci95,
    beta_egger = if (nrow(e)) e$b else NA, p_egger = if (nrow(e)) e$pval else NA,
    beta_wm = if (nrow(w)) w$b else NA, p_wm = if (nrow(w)) w$pval else NA,
    Q_pval = if (nrow(q_row)) q_row$Q_pval else NA,
    egger_intercept_p = if (nrow(eg_int)) eg_int$pval else NA,
    presso_global_p = p_presso,
    direction_consistent = sign_consistent,
    causal_gate = gate_pass, stringsAsFactors = FALSE))
}
print(gate)
write.csv(gate, file.path(OUT, "mr_gate_summary.csv"), row.names = FALSE)

prior <- lapply(seq_len(nrow(gate)), function(k) {
  list(exposure = gate$exposure[k], nsnps = gate$nsnps[k],
       beta = gate$beta_ivw[k], se = gate$se_ivw[k],
       p = gate$p_ivw[k], or = gate$or_ivw[k])
})
jp <- list(generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           outcome = list(gwas = "FinnGen R13 DM_RETINOPATHY", ncase = N_CASE, ncontrol = N_CTRL),
           exposures = list(LDL="ieu-a-300", HDL="ieu-a-299", TC="ieu-a-301", TG="ieu-a-302",
                            ApoA1="ukb-d-30630_raw", ApoB="ukb-d-30640_raw"),
           mr_methods = c("IVW","MR-Egger","Weighted median","MR-PRESSO"),
           causal_gate = gate$exposure[gate$causal_gate],
           causal_prior = prior)
write(toJSON(jp, auto_unbox = TRUE, pretty = TRUE), file.path(OUT, "gate_prior.json"))

cat("\n==== Phase 1 pipeline done ====\n")
cat("Time:", format(Sys.time()), "\n")
cat("causal_gate passed:", paste(gate$exposure[gate$causal_gate], collapse = ", "), "\n")
sink()
