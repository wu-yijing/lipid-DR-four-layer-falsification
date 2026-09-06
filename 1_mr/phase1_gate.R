# Phase 1.6 (lightweight): build causal gate + prior from existing MR outputs
# Reuses mr_main.csv / mr_heterogeneity.csv / mr_egger_intercept.csv / mr_presso_all.csv
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
suppressMessages(library(jsonlite))
cat("==== Phase 1.6 gate + prior ====\n")

res  <- read.csv(file.path(OUT, "mr_main.csv"))
het  <- read.csv(file.path(OUT, "mr_heterogeneity.csv"))
ple  <- read.csv(file.path(OUT, "mr_egger_intercept.csv"))
pa   <- tryCatch(read.csv(file.path(OUT, "mr_presso_all.csv")), error = function(e) data.frame())
pa_global <- if (nrow(pa) && "global_p" %in% names(pa)) pa else data.frame()

ivw <- res[res$method == "Inverse variance weighted", ]
eg  <- res[res$method == "MR Egger", ]
wm  <- res[res$method == "Weighted median", ]

gate <- data.frame()
for (nm in unique(ivw$exposure)) {
  i <- ivw[ivw$exposure == nm, ]; e <- eg[eg$exposure == nm, ]; w <- wm[wm$exposure == nm, ]
  signs <- c(sign(i$b), if (nrow(e)) sign(e$b), if (nrow(w)) sign(w$b))
  dir_consistent <- length(unique(signs[!is.na(signs)])) == 1 && !any(is.na(signs))
  q_row <- het[het$exposure == nm & het$method == "Inverse variance weighted", ]
  eg_int <- ple[ple$exposure == nm, ]
  p_presso <- NA
  if (nrow(pa_global)) {
    pr <- pa_global[pa_global$exposure == nm, ]
    if (nrow(pr)) p_presso <- pr$global_p[1]
  }
  pass_signal  <- i$pval < 0.05
  pass_dir     <- isTRUE(dir_consistent)
  pass_pleio   <- if (nrow(eg_int)) eg_int$pval > 0.05 else NA
  pass_presso  <- if (!is.na(p_presso)) p_presso > 0.05 else TRUE
  het_caveat   <- if (nrow(q_row)) q_row$Q_pval < 0.05 else NA

  # Gate = significant IVW signal + directionally consistent + no directional pleiotropy.
  # Heterogeneity (Cochran Q) is near-universal in cis/ polygenic MR and is reported
  # as a caveat column, NOT a hard disqualifier (would zero out all 6 otherwise).
  gate_pass <- isTRUE(pass_signal) && isTRUE(pass_dir) &&
               (is.na(pass_pleio) || pass_pleio) && pass_presso

  gate <- rbind(gate, data.frame(
    exposure = nm,
    nsnps = i$nsnp,
    beta_ivw = round(i$b, 5), se_ivw = round(i$se, 5), p_ivw = i$pval,
    or_ivw = round(i$or, 4), or_lci = round(i$or_lci95, 4), or_uci = round(i$or_uci95, 4),
    beta_egger = if (nrow(e)) round(e$b, 5) else NA, p_egger = if (nrow(e)) e$pval else NA,
    beta_wm = if (nrow(w)) round(w$b, 5) else NA, p_wm = if (nrow(w)) w$pval else NA,
    Q_pval = if (nrow(q_row)) q_row$Q_pval else NA,
    egger_intercept_p = if (nrow(eg_int)) eg_int$pval else NA,
    presso_global_p = p_presso,
    heterogeneity_caveat = if (is.na(het_caveat)) NA else het_caveat,
    direction_consistent = dir_consistent,
    causal_gate = gate_pass, stringsAsFactors = FALSE))
}
print(gate)
write.csv(gate, file.path(OUT, "mr_gate_summary.csv"), row.names = FALSE)

prior <- lapply(seq_len(nrow(gate)), function(k) {
  list(exposure = gate$exposure[k], nsnps = gate$nsnps[k],
       beta = gate$beta_ivw[k], se = gate$se_ivw[k],
       p = gate$p_ivw[k], or = gate$or_ivw[k])
})
jp <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  note = "Phase 1 MR gating. Outcome = FinnGen R13 DM_RETINOPATHY (15353 cases / 88091 controls). Gate = IVW p<0.05 AND direction consistent (IVW/Egger/WM) AND Egger intercept p>0.05. Cochran Q heterogeneity reported as caveat (universal here, not a hard gate). MR-PRESSO unstable at Nb=1000 -> global_p NA (not used as hard gate; Egger intercept is the primary pleiotropy test).",
  outcome = list(gwas = "FinnGen R13 DM_RETINOPATHY", ncase = 15353, ncontrol = 88091),
  exposures = list(LDL="ieu-a-300", HDL="ieu-a-299", TC="ieu-a-301", TG="ieu-a-302",
                   ApoA1="ukb-d-30630_raw", ApoB="ukb-d-30640_raw"),
  mr_methods = c("IVW","MR-Egger","Weighted median","MR-PRESSO(attempted)"),
  causal_gate_passed = gate$exposure[gate$causal_gate],
  causal_prior = prior)
write(toJSON(jp, auto_unbox = TRUE, pretty = TRUE), file.path(OUT, "gate_prior.json"))

cat("\ncausal_gate passed:", paste(gate$exposure[gate$causal_gate], collapse = ", "), "\n")
cat("==== done ====\n")
