OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(OUT, "phase1_mvmr.log"), split = TRUE)
suppressMessages(library(TwoSampleMR))

d <- read.csv(file.path(OUT, "harmonised_data.csv"))
d <- d[d$mr_keep, ]
cat("harmonised kept rows:", nrow(d), "\n")

build_mvdat <- function(sub, exp_list) {
  snps_by_exp <- lapply(exp_list, function(e) unique(sub$SNP[sub$exposure == e]))
  common <- Reduce(intersect, snps_by_exp)
  cat("  common instruments across", paste(exp_list, collapse = "+"), ":", length(common), "\n")
  if (length(common) < (length(exp_list) + 2)) {
    cat("  !! too few common SNPs for stable MVMR\n")
  }
  sub <- sub[sub$SNP %in% common & sub$exposure %in% exp_list, ]
  oc <- sub[!duplicated(sub$SNP), c("SNP", "beta.outcome", "se.outcome", "outcome", "id.outcome")]
  wbeta <- wse <- wpval <- matrix(NA_real_, nrow = length(common), ncol = length(exp_list))
  rownames(wbeta) <- common
  for (i in seq_along(exp_list)) {
    e <- exp_list[i]
    rows <- sub[sub$exposure == e, ]
    m <- match(common, rows$SNP)
    wbeta[, i] <- rows$beta.exposure[m]
    wse[, i] <- rows$se.exposure[m]
    wpval[, i] <- rows$pval.exposure[m]
  }
  colnames(wbeta) <- colnames(wse) <- colnames(wpval) <- exp_list
  stopifnot(all(oc$SNP == rownames(wbeta)))
  list(
    outcome_beta = oc$beta.outcome,
    outcome_se = oc$se.outcome,
    exposure_beta = wbeta,
    exposure_se = wse,
    exposure_pval = wpval,
    expname = data.frame(id.exposure = exp_list, exposure = exp_list, stringsAsFactors = FALSE),
    outname = data.frame(id.outcome = "DR", outcome = "Diabetic retinopathy (FinnGen R13)", stringsAsFactors = FALSE)
  )
}

models <- list(
  core = c("TG", "LDL", "HDL"),
  full = c("TG", "LDL", "HDL", "TC", "ApoA1", "ApoB")
)

allres <- list()
for (nm in names(models)) {
  cat("\n==== Model:", nm, "====\n")
  mvdat <- build_mvdat(d, models[[nm]])
  res <- mv_multiple(mvdat, intercept = FALSE)
  r <- res$result
  r$or <- exp(r$b)
  r$or_lci95 <- exp(r$b - 1.96 * r$se)
  r$or_uci95 <- exp(r$b + 1.96 * r$se)
  print(r[, c("exposure", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")])
  write.csv(r, file.path(OUT, paste0("mvmr_", nm, ".csv")), row.names = FALSE)
  allres[[nm]] <- r
}

cat("\n==== done ====\n")
sink()
