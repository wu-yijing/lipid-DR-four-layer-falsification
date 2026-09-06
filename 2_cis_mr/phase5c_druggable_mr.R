# Improvement 3, step 3: drug-target cis-MR (lipid level -> FinnGen R13 DR), per target
suppressMessages({ library(TwoSampleMR) })
set.seed(20260905)
P1 <- "E:/workbuddy/2026-09-04-18-19-16/phase1_out"
sink(file.path(P1, "druggable_mr.log"), split = TRUE)

iv <- read.csv(file.path(P1, "iv_druggable.csv"), stringsAsFactors = FALSE)
oc <- read.delim(file.path(P1, "finngen_drug_subset.tsv"), stringsAsFactors = FALSE)
names(oc)[names(oc) == "rsids"] <- "SNP"

targets <- unique(iv$target)
res <- data.frame()
for (tg_ in targets) {
  e <- iv[iv$target == tg_, ]
  e_dat <- data.frame(
    SNP = e$rsid, beta.exposure = e$beta, se.exposure = e$se,
    effect_allele.exposure = toupper(e$ea), other_allele.exposure = toupper(e$nea),
    exposure = paste0(tg_, "_", e$lipid_trait[1]), id.exposure = paste0(tg_, "_", e$lipid_trait[1]),
    pval.exposure = e$p, eaf.exposure = e$eaf)
  o <- oc[oc$SNP %in% e$rsid, ]
  o_dat <- data.frame(
    SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$se,
    effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
    outcome = "FinnGen_R13_DR", id.outcome = "FinnGen_R13_DR",
    pval.outcome = o$pval, eaf.outcome = o$af_alt,
    ncase.outcome = 15353, ncontrol.outcome = 88091)
  dat <- harmonise_data(e_dat, o_dat, action = 2)
  keep <- dat[dat$mr_keep, ]
  nkept <- nrow(keep)
  if (nkept == 0) {
    res <- rbind(res, data.frame(target = tg_, lipid_trait = e$lipid_trait[1],
      n_iv = nrow(dat), n_harmonised = 0, method = NA, b = NA, se = NA, p = NA,
      OR = NA, or_lci95 = NA, or_uci95 = NA))
    cat(sprintf("[%s] 0 harmonised SNPs - skipped\n", tg_)); next
  }
  if (nkept == 1) {
    r <- data.frame(method = "Wald ratio", b = keep$beta.outcome / keep$beta.exposure,
      se = sqrt(keep$se.outcome^2 / keep$beta.exposure^2))
    r$p <- 2 * pnorm(-abs(r$b / r$se)); r$nsnp <- 1
  } else {
    mi <- mr(keep, method_list = c("mr_ivw"))
    r <- data.frame(method = mi$method[1], b = mi$b[1], se = mi$se[1], p = mi$pval[1],
                    nsnp = mi$nsnp[1])
  }
  res <- rbind(res, data.frame(target = tg_, lipid_trait = e$lipid_trait[1],
    n_iv = nrow(dat), n_harmonised = nkept, method = r$method,
    b = r$b, se = r$se, p = r$p, OR = exp(r$b),
    or_lci95 = exp(r$b - 1.96 * r$se), or_uci95 = exp(r$b + 1.96 * r$se)))
  cat(sprintf("[%s] n_harm=%d %s OR=%.3f (%.3f-%.3f) p=%.3g\n",
      tg_, nkept, r$method, exp(r$b), exp(r$b - 1.96*r$se), exp(r$b + 1.96*r$se), r$p))
  flush.console()
  Sys.sleep(0.5)
}
write.csv(res, file.path(P1, "druggable_mr_results.csv"), row.names = FALSE)
print(res)
cat("\n==== ALL DONE ====\n")
sink()
