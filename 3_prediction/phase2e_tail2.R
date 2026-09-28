## phase2e_tail2.R -- completes what phase2e_redo2.R could not.
## (a) three-row AUC decomposition (values from the completed run in
##     phase2_out/phase2e_redo2.log; provenance stated inline)
## (b) pmsampsize grid with a corrected call (pmsampsize 1.1.3 accepts
##     cstatistic OR csquared, not both)
## Pure ASCII, no sink() (Rscript + sink(split=TRUE) aborted on this host).
OUT <- "E:/workbuddy/2026-09-04-18-19-16/phase2_out"
suppressMessages(library(pmsampsize))

dec <- data.frame(
  row = c("1_CCcoef_to_CCsubset_apparent", "2_CCcoef_to_CCsubset_corrected",
          "3_CCcoef_to_MIcohort", "4_MIcoef_to_MIcohort_apparent",
          "5_MIcoef_to_MIcohort_corrected", "6_IPW_CC_subset"),
  AUC = c(0.6923, 0.6678, 0.5666, 0.6266, 0.6083, 0.6920),
  source = c(rep("phase2e_redo2.log", 5), "phase2e_ipw.log"))
write.csv(dec, file.path(OUT, "phase2e_decomposition.csv"), row.names = FALSE)
cat("--- decomposition (provenance: phase2e_redo2.log / phase2e_ipw.log) ---\n")
print(dec, row.names = FALSE)

p <- 0.1398779; s <- 0.90
cat(sprintf("\npmsampsize %s | prevalence=%.7f | shrinkage=%.2f\n",
            as.character(packageVersion("pmsampsize")), p, s))
grid <- expand.grid(S = c(8, 9), C = c(0.67, 0.66, 0.65, 0.61, 0.60, 0.55))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  S <- grid$S[i]; C <- grid$C[i]
  r <- tryCatch(pmsampsize(type = "b", cstatistic = C, parameters = S,
                           prevalence = p, shrinkage = s), error = function(e) NULL)
  if (is.null(r)) return(data.frame(S = S, C = C, n_required = NA, events_required = NA))
  data.frame(S = S, C = C, n_required = r$sample_size,
             events_required = ceiling(r$sample_size * p))
}))
res$meets_n1966 <- res$n_required <= 1966
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "phase2e_pmsampsize_grid.csv"), row.names = FALSE)

be <- do.call(rbind, lapply(seq(0.55, 0.70, by = 0.005), function(C) {
  r <- tryCatch(pmsampsize(type = "b", cstatistic = C, parameters = 8,
                           prevalence = p, shrinkage = s), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(C = C, n_required = r$sample_size, meets_n1966 = r$sample_size <= 1966)
}))
write.csv(be, file.path(OUT, "phase2e_pmsampsize_breakeven.csv"), row.names = FALSE)
cat(sprintf("\nbreak-even (S=8, n=1966): highest C still requiring <=1966 = %.3f\n",
            max(be$C[be$meets_n1966])))

nh <- pmsampsize(type = "b", cstatistic = 0.733, parameters = 10,
                 prevalence = 0.3006757, shrinkage = s)
cat(sprintf("NHANES full model (C=0.733, S=10, prevalence=0.3007): n=%d events=%d\n",
            nh$sample_size, ceiling(nh$sample_size * 0.3006757)))
write.csv(data.frame(
  model = c("chinese_lipid_S8_C0.67", "nhanes_full_S10_C0.733"),
  n_required = c(res$n_required[res$S == 8 & res$C == 0.67], nh$sample_size),
  events_required = c(res$events_required[res$S == 8 & res$C == 0.67],
                      ceiling(nh$sample_size * 0.3006757))),
  file.path(OUT, "phase2e_pmsampsize_headline.csv"), row.names = FALSE)
cat("\n==== ALL DONE ====\n")
