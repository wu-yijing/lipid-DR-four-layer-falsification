# Layer-1 sensitivity to the estimator: add tool-variable-direction-robust estimators
# (MR-RAPS, MR-Lasso, contamination mixture, weighted mode, simple median) alongside
# the estimators already reported (IVW, MR-Egger, weighted median, MR-PRESSO, leave-one-out).
suppressMessages({library(TwoSampleMR); library(MendelianRandomization)})
harm <- Sys.getenv("LDF_HARM", "E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv")
outdir <- Sys.getenv("LDF_OUT", "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")
d <- read.csv(harm, stringsAsFactors = FALSE)
d <- d[d$mr_keep %in% c(TRUE, "TRUE"), ]
order_exp <- c("TG", "LDL", "HDL", "TC", "ApoA1", "ApoB")
res <- list()
add <- function(exp, method, nsnp, b, se, p) {
  res[[length(res) + 1]] <<- data.frame(exposure = exp, method = method, nsnp = nsnp,
    b = b, se = se, p = p, OR = exp(b), lo = exp(b - 1.959964 * se), hi = exp(b + 1.959964 * se),
    stringsAsFactors = FALSE)
}
for (e in order_exp) {
  s <- d[d$exposure == e, ]
  if (nrow(s) < 3) next
  bx <- s$beta.exposure; by <- s$beta.outcome; sx <- s$se.exposure; sy <- s$se.outcome
  r <- tryCatch(suppressWarnings(mr_raps(b_exp = bx, b_out = by, se_exp = sx, se_out = sy)),
                error = function(err) NULL)
  if (!is.null(r)) add(e, "MR-RAPS", r$nsnp, r$b, r$se, r$pval)
  mr <- mr_input(bx = bx, bxse = sx, by = by, byse = sy)
  l <- tryCatch(MendelianRandomization::mr_lasso(mr), error = function(err) NULL)
  if (!is.null(l)) add(e, "MR-Lasso", length(bx), l@Estimate, l@StdError, l@Pvalue)
  cm <- tryCatch(MendelianRandomization::mr_conmix(mr), error = function(err) NULL)
  if (!is.null(cm)) {
    se_cm <- tryCatch((cm@CIUpper - cm@CILower) / (2 * 1.959964), error = function(err) NA_real_)
    if (is.na(se_cm)) se_cm <- tryCatch(cm@StdError, error = function(err) NA_real_)
    add(e, "Contamination mixture", length(bx), cm@Estimate, se_cm, cm@Pvalue)
  }
  mo <- tryCatch(MendelianRandomization::mr_weighted_mode(mr, distribution = "normal", stderror = "delta"), error = function(err) NULL)
  if (!is.null(mo)) {
    se_mo <- tryCatch(mo@StdError, error = function(err) NA_real_)
    if (is.na(se_mo)) se_mo <- tryCatch((mo@CIUpper - mo@CILower) / (2 * 1.959964), error = function(err) NA_real_)
    add(e, "Weighted mode", length(bx), mo@Estimate, se_mo, mo@Pvalue)
  }
  sm <- tryCatch(MendelianRandomization::mr_simple_median(mr, distribution = "normal", stderror = "delta"), error = function(err) NULL)
  if (!is.null(sm)) {
    se_sm <- tryCatch(sm@StdError, error = function(err) NA_real_)
    if (is.na(se_sm)) se_sm <- tryCatch((sm@CIUpper - sm@CILower) / (2 * 1.959964), error = function(err) NA_real_)
    add(e, "Simple median", length(bx), sm@Estimate, se_sm, sm@Pvalue)
  }
}
tab <- do.call(rbind, res)
tab <- tab[order(match(tab$exposure, order_exp), tab$method), ]
num <- c("b", "se", "p", "OR", "lo", "hi")
tab[num] <- lapply(tab[num], function(x) round(as.numeric(x), 6))
write.csv(tab, file.path(outdir, "layer1_robust_estimators.csv"), row.names = FALSE)
cat("rows:", nrow(tab), "\n")
print(tab, row.names = FALSE)
