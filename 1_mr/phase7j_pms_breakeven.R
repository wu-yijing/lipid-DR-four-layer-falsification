## Locate the break-even C-statistic at which the required sample size equals the
## observed n = 1,966 (gradeable/MI cohort), for S = 8 and S = 9 parameters.
suppressMessages(library(pmsampsize))

p <- 0.1398779
s <- 0.90
n_obs <- 1966

cat(sprintf("%-6s %-8s %-12s %-10s %s\n", "S", "C", "SampSize", "Events", "meets n=1966?"))
for (S in c(8, 9)) {
  for (C in seq(0.58, 0.68, by = 0.01)) {
    out <- tryCatch(
      pmsampsize(type = "b", cstatistic = C, parameters = S,
                 prevalence = p, shrinkage = s),
      error = function(e) NULL
    )
    if (is.null(out)) next
    ss <- out$sample_size
    cat(sprintf("%-6d %-8.2f %-12d %-10d %s\n", S, C, ss, ceiling(ss * p),
                ifelse(ss <= n_obs, "YES", "NO")))
  }
  cat("\n")
}
cat("==== ALL DONE (break-even) ====\n")
