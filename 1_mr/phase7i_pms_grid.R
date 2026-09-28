## Verify pmsampsize across parameter-count and C-statistic assumptions.
## pmsampsize 1.1.3 signature uses `cstatistic` (not `csquared`).
suppressMessages(library(pmsampsize))

p <- 0.1398779   # gradeable-cohort prevalence
s <- 0.90        # target shrinkage

cat("pmsampsize version:", as.character(packageVersion("pmsampsize")), "\n")
cat("prevalence =", p, " shrinkage =", s, "\n\n")

cat(sprintf("%-6s %-8s %-10s %-10s\n", "S", "C", "SampSize", "Events"))
for (S in c(8, 9)) {
  for (C in c(0.67, 0.60, 0.55)) {
    out <- tryCatch(
      pmsampsize(type = "b", cstatistic = C, parameters = S,
                 prevalence = p, shrinkage = s),
      error = function(e) e
    )
    if (inherits(out, "error")) {
      cat(sprintf("%-6d %-8.2f %-10s %-10s\n", S, C, "ERROR", "-"))
    } else {
      ss <- out$sample_size
      cat(sprintf("%-6d %-8.2f %-10d %-10d\n", S, C, ss, ceiling(ss * p)))
    }
  }
}

cat("\n-- reference: NHANES full model (prevalence 0.3006757, C=0.735, S=10) --\n")
o <- pmsampsize(type = "b", cstatistic = 0.735, parameters = 10,
                prevalence = 0.3006757, shrinkage = s)
cat("SampSize =", o$sample_size, " Events =", ceiling(o$sample_size * 0.3006757), "\n")

cat("\n-- observed in this study --\n")
cat("gradeable/MI cohort : n = 1966, events = 275\n")
cat("complete cases      : n = 1207, events = 121\n")
cat("\n==== ALL DONE (pms grid) ====\n")
