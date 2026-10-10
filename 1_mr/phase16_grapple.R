# GRAPPLE (Genome-wide mR Analysis under Pervasive PLEiotropy): robust estimation of the
# causal effect under a random-effects model of pleiotropic effects, per exposure.
# Uses the already harmonised instrument-level input, so no PLINK/clumping is needed.
harm <- Sys.getenv("LDF_HARM", "E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv")
outdir <- Sys.getenv("LDF_OUT", "E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")
ok <- suppressMessages(requireNamespace("grapple", quietly = TRUE)) ||
      suppressMessages(requireNamespace("GRAPPLE", quietly = TRUE))
if (!ok) { cat("GRAPPLE not installed\n"); quit(status = 0) }
ns <- if (requireNamespace("grapple", quietly=TRUE)) "grapple" else "GRAPPLE"
d <- read.csv(harm, stringsAsFactors = FALSE)
d <- d[d$mr_keep %in% c(TRUE, "TRUE"), ]
order_exp <- c("TG", "LDL", "HDL", "TC", "ApoA1", "ApoB")
res <- list()
f <- get("grappleRobustEst", envir = asNamespace(ns))
for (e in order_exp) {
  s <- d[d$exposure == e, ]
  dat <- data.frame(SNP = s$SNP, gamma_exp1 = s$beta.exposure, se_exp1 = s$se.exposure,
                    gamma_out1 = s$beta.outcome, se_out1 = s$se.outcome, stringsAsFactors = FALSE)
  r <- tryCatch(f(dat, p.thres = NULL, loss.function = "tukey"), error = function(err) { cat(e, "err:", conditionMessage(err), "\n"); NULL })
  if (is.null(r)) next
  res[[length(res) + 1]] <- data.frame(exposure = e, method = "GRAPPLE (Tukey, random effects)",
    nsnp = nrow(dat), b = as.numeric(r$beta.hat)[1], se = sqrt(as.numeric(r$beta.variance)[1, 1]),
    p = as.numeric(r$beta.p.vaue)[1], tau2 = as.numeric(r$tau2.hat)[1], stringsAsFactors = FALSE)
}
if (length(res)) {
  tab <- do.call(rbind, res)
  z <- 1.959963984540054
  tab$OR <- exp(tab$b); tab$lo <- exp(tab$b - z * tab$se); tab$hi <- exp(tab$b + z * tab$se)
  write.csv(tab, file.path(outdir, "layer1_grapple.csv"), row.names = FALSE)
  print(tab, row.names = FALSE)
} else cat("no GRAPPLE results\n")
