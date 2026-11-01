## Layer-2 calibration: apparent and optimism-corrected slope and intercept
##
## STATUS: NOT RUN in this repository. It needs the individual-level development
## set (phase2_dev/dev_set_v2_dedup.csv), which is not deposited here: the
## repository holds aggregate tables only. The manuscript therefore reports that
## calibration of the layer-2 model is not available (Discussion 4.5) rather than
## quoting the apparent slope of 1.00 and intercept of 0.00, which are identities
## when a model is fitted and evaluated on the same records.
##
## Run by the data holder as:
##   Rscript 3_prediction/revision_20261102_layer2_calibration.R
##   LDF_DEV=/path/to/dev_set_v2_dedup.csv LDF_OUT=/path/to/out \
##     Rscript 3_prediction/revision_20261102_layer2_calibration.R
##
## Model, outcome and complete-case definition follow phase2_model_v2_dedup.R so
## that the calibration refers to the coefficients of Table 3.

suppressMessages(library(stats))

dev_path <- Sys.getenv("LDF_DEV",
  unset = "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2_dedup.csv")
out_dir <- Sys.getenv("LDF_OUT", unset = "3_prediction/_out_layer2_calibration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261102)
B <- as.integer(Sys.getenv("LDF_B", unset = "1000"))

cat("==== load ====\n")
d <- read.csv(dev_path, stringsAsFactors = FALSE)
vars <- c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
sub <- d[d$gradeable == 1 & complete.cases(d[, c(vars, "any_dr")]), ]
sub$tg_log <- log2(sub$tg)
sub$age10  <- sub$age / 10
f <- any_dr ~ age10 + male + tg_log + tc + hdl_c + ldl_c + apo_a1 + apo_b
cat(sprintf("complete cases n=%d  events=%d (%.1f%%)\n",
            nrow(sub), sum(sub$any_dr), 100 * mean(sub$any_dr)))

## --- calibration summaries, defined as in Van Calster et al. ---
## slope     : coefficient of the linear predictor in glm(y ~ lp)
## intercept : intercept of glm(y ~ offset(lp)), "calibration-in-the-large"
cal <- function(y, lp) {
  s <- try(glm(y ~ lp, family = binomial), silent = TRUE)
  i <- try(glm(y ~ offset(lp), family = binomial), silent = TRUE)
  c(slope     = if (inherits(s, "try-error")) NA_real_ else unname(coef(s)[2]),
    intercept = if (inherits(i, "try-error")) NA_real_ else unname(coef(i)[1]))
}

fit <- glm(f, data = sub, family = binomial)
lp  <- predict(fit, type = "link")
app <- cal(sub$any_dr, lp)
cat(sprintf("apparent   slope=%.3f  intercept=%.3f  (both identities by construction)\n",
            app["slope"], app["intercept"]))

cat(sprintf("\n==== optimism correction, B=%d ====\n", B))
opt <- matrix(NA_real_, nrow = B, ncol = 2,
              dimnames = list(NULL, c("slope", "intercept")))
for (b in seq_len(B)) {
  idx <- sample.int(nrow(sub), replace = TRUE)
  bs  <- sub[idx, ]
  fb  <- try(glm(f, data = bs, family = binomial), silent = TRUE)
  if (inherits(fb, "try-error")) next
  ## performance of the bootstrap model on the bootstrap sample
  cb <- cal(bs$any_dr, predict(fb, newdata = bs, type = "link"))
  ## same model applied to the original sample
  co <- cal(sub$any_dr, predict(fb, newdata = sub, type = "link"))
  opt[b, ] <- cb - co                      # optimism for this replicate
}
ok  <- stats::complete.cases(opt)
corr <- app - colMeans(opt[ok, , drop = FALSE])

res <- data.frame(
  quantity = c("apparent_slope", "apparent_intercept",
               "optimism_slope", "optimism_intercept",
               "corrected_slope", "corrected_intercept",
               "n", "events", "B_requested", "B_used"),
  value = c(app["slope"], app["intercept"],
            mean(opt[ok, "slope"]), mean(opt[ok, "intercept"]),
            corr["slope"], corr["intercept"],
            nrow(sub), sum(sub$any_dr), B, sum(ok)),
  row.names = NULL)

write.csv(res, file.path(out_dir, "layer2_calibration_optimism.csv"), row.names = FALSE)
write.csv(opt[ok, , drop = FALSE], file.path(out_dir, "layer2_calibration_replicates.csv"),
          row.names = FALSE)

cat(sprintf("corrected  slope=%.3f  intercept=%.3f  (B used=%d)\n",
            corr["slope"], corr["intercept"], sum(ok)))
cat("written to:", out_dir, "\n")

## NOTE the multiple-imputation variant is not implemented here: it needs the
## long-format imputed file and Rubin's rules applied to the calibration slope
## itself, which is a non-linear functional of the fitted probabilities.
