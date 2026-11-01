## =====================================================================
## revision_20261031: Supplementary Table S1 under the primary outcome
## definition of Methods 2.2 (ungradable photographs = missing outcome).
##
## Why: Table S1 as published in the revision was still computed on the
## 1,119-participant diabetes subgroup with ungradable photographs counted
## as non-cases (the original definition). Every other layer-3 table
## (S19, S21, S31, S41, S46) now uses the primary frame (884 with
## ascertainable retinal status), so S1 was the only table on the old
## footing.
##
## Method is identical to phase4c_rc_tyg_v2.R (same metric definitions,
## same complete-case filter, same rank-based unweighted AUC); only the
## row selection changes. Variant A is reproduced and asserted against
## the published S1 values before variant C is read.
##
## Run: Rscript 4_external/revision_20261031_tableS1_primary_frame.R
## =====================================================================

SRC <- Sys.getenv("LDF_SRC", "E:/workbuddy/2026-09-04-18-19-16")
this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- if (is.na(this)) "." else normalizePath(file.path(dirname(this), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "output_tables"))) REPO <- "."
OUT <- file.path(REPO, "output_tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

nh0 <- read.csv(file.path(SRC, "phase2_ext/NHANES_B.csv"))
ungradable <- is.na(nh0$dr_level4)

prep <- function(nh) {
  nh$tc_mmol  <- nh$tc  / 38.67
  nh$hdl_mmol <- nh$hdl / 38.67
  nh$ldl_mmol <- nh$ldl / 38.67
  nh$tg_mmol  <- nh$tg  / 88.57
  nh$rc_mmol  <- nh$tc_mmol - nh$ldl_mmol - nh$hdl_mmol        # RC = TC - LDL - HDL
  nh$tyg      <- log(nh$tg * nh$glucose / 2)                    # ln[TG(mg/dL) x FPG(mg/dL)/2]
  nh$tyg_bmi  <- nh$tyg * nh$bmi
  nh$age10    <- nh$age / 10
  nh
}

mets <- c("tc_mmol", "hdl_mmol", "ldl_mmol", "tg_mmol", "rc_mmol", "tyg", "tyg_bmi")
mlab <- c("TC", "HDL-C", "LDL-C", "TG", "Remnant cholesterol", "TyG index", "TyG x BMI")

run <- function(variant) {
  nh <- nh0
  if (variant == "B") nh$any_dr[ungradable] <- 1
  if (variant == "C") nh <- nh[!ungradable, ]
  nh <- prep(nh)
  dm <- nh[nh$diabetes == 1, ]
  res <- data.frame()
  for (i in seq_along(mets)) {
    x <- mets[i]
    dx <- dm[!is.na(dm[[x]]) & !is.na(dm$any_dr) & !is.na(dm$age10) & !is.na(dm$male), ]
    res <- rbind(res, data.frame(variant = variant, metric = mlab[i],
                                 n = nrow(dx), events = sum(dx$any_dr),
                                 auc = auc_fun(dx[[x]], dx$any_dr)))
  }
  attr(res, "subgroup") <- nrow(dm)
  res
}

A <- run("A"); B <- run("B"); C <- run("C")
cat(sprintf("diabetes subgroup: A=%d  B=%d  C=%d (retinal status ascertainable)\n",
            attr(A, "subgroup"), attr(B, "subgroup"), attr(C, "subgroup")))

## ---- assert the published variant A (Table S1, old footing) ----------
published_n   <- c(1041, 1041, 494, 525, 494, 523, 502)
published_auc <- c(0.501, 0.518, 0.495, 0.491, 0.467, 0.543, 0.474)
ok <- all(A$n == published_n) && max(abs(A$auc - published_auc)) < 0.0005
cat("variant A reproduces published S1:", if (ok) "PASS" else "FAIL", "\n")
if (!ok) {
  print(cbind(A[, c("metric","n","auc")], published_n, published_auc))
  stop("baseline not reproduced - do not use variant C")
}

print(C, row.names = FALSE)
cat("\nChange from A to C (n, AUC):\n")
print(data.frame(metric = C$metric, n_A = A$n, n_C = C$n,
                 d_n = C$n - A$n, auc_A = round(A$auc, 3),
                 auc_C = round(C$auc, 3), d_auc = round(C$auc - A$auc, 3)),
      row.names = FALSE)

out <- rbind(A, B, C)
out$auc <- round(out$auc, 4)
write.csv(out, file.path(OUT, "tableS1_primary_frame.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT, "tableS1_primary_frame.csv"), "\n")
