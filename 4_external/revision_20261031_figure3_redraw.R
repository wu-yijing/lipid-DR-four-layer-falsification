## Figure 3 re-drawn for the primary analysis (ungradable photographs = missing outcome)
## Same frozen inputs and score definitions as 06_table4_temporal.R (pipeB), so that the
## panels, Table 4 and the legend all describe one prediction set.
## Outputs: Figure3_new.tiff / .png plus a CSV of every plotted point.

SRC <- "E:/workbuddy/2026-09-04-18-19-16"
OUT <- "E:/workbuddy/_r1_review_20261029/pipeB/out_B/fig3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setlocale("LC_ALL", "English")

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average"); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
calib_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y) & p > 1e-8 & p < 1 - 1e-8
  f <- suppressWarnings(glm(y[ok] ~ qlogis(p[ok]), family = binomial()))
  c(intercept = unname(coef(f)[1]), slope = unname(coef(f)[2]))
}
deciles <- function(p, y, g = 10) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  br <- unique(as.numeric(quantile(p, probs = seq(0, 1, length.out = g + 1))))
  grp <- cut(p, breaks = br, include.lowest = TRUE, labels = FALSE)
  data.frame(bin = seq_along(tapply(p, grp, mean)),
             pred_mean = as.numeric(tapply(p, grp, mean)),
             obs_rate  = as.numeric(tapply(y, grp, mean)),
             n         = as.numeric(tapply(y, grp, length)))
}

nh <- read.csv(file.path(SRC, "phase2_ext/NHANES_B.csv"))
nh$tc_mmol <- nh$tc / 38.67; nh$hdl_mmol <- nh$hdl / 38.67
nh$ldl_mmol <- nh$ldl / 38.67; nh$tg_mmol <- nh$tg / 88.57
nh$apob_gl <- nh$apob / 100; nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age / 10
nh$gradable <- !is.na(nh$dr_level4)

dev <- read.csv(file.path(SRC, "phase2_dev/dev_set_v2.csv"))
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age / 10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
m_std  <- glm(any_dr ~ age10 + male + tc + hdl_c + tg_log + ldl_c, data = d, family = binomial())
L <- read.csv(file.path(SRC, "phase2_out_v2/phase2_lasso.csv"))
lc <- setNames(L$coef, L$var)
allt <- c("(Intercept)","age10","male","tg_log","tc","hdl_c","ldl_c","apo_a1","apo_b")
miss <- setdiff(allt, names(lc)); if (length(miss)) lc[miss] <- 0

## primary frame: diabetes, ascertainable retinal status (same as Table 4, scenario B)
dm <- nh[nh$diabetes == 1 & nh$gradable, ]
cf <- coef(m_core); cf2 <- coef(m_std)
score <- list(
  core = plogis(cf[1] + cf["age10"] * dm$age10 + cf["male"] * dm$male +
                cf["tc"] * dm$tc_mmol + cf["hdl_c"] * dm$hdl_mmol),
  std  = plogis(cf2[1] + cf2["age10"] * dm$age10 + cf2["male"] * dm$male +
                cf2["tc"] * dm$tc_mmol + cf2["hdl_c"] * dm$hdl_mmol +
                cf2["tg_log"] * dm$tg_log + cf2["ldl_c"] * dm$ldl_mmol),
  lasso = plogis(lc["(Intercept)"] + lc["age10"] * dm$age10 + lc["male"] * dm$male +
                 lc["tg_log"] * dm$tg_log + lc["tc"] * dm$tc_mmol +
                 lc["hdl_c"] * dm$hdl_mmol + lc["ldl_c"] * dm$ldl_mmol +
                 lc["apo_b"] * dm$apob_gl)
)

res <- list()
for (nm in names(score)) {
  p <- score[[nm]]; y <- dm$any_dr
  ok <- !is.na(p)
  ct <- deciles(p, y)
  ca <- calib_fun(p, y)
  res[[nm]] <- list(p = p, y = y, ok = ok, ct = ct, auc = auc_fun(p, y), slope = ca["slope"],
                    intercept = ca["intercept"],
                    mean_pred = mean(p[ok]), mean_obs = mean(y[ok]),
                    dev_lo = min(ct$obs_rate - ct$pred_mean), dev_hi = max(ct$obs_rate - ct$pred_mean),
                    n = sum(ok), events = sum(y[ok]))
  cat(sprintf("%-5s n=%3d ev=%3d AUC=%.4f slope=%+.3f int=%+.3f meanpred=%.4f meanobs=%.4f dev=[%+.3f, %+.3f]\n",
              nm, sum(ok), sum(y[ok]), res[[nm]]$auc, res[[nm]]$slope, res[[nm]]$intercept,
              res[[nm]]$mean_pred, res[[nm]]$mean_obs, res[[nm]]$dev_lo, res[[nm]]$dev_hi))
}

pts <- do.call(rbind, lapply(names(res), function(nm)
  cbind(version = nm, res[[nm]]$ct)))
write.csv(pts, file.path(OUT, "figure3_deciles.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(names(res), function(nm) data.frame(
  version = nm, n = res[[nm]]$n, events = res[[nm]]$events, auc = res[[nm]]$auc,
  slope = res[[nm]]$slope, intercept = res[[nm]]$intercept,
  mean_pred = res[[nm]]$mean_pred, mean_obs = res[[nm]]$mean_obs,
  dev_min = res[[nm]]$dev_lo, dev_max = res[[nm]]$dev_hi)))
write.csv(summ, file.path(OUT, "figure3_summary.csv"), row.names = FALSE)

## ---------------- draw ----------------
titles <- c(core = "core (age, sex, TC, HDL-C)", std = "std (+TG, LDL-C)", lasso = "LASSO (+ApoB)")
tags   <- c(core = "(a) core", std = "(b) std", lasso = "(c) LASSO")
cols   <- c(core = "#2C7FB8", std = "#E69F00", lasso = "#D55E00")
lim <- ceiling(max(unlist(lapply(res, function(z) c(z$ct$pred_mean, z$ct$obs_rate)))) * 1.1 * 20) / 20
lim <- max(lim, 0.65); lim <- 0.65

draw <- function(file, type) {
  if (type == "tiff") tiff(file, width = 6.09, height = 3.15, units = "in", res = 300,
                           compression = "lzw", family = "sans", pointsize = 9)
  else png(file, width = 6.09, height = 3.15, units = "in", res = 300, family = "sans",
           type = "cairo", pointsize = 9)
  par(mfrow = c(1, 3), mar = c(4.5, 3.7, 1.9, 0.6), mgp = c(2.0, 0.5, 0), las = 1,
      cex.axis = 0.80, cex.lab = 0.90, tcl = -0.25, bty = "n")
  for (nm in c("core", "std", "lasso")) {
    z <- res[[nm]]; ct <- z$ct
    plot(NA, xlim = c(0, lim), ylim = c(0, lim), xlab = "", ylab = "", axes = FALSE)
    abline(0, 1, lty = 2, col = "grey55", lwd = 1.1)
    lines(ct$pred_mean, ct$obs_rate, col = cols[nm], lwd = 1.4)
    points(ct$pred_mean, ct$obs_rate, pch = 19, cex = 0.80, col = cols[nm])
    axis(1, at = c(seq(0, 0.6, 0.2), lim), labels = c("0.0","0.2","0.4","0.6",""),
         col = "grey25", col.axis = "black", lwd = 0.9)
    axis(2, at = c(seq(0, 0.6, 0.1), lim), labels = c(sprintf("%.1f", seq(0, 0.6, 0.1)), ""),
         col = "grey25", col.axis = "black", lwd = 0.9)
    title(main = titles[nm], cex.main = 0.90, font.main = 1, line = 0.5)
    title(xlab = "Predicted risk", cex.lab = 0.88, line = 1.75)
    if (nm == "core") title(ylab = "Observed any-DR rate", cex.lab = 0.88, line = 2.2)
    text(lim * 0.97, lim * 0.03, adj = c(1, 0), cex = 0.70, col = "black",
         labels = sprintf("AUC %.3f\nslope %+.2f\nintercept %+.2f", z$auc, z$slope, z$intercept))
    mtext(tags[nm], side = 1, line = 3.15, adj = 0.5, font = 2, cex = 0.92)
  }
  dev.off()
}
draw(file.path(OUT, "Figure3_new.tiff"), "tiff")
draw(file.path(OUT, "Figure3_new.png"), "png")
cat("written Figure3_new.tiff / .png\n")
