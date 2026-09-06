# Figure 2: MR forest plot - univariable (Panel A) vs multivariable (Panel B)
# Publication style: 300 dpi, sans-serif, colorblind-friendly
png("E:/workbuddy/2026-09-04-18-19-16/phase6_out/fig2_mr_forest.png",
    width = 2600, height = 2100, res = 300)
par(mfrow = c(2, 1), mar = c(4.5, 10, 2.5, 14), family = "sans")

BLUE <- "#2C7FB8"; RED <- "#D7301F"; GREY <- "#636363"

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 1))
  formatC(p, format = "f", digits = 3)
}

## ---- Panel A: univariable IVW ----
main <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase1_out/mr_main.csv")
ivw <- main[grepl("Inverse variance", main$method, ignore.case = TRUE),
            c("exposure", "or", "or_lci95", "or_uci95", "pval")]
# fixed order, key lipid on top
order_a <- c("TG", "LDL", "HDL", "TC", "ApoA1", "ApoB")
lab_a <- c(TG = "TG", LDL = "LDL-C", HDL = "HDL-C", TC = "Total cholesterol",
           ApoA1 = "ApoA-I", ApoB = "ApoB")
ivw <- ivw[match(order_a, ivw$exposure), ]
ivw$lab <- lab_a[ivw$exposure]

xlim <- c(0.55, 1.75)
ylim <- c(0.5, nrow(ivw) + 1.2)
plot(NA, NA, xlim = xlim, ylim = ylim, log = "x",
     xlab = "OR for diabetic retinopathy (per 1-SD genetically predicted lipid)",
     ylab = "", yaxt = "n", main = "A  Univariable MR (IVW)", cex.main = 1.05,
     cex.lab = 0.9)
abline(v = 1, lty = 2, col = GREY, lwd = 1)
for (i in seq_len(nrow(ivw))) {
  y <- nrow(ivw) - i + 1
  sig <- ivw$pval[i] < 0.05
  segments(ivw$or_lci95[i], y, ivw$or_uci95[i], y, lwd = 2,
           col = if (sig) RED else BLUE)
  points(ivw$or[i], y, pch = if (sig) 16 else 21, cex = 1.15,
         col = if (sig) RED else BLUE,
         bg = if (sig) RED else "white")
  text(1.78, y, pos = 4, cex = 0.82, family = "sans",
       labels = sprintf("%.2f (%.2f-%.2f)  p=%s",
                        ivw$or[i], ivw$or_lci95[i], ivw$or_uci95[i], fmt_p(ivw$pval[i])),
       col = if (sig) RED else "black", xpd = TRUE)
  text(0.53, y, pos = 2, cex = 0.95, labels = ivw$lab[i], xpd = TRUE, family = "sans")
}

## ---- Panel B: multivariable MR (core: TG + LDL + HDL mutually adjusted) ----
mv <- read.csv("E:/workbuddy/2026-09-04-18-19-16/phase1_out/mvmr_core_TG_LDL_HDL.csv")
order_b <- c("TG", "LDL", "HDL")
mv$short <- trimws(sub(" \\|\\|.*", "", mv$exposure))   # strip " || id:..." suffix
full_b <- c(TG = "Triglycerides", LDL = "LDL cholesterol", HDL = "HDL cholesterol")
mv <- mv[match(unname(full_b), mv$short), ]
mv$exposure <- names(full_b)
mv$lab <- lab_a[mv$exposure]
stopifnot(!any(is.na(mv$pval)))

plot(NA, NA, xlim = xlim, ylim = c(0.5, length(order_a) + 1.2), log = "x",
     xlab = "OR for diabetic retinopathy (per 1-SD, mutually adjusted)",
     ylab = "", yaxt = "n", main = "B  Multivariable MR (TG + LDL-C + HDL-C)",
     cex.main = 1.05, cex.lab = 0.9)
abline(v = 1, lty = 2, col = GREY, lwd = 1)
for (i in seq_len(nrow(mv))) {
  y <- length(order_a) - i + 1   # align rows with Panel A
  sig <- mv$pval[i] < 0.05
  segments(mv$or_lci95[i], y, mv$or_uci95[i], y, lwd = 2,
           col = if (sig) RED else BLUE)
  points(mv$or[i], y, pch = if (sig) 16 else 21, cex = 1.15,
         col = if (sig) RED else BLUE, bg = if (sig) RED else "white")
  text(1.78, y, pos = 4, cex = 0.82, family = "sans",
       labels = sprintf("%.2f (%.2f-%.2f)  p=%s",
                        mv$or[i], mv$or_lci95[i], mv$or_uci95[i], fmt_p(mv$pval[i])),
       col = "black", xpd = TRUE)
  text(0.53, y, pos = 2, cex = 0.95, labels = mv$lab[i], xpd = TRUE, family = "sans")
}
# annotate attenuated TG (annotation above the TG row, inside plot area)
tg_row <- mv[mv$exposure == "TG", ]
y_tg <- length(order_a) - match("TG", order_b) + 1
arrows(1.28, y_tg + 0.40, 1.02, y_tg + 0.10, length = 0.08, col = RED, lwd = 1.2)
text(1.29, y_tg + 0.48, "signal attenuates to null after adjustment",
     cex = 0.78, col = RED, adj = c(0, 0), xpd = TRUE, family = "sans")

dev.off()
cat("fig2_mr_forest.png written\n")

# also export the plotted numbers for the caption/tables
write.csv(ivw, "E:/workbuddy/2026-09-04-18-19-16/phase6_out/fig2_panelA_univariable.csv", row.names = FALSE)
write.csv(mv,  "E:/workbuddy/2026-09-04-18-19-16/phase6_out/fig2_panelB_mvmr.csv", row.names = FALSE)
cat("panel CSVs written\n")
