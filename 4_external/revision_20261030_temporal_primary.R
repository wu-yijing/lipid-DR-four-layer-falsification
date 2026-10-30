# -*- coding: utf-8 -*-
# Revision round 1 (2026-10-30). Temporal validation under the primary outcome definition
# (ungradable photographs treated as missing), adding the transported institutional score to
# the two models that 4_external/revision_20261030_nhanes_design.R covers.
## Temporal validation under the primary outcome definition (ungradable = missing),
## adding the transported Chinese core score to the two models the previous script covered.
suppressMessages({ library(survey) })
SRC <- "E:/workbuddy/2026-09-04-18-19-16"
OUT <- "E:/workbuddy/_r1_review_20261029/pipeB/out_B"

auc_fun <- function(p, y) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0); if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method = "average"); (sum(r[y == 1]) - n1*(n1+1)/2)/(n1*n0)
}

nh <- read.csv(file.path(SRC, "phase2_ext/NHANES_B.csv"))
nh$tc_mmol <- nh$tc/38.67; nh$hdl_mmol <- nh$hdl/38.67
nh$ldl_mmol <- nh$ldl/38.67; nh$tg_mmol <- nh$tg/88.57
nh$tg_log <- log2(nh$tg_mmol); nh$age10 <- nh$age/10
nh$cycle <- nh[[grep("cycle", names(nh), ignore.case = TRUE, value = TRUE)[1]]]
nh$gradable <- !is.na(nh$dr_level4)

dev <- read.csv(file.path(SRC, "phase2_dev/dev_set_v2.csv"))
d <- dev[dev$gradeable == 1 & complete.cases(dev[, c("age","male","tg","tc","hdl_c","ldl_c","apo_a1","apo_b","any_dr")]), ]
d$tg_log <- log2(d$tg); d$age10 <- d$age/10
m_core <- glm(any_dr ~ age10 + male + tc + hdl_c, data = d, family = binomial())
cf <- coef(m_core)

transport <- function(dat) {
  plogis(cf[1] + cf["age10"]*dat$age10 + cf["male"]*dat$male + cf["tc"]*dat$tc_mmol + cf["hdl_c"]*dat$hdl_mmol)
}

fml_full <- any_dr ~ age10 + male + hba1c_pct + dm_duration + egfr_ckdepi2021 + bmi +
            smoker_100 + insulin_use + oral_med + sys_bp
fml_lip  <- any_dr ~ age10 + male + tc_mmol + hdl_mmol

temporal <- function(dat, tag) {
  rs <- list()
  for (pair in list(c("2005-2006","2007-2008"), c("2007-2008","2005-2006"))) {
    tr <- dat[dat$cycle == pair[1], ]; te <- dat[dat$cycle == pair[2], ]
    f1 <- suppressWarnings(glm(fml_full, data = tr, family = binomial()))
    f2 <- suppressWarnings(glm(fml_lip,  data = tr, family = binomial()))
    pf <- plogis(predict(f1, newdata = te, type = "link"))
    pl <- plogis(predict(f2, newdata = te, type = "link"))
    pt <- transport(te)
    rs[[length(rs)+1]] <- data.frame(frame = tag, train = pair[1], test = pair[2],
      n_test = nrow(te), events = sum(te$any_dr),
      auc_full = auc_fun(pf, te$any_dr), auc_lipid = auc_fun(pl, te$any_dr),
      auc_transported = auc_fun(pt, te$any_dr), row.names = NULL)
  }
  do.call(rbind, rs)
}

tt <- rbind(temporal(nh[nh$diabetes == 1, ], "A_submitted"),
            temporal(nh[nh$diabetes == 1 & nh$gradable, ], "B_primary"))
print(tt, digits = 4)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(tt, file.path(OUT, "temporal_three_scores.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "temporal_three_scores.csv"), "\n")
