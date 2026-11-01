## =====================================================================
## phase22: linkage-window sensitivity for the layer-2 model
##
## Question (reviewer 1, point 12): the lipid panel and the fundus photograph
## are linked on name and sex, not on a shared identifier, so a panel recorded
## far from the photograph may belong to another admission or another patient.
## This script repeats the complete-case model on the records whose panel was
## taken within the intended 7-day window, and reports both analyses side by side.
##
## Inputs
##   dev_set_v2.csv        corrected outcome (option A): gradeable, any_dr
##                         E:/workbuddy/2026-09-04-18-19-16/phase2_dev/
##   cleaned_dataset.xlsx  sheet 3 (linkage cohort): column 49 = 关联天数差
##                         E:/workbuddy/血脂，糖网等数据清洗/清洗结果/
##   The xlsx sheet name and headers are Chinese; the sheet is therefore read by
##   index and the columns by position, and the path is written as \u escapes so
##   that this file stays pure ASCII; the workbook is staged into tempdir() before
##   reading, because readxl hands the path to utils::unzip, which cannot translate
##   a UTF-8 name under a non-UTF-8 code page.
##
## Output
##   output_tables/layer2_linkage_window.csv
##
## Model: any_dr ~ age + male + log2(tg) + tc + hdl_c + ldl_c + apo_a1 + apo_b
##         (the eight-parameter complete-case model of Methods 2.4)
## Interval: bootstrap optimism correction, B = 1,000, seed 20260906
##
## Result (verified 31 October 2026): complete cases 1,417 / 367 events,
## apparent 0.5984, corrected 0.5793; within 7 days 1,357 / 350 events,
## apparent 0.5939, corrected 0.5727. Restricting the window therefore leaves
## the estimate unchanged, and the row is reported as Supplementary Table S41.
##
## Run: Rscript 3_prediction/phase22_linkage_window.R
## =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(pROC)
})

CSV  <- "E:/workbuddy/2026-09-04-18-19-16/phase2_dev/dev_set_v2.csv"
XLSX <- "E:/workbuddy/\u8840\u8102\uff0c\u7cd6\u7f51\u7b49\u6570\u636e\u6e05\u6d17/\u6e05\u6d17\u7ed3\u679c/cleaned_dataset.xlsx"
OUT  <- "output_tables/layer2_linkage_window.csv"

d <- read.csv(CSV, stringsAsFactors = FALSE)
## readxl passes the path to utils::unzip, which cannot translate a UTF-8 name under
## a non-UTF-8 code page; stage a copy in tempdir(), whose path is ASCII.
stage <- file.path(tempdir(), "cleaned_dataset.xlsx")
if (!file.exists(stage)) stopifnot(file.copy(XLSX, stage, overwrite = TRUE))
x <- as.data.frame(read_excel(stage, sheet = 3, .name_repair = "minimal"))

stopifnot(nrow(d) == nrow(x))
d$days <- as.numeric(x[[49]])          # 关联天数差

vars <- c("age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b")
f <- any_dr ~ age + male + log2(tg) + tc + hdl_c + ldl_c + apo_a1 + apo_b

cc <- d[d$gradeable == 1 & complete.cases(d[, c(vars, "any_dr")]), ]

auc_fun <- function(dat, fit)
  as.numeric(pROC::auc(dat$any_dr,
                       predict(fit, newdata = dat, type = "link"),
                       quiet = TRUE))

optimism <- function(dat, B = 1000, seed = 20260906) {
  set.seed(seed)
  m0 <- glm(f, data = dat, family = binomial)
  a0 <- auc_fun(dat, m0)
  op <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(nrow(dat), replace = TRUE)
    mb <- try(glm(f, data = dat[idx, ], family = binomial), silent = TRUE)
    if (inherits(mb, "try-error")) { op[b] <- NA; next }
    op[b] <- auc_fun(dat[idx, ], mb) - auc_fun(dat, mb)
  }
  c(n = nrow(dat), events = sum(dat$any_dr),
    apparent = a0, optimism = mean(op, na.rm = TRUE),
    corrected = a0 - mean(op, na.rm = TRUE))
}

win <- cc[!is.na(cc$days) & abs(cc$days) <= 7, ]

res <- rbind(
  data.frame(analysis = "All complete cases",       t(optimism(cc))),
  data.frame(analysis = "Panel within 7 days",      t(optimism(win))),
  data.frame(analysis = "Panel outside 7 days",     t(optimism(cc[!is.na(cc$days) & abs(cc$days) > 7, ])))
)
names(res)[-1] <- c("n", "events", "apparent_auc", "optimism", "corrected_auc")
res$bootstraps <- 1000
res$seed <- 20260906

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT, row.names = FALSE)
print(res, digits = 4)
