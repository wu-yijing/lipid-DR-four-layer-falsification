# OpenGWAS 端到端 smoke test —— 验证 token 鉴权 + 拉数
# 关键: ieugwasr 1.1.0.9000 的 get_opengwas_jwt() 读取环境变量 OPENGWAS_JWT
#        (不是 OPENGWAS_TOKEN, 也不是 ~/.opengwas_token 文件)

sink("E:/workbuddy/2026-09-04-18-19-16/smoke_test.log", type = "output")

# 1) 读取 token 并设置正确环境变量
tok_path <- "E:/workbuddy/2026-06-23-15-36-17/data/opengwas_token.txt"
tok <- trimws(readLines(tok_path, warn = FALSE)[1])
Sys.setenv(OPENGWAS_JWT = tok)
cat("== token 读取成功, 长度:", nchar(tok), "\n")
cat("== OPENGWAS_JWT 已设置, 长度:", nchar(Sys.getenv("OPENGWAS_JWT")), "\n")
cat("== ieugwasr 默认 API base:", getOption("ieugwasr_api"), "\n")

library(ieugwasr)

# 2) 鉴权测试：user() 需要有效 token, 无效会 401
cat("\n== [1/3] 鉴权测试 ieugwasr::user() ==\n")
usr <- tryCatch(user(), error = function(e) {
  cat("!! user() 失败:", conditionMessage(e), "\n"); NULL
})
if (!is.null(usr)) {
  cat(">> 鉴权成功, 账号信息:\n")
  print(usr)
} else {
  cat(">> 鉴权失败\n")
  sink(); quit(status = 1)
}

# 3) 定位血脂 trait (从元数据自动匹配, 避免硬编码错误 ID)
cat("\n== [2/3] 定位血脂 trait (traits 元数据) ==\n")
tdef <- tryCatch(traits(), error = function(e) {
  cat("!! traits() 失败:", conditionMessage(e), "\n"); NULL
})
lipid_id <- NULL
if (!is.null(tdef) && nrow(tdef) > 0) {
  cats <- tolower(paste(tdef$trait, tdef$id, sep = " "))
  hit <- grep("cholesterol|ldl|hdl|triglyceride|lipid", cats)
  if (length(hit) > 0) {
    lipid_id <- tdef$id[hit[1]]
    cat(">> 命中血脂 trait:", tdef$trait[hit[1]], " (", lipid_id, ")\n")
  } else {
    cat("!! 未匹配, 回退 ieu-a-300\n"); lipid_id <- "ieu-a-300"
  }
} else {
  lipid_id <- "ieu-a-300"
}

# 4) 拉数测试：取该 trait 的 top hits (p < 1e-5)
cat("\n== [3/3] 拉数测试 tophits() ==\n")
th <- tryCatch(tophits(id = lipid_id, pval = 1e-5), error = function(e) {
  cat("!! tophits() 失败:", conditionMessage(e), "\n"); NULL
})
if (!is.null(th) && nrow(th) > 0) {
  cat(">> 拉到", nrow(th), "个显著位点 (p<1e-5)\n")
  cols <- intersect(c("rsid","chromosome","position","beta","se","p","id"), names(th))
  cat(">> 示例前 5 行:\n")
  print(head(th[ , cols], 5))
  cat("\n>>> 端到端 smoke test: 全部通过 (鉴权 OK + 拉数 OK) <<<\n")
} else {
  cat("!! 未拉到数据 (trait id 或端点问题), 但鉴权已通过\n")
}
sink()
