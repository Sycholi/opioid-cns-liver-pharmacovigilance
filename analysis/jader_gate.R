# jader_gate.R —— JADER(PMDA 202608) 阿片×意识抑制门控 + 肝病人群分层
# 表: demo/drug/reac/hist, CP932编码, 識別番号+報告回数去重(保留最新報告回数)
library(data.table)
setDTthreads(8)
zip <- path.expand("~/nurse_faers/pmda_case_report_202608.zip")
outdir <- path.expand("~/nurse_faers/jader_gate")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_j <- function(file, select_cols) {
  cmd <- sprintf("unzip -p %s %s | iconv -f CP932 -t UTF-8", shQuote(zip), file)
  fread(cmd = cmd, encoding = "UTF-8", select = select_cols, showProgress = FALSE)
}
# ---- 1. 病例表去重: 每識別番号保留最大報告回数 ----
demo <- read_j("demo202608.csv", c("識別番号","報告回数"))
demo[, rk := suppressWarnings(as.integer(報告回数))]
setorder(demo, 識別番号, rk)
keep <- demo[, .SD[.N], by = 識別番号][, .(識別番号, 報告回数)]
n_all <- nrow(keep)
message("JADER 总病例(去重): ", n_all)
kd <- function(dt) dt[keep, on = c("識別番号","報告回数"), nomatch = 0]

# ---- 2. 药物表: 被疑薬(相当于PS)中的阿片 ----
opi_rx <- "(モルヒネ|フェンタニル|オキシコドン|トラマドール|ヒドロモルフォン|コデイン|オキシコンティン|デュロテップ|ブプレノルフィン)"
drug <- read_j("drug202608.csv", c("識別番号","報告回数","医薬品の関与","医薬品（一般名）","医薬品（販売名）"))
drug <- kd(drug)
sus <- drug[医薬品の関与 == "被疑薬"]
opi_ids <- sus[grepl(opi_rx, `医薬品（一般名）`) | grepl(opi_rx, `医薬品（販売名）`), unique(識別番号)]
message("被疑药报告: ", uniqueN(sus$識別番号), "; 阿片被疑药病例: ", length(opi_ids))

# ---- 3. 副作用表: 意识抑制PT集 ----
pt_cns <- c("肝性脳症","傾眠","意識障害","錯乱状態","昏睡")
reac <- read_j("reac202608.csv", c("識別番号","報告回数","有害事象"))
reac <- kd(reac)
ev_he  <- reac[有害事象 == "肝性脳症", unique(識別番号)]
ev_cns <- reac[有害事象 %chin% pt_cns, unique(識別番号)]
message("肝性脳症病例: ", length(ev_he), "; CNS5病例: ", length(ev_cns))

# ---- 4. 原疾患: 肝病人群 ----
liv_rx <- "(肝細胞がん|肝細胞癌|肝癌|肝がん|肝硬変|肝不全)"
hist <- read_j("hist202608.csv", c("識別番号","報告回数","原疾患等"))
hist <- kd(hist)
liver_ids <- hist[grepl(liv_rx, 原疾患等), unique(識別番号)]
message("肝病原疾患病例: ", length(liver_ids))

# ---- 5. ROR: 全库与肝病层内 ----
mk <- function(ev, expv, n) {
  a <- length(intersect(expv, ev)); c <- length(setdiff(ev, expv))
  b <- length(setdiff(expv, ev));   d <- n - a - b - c
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  ror <- (a*d)/(b*c); se <- sqrt(1/a + 1/b + 1/c + 1/d)
  data.table(a=a, b=b, c=c, d=d, ROR=ror,
             ROR_lo=exp(log(ror)-1.96*se), ROR_hi=exp(log(ror)+1.96*se))
}
tab <- rbindlist(list(
  data.table(layer="全库", outcome="CNS5",  mk(ev_cns, opi_ids, n_all)),
  data.table(layer="全库", outcome="肝性脳症", mk(ev_he,  opi_ids, n_all)),
  data.table(layer="肝病层内", outcome="CNS5", mk(ev_cns[ev_cns %in% liver_ids],
        opi_ids[opi_ids %in% liver_ids], length(liver_ids))),
  data.table(layer="肝病层内", outcome="肝性脳症", mk(ev_he[ev_he %in% liver_ids],
        opi_ids[opi_ids %in% liver_ids], length(liver_ids)))))
sub <- data.table(item = c("总病例","阿片被疑药","CNS5","肝性脳症","肝病人群",
                           "肝病×阿片","肝病×阿片×CNS5","肝病×阿片×肝性脳症"),
                  n = c(n_all, length(opi_ids), length(ev_cns), length(ev_he),
                        length(liver_ids), length(intersect(opi_ids, liver_ids)),
                        length(intersect(intersect(opi_ids, liver_ids), ev_cns)),
                        length(intersect(intersect(opi_ids, liver_ids), ev_he))))
print(tab); print(sub)
fwrite(tab, file.path(outdir, "jader_ror.csv"))
fwrite(sub, file.path(outdir, "jader_counts.csv"))
saveRDS(list(opi_ids=opi_ids, ev_he=ev_he, ev_cns=ev_cns,
             liver_ids=liver_ids), file.path(outdir, "jader_intermediates.rds"))
message("JADER门控完成")
