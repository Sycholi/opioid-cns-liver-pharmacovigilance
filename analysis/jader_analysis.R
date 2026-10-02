# jader_analysis.R —— JADER 正式分析(药物集与FAERS对齐: 吗啡/氢吗啡酮/芬太尼族/羟考酮/曲马多)
# 全库+肝病原疾患层: 总体ROR, 分药ROR, 转归构成, Weibull时间起点
library(data.table)
set.seed(20260917)
setDTthreads(8)
zip <- path.expand("~/nurse_faers/pmda_case_report_202608.zip")
outdir <- path.expand("~/nurse_faers/analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_j <- function(file, select_cols) {
  cmd <- sprintf("unzip -p %s %s | iconv -f CP932 -t UTF-8", shQuote(zip), file)
  fread(cmd = cmd, encoding = "UTF-8", select = select_cols, showProgress = FALSE)
}
demo <- read_j("demo202608.csv", c("識別番号","報告回数"))
demo[, rk := suppressWarnings(as.integer(報告回数))]
setorder(demo, 識別番号, rk)
keep <- demo[, .SD[.N], by = 識別番号][, .(識別番号, 報告回数)]
n_all <- nrow(keep)
kd <- function(dt) dt[keep, on = c("識別番号","報告回数"), nomatch = 0]

cls_map <- list(
  morphine = "モルヒネ",
  hydromorphone = "ヒドロモルフォン",
  fentanyl = "(フェンタニル|デュロテップ)",
  oxycodone = "(オキシコドン|オキシコンティン)",
  tramadol = "(トラマドール|トラマセット)")
any_opi <- paste0("(", paste(cls_map, collapse = "|"), ")")

drug <- kd(read_j("drug202608.csv", c("識別番号","報告回数","医薬品の関与",
        "医薬品（一般名）","医薬品（販売名）","投与開始日")))
sus <- drug[医薬品の関与 == "被疑薬"]
opi_rows <- sus[grepl(any_opi, `医薬品（一般名）`) | grepl(any_opi, `医薬品（販売名）`)]
opi_rows[, cls := NA_character_]
for (cl in names(cls_map))
  opi_rows[is.na(cls) & (grepl(cls_map[[cl]], `医薬品（一般名）`) |
            grepl(cls_map[[cl]], `医薬品（販売名）`)), cls := cl]
opi_ids <- unique(opi_rows$識別番号)
message("JADER 病例 ", n_all, "; 对齐药物集阿片被疑药 ", length(opi_ids))

pt_cns <- c("肝性脳症","傾眠","意識障害","錯乱状態","昏睡")
reac <- kd(read_j("reac202608.csv", c("識別番号","報告回数","有害事象","転帰","有害事象の発現日")))
ev_he  <- unique(reac[有害事象 == "肝性脳症"]$識別番号)
ev_cns <- unique(reac[有害事象 %chin% pt_cns]$識別番号)
message("CNS5 ", length(ev_cns), "; 肝性脳症 ", length(ev_he))

liv_rx <- "(肝細胞がん|肝細胞癌|肝癌|肝がん|肝硬変|肝不全)"
hist <- kd(read_j("hist202608.csv", c("識別番号","報告回数","原疾患等")))
liver_ids <- unique(hist[grepl(liv_rx, 原疾患等)]$識別番号)
message("肝病原疾患 ", length(liver_ids))

mk <- function(ev, expv, n) {
  a <- length(intersect(ev, expv)); cc <- length(setdiff(ev, expv))
  b <- length(setdiff(expv, ev));  d <- n - a - b - cc
  if (a < 3) return(data.table(a=a, b=b, c=cc, d=d, ROR=NA, lo=NA, hi=NA))
  a<-as.numeric(a); b<-as.numeric(b); cc<-as.numeric(cc); d<-as.numeric(d)
  ror <- (a*d)/(b*cc); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(a=a, b=b, c=cc, d=d, ROR=ror,
             lo=exp(log(ror)-1.96*se), hi=exp(log(ror)+1.96*se))
}
ror_tab <- rbindlist(list(
  data.table(layer="全库", drug="全阿片类(对齐集)", outcome="CNS5", mk(ev_cns, opi_ids, n_all)),
  data.table(layer="全库", drug="全阿片类(对齐集)", outcome="肝性脳症", mk(ev_he, opi_ids, n_all)),
  data.table(layer="肝病层", drug="全阿片类(对齐集)", outcome="CNS5",
             mk(ev_cns[ev_cns %in% liver_ids], opi_ids[opi_ids %in% liver_ids], length(liver_ids)))))
for (cl in names(cls_map)) {
  ex <- unique(opi_rows[cls == cl]$識別番号)
  ror_tab <- rbind(ror_tab,
    data.table(layer="全库", drug=cl, outcome="CNS5", mk(ev_cns, ex, n_all)))
}

# 转归构成(阿片×CNS5)
tar <- intersect(opi_ids, ev_cns)
oc_tab <- reac[識別番号 %chin% tar, .N, by = 転帰][order(-N)]
oc_tab[, pct := round(100*N/sum(N), 1)]

# Weibull 时间起点: 阿片投与開始日 -> CNS5発現日
weibull_fit <- function(t) {
  t <- t[t > 0]; n <- length(t)
  if (n < 20) return(list(n=n, shape=NA, scale=NA, median=NA, q1=NA, q3=NA, p_le7d=NA, p_le30d=NA))
  ll <- function(p) { k <- exp(p[1]); lam <- exp(p[2])
    -sum(log(k) - log(lam) + (k-1)*(log(t)-log(lam)) - (t/lam)^k) }
  o <- optim(log(c(1, max(median(t), 1))), ll)
  k <- exp(o$par[1]); lam <- exp(o$par[2])
  list(n=n, shape=k, scale=lam, median=qweibull(0.5,k,lam),
       q1=qweibull(0.25,k,lam), q3=qweibull(0.75,k,lam),
       p_le7d=mean(t<=7), p_le30d=mean(t<=30))
}
to_date <- function(x) as.Date(substr(as.character(x), 1, 8), format = "%Y%m%d", optional = TRUE)
start_opi <- opi_rows[!is.na(投与開始日) & 投与開始日 != "",
                      .(start = min(suppressWarnings(as.numeric(投与開始日)))), by = 識別番号]
ev_cns_dt <- reac[有害事象 %chin% pt_cns & !is.na(有害事象の発現日) & 有害事象の発現日 != "",
                  .(edate = min(suppressWarnings(as.numeric(有害事象の発現日)))), by = 識別番号]
tto <- merge(start_opi, ev_cns_dt, by = "識別番号")
tto[, d := as.numeric(to_date(edate) - to_date(start))]
tto <- tto[!is.na(d) & d >= 0 & d <= 730]
tto_out <- as.data.table(c(list(group = "JADER全库·阿片×CNS5"), as.list(weibull_fit(tto$d))))

fwrite(ror_tab, file.path(outdir, "jader_ror.csv"))
fwrite(oc_tab,  file.path(outdir, "jader_outcome_cns5.csv"))
fwrite(tto_out, file.path(outdir, "jader_tto.csv"))
message("JADER 分析完成")
print(ror_tab); print(oc_tab); print(tto_out)
