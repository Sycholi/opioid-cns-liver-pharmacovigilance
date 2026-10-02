# sensitivity_analysis.R —— 标配内容补齐
# A. 多指标一致性: ROR / PRR+chi2(Evans) / IC025(Bate1998 beta-binomial, PhViD解析式)
# B. FAERS 基线特征表(年龄/性别/报告人/国家/年份) + 报告人/性别/年龄/时期分层
# C. JADER 基线特征表 + 分层(職種/性别/年龄)
# TEST=1 时: A 全量, B 仅读前3个季度做解析诊断, C 仅打印職種原始分布, 不写 B/C 结果
library(data.table)
set.seed(20260917); setDTthreads(8)
TEST <- identical(Sys.getenv("TEST", "0"), "1")
zipdir <- path.expand("~/nurse_faers/ascii_zip")
outdir <- path.expand("~/nurse_faers/analysis")
jzip   <- path.expand("~/nurse_faers/pmda_case_report_202608.zip")

## ---------- 指标(向量化) ----------
meas <- function(a, b, c, d) {
  a<-as.numeric(a); b<-as.numeric(b); c<-as.numeric(c); d<-as.numeric(d)
  N <- a+b+c+d
  ror <- a*d/(b*c); se <- sqrt(1/a+1/b+1/c+1/d)
  prr <- (a/(a+b))/(c/(c+d)); sep <- sqrt(1/a-1/(a+b)+1/c-1/(c+d))
  chi2 <- N*(a*d-b*c)^2/((a+b)*(c+d)*(a+c)*(b+d))
  # IC: Bate 1998 beta-binomial, 同 PhViD BCPNN(MC=FALSE) 解析式
  p1<-1+a+b; p2<-1+N-a-b; q1<-1+a+c; q2<-1+N-a-c; r1<-1+a
  r2b <- N-a-1+(2+N)^2/(q1*p1)
  IC  <- (digamma(r1)-digamma(r1+r2b)
        -(digamma(p1)-digamma(p1+p2)+digamma(q1)-digamma(q1+q2)))/log(2)
  VIC <- (trigamma(r1)-trigamma(r1+r2b)
        +(trigamma(p1)-trigamma(p1+p2)+trigamma(q1)-trigamma(q1+q2)))/log(2)^2
  ic025 <- IC+qnorm(0.025)*sqrt(pmax(VIC,0))
  data.table(N=N, ROR=ror, ROR_lo=exp(log(ror)-1.96*se), ROR_hi=exp(log(ror)+1.96*se),
             PRR=prr, PRR_lo=exp(log(prr)-1.96*sep), PRR_hi=exp(log(prr)+1.96*sep),
             chi2=chi2, IC=IC, IC025=ic025,
             sig_ROR=exp(log(ror)-1.96*se)>1,
             sig_PRR=prr>=2 & chi2>=4 & a>=3,
             sig_IC=ic025>0)
}

## ---------- A. 多指标一致性 ----------
for (db in c("faers","jader")) {
  tab <- fread(file.path(outdir, sprintf("%s_ror.csv", db)))
  out <- cbind(tab, meas(tab$a, tab$b, tab$c, tab$d))
  fwrite(out, file.path(outdir, sprintf("%s_measures.csv", db)))
  message("\n== ", db, " 多指标一致性 ==")
  print(out[, .(layer, drug, outcome, a=a, ROR=sprintf("%.2f(%.2f-%.2f)",ROR,ROR_lo,ROR_hi),
                PRR=sprintf("%.2f(%.2f-%.2f)",PRR,PRR_lo,PRR_hi),
                chi2=sprintf("%.0f",chi2), IC025=sprintf("%.2f",IC025),
                信号一致=sig_ROR==sig_IC & sig_IC==sig_PRR)])
  message(db, " 三指标判定不一致的行数: ",
          out[sig_ROR!=sig_IC | sig_IC!=sig_PRR, .N])
}

## ---------- B. FAERS ----------
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
read_q <- function(zip, tag, want) {
  lst <- suppressWarnings(system(paste("unzip -Z1", shQuote(zip)), intern = TRUE))
  fn <- lst[grepl(paste0("(^|/)", tag, "[0-9].*\\.txt$"), lst, ignore.case = TRUE)][1]
  if (is.na(fn)) return(NULL)
  cmd <- paste("unzip -p", shQuote(zip), shQuote(fn))
  hdr <- fread(cmd = cmd, sep = "$", header = TRUE, quote = "", nrows = 0,
               fill = TRUE, showProgress = FALSE)
  nm <- norm(names(hdr))
  idx <- vapply(want, function(p) { i <- grep(paste0("^", p, "$"), nm)[1]
                                    if (is.na(i)) NA_integer_ else i }, integer(1))
  keep <- !is.na(idx); if (!keep[[1]]) return(NULL)
  dt <- tryCatch(
    fread(cmd = cmd, sep = "$", header = TRUE, quote = "", na.strings = c("", "NA"),
          fill = Inf, showProgress = FALSE, select = as.integer(idx[keep])),
    error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  setnames(dt, names(idx[keep])); dt[]
}

zips <- sort(list.files(zipdir, pattern = "\\.zip$", full.names = TRUE))
if (TEST) zips <- zips[1:3]
dwant <- c(primaryid="primaryid|isr", caseid="caseid|case", caseversion="caseversion",
           fda_dt="fdadt", age="age", age_cod="agecod", sex="sex|gndrcod",
           occp="occpcod", country="reportercountry")
demo <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "DEMO", dwant)
  if (is.null(d)) return(NULL)
  if (!"caseversion" %in% names(d)) d[, caseversion := 0L]
  has <- function(cl) if (cl %chin% names(d)) d[[cl]] else NA_character_
  data.table(
    primaryid = as.character(d$primaryid),
    caseid = as.character(d$caseid),
    caseversion = suppressWarnings(as.integer(d$caseversion)),
    fda_dt = suppressWarnings(as.integer(d$fda_dt)),
    age = suppressWarnings(as.numeric(has("age"))),
    age_cod = toupper(has("age_cod")),
    sex = toupper(has("sex")),
    occp = toupper(has("occp")),
    country = toupper(has("country")))
}), use.names = TRUE, fill = TRUE)
setorder(demo, caseid, caseversion, fda_dt, primaryid)
demo <- demo[, .SD[.N], by = caseid]      # 与主分析完全一致的去重
message("DEMO 去重后病例(测试=", TEST, "): ", nrow(demo))

if (TEST) {
  message("occp 分布: ");   print(demo[, .N, by = occp][order(-N)][1:8])
  message("sex 分布: ");    print(demo[, .N, by = sex])
  message("age_cod 分布: ");print(demo[, .N, by = age_cod])
  message("country 前8: "); print(demo[, .N, by = country][order(-N)][1:8])
}

if (!TEST) {
ck <- readRDS(file.path(outdir, "ckpt_ids.rds"))
keep_ids <- ck$keep_ids; exp_ids <- ck$exp_ids; ev_cns <- ck$ev_cns
liver_ids <- ck$liver_ids
demo <- demo[primaryid %chin% keep_ids]
demo[, age_y := {
  f <- c(YR=1, YEAR=1, MON=1/12, MONTH=1/12, WK=1/52.17857, WEEK=1/52.17857,
         DY=1/365.25, DAY=1/365.25, HR=1/8766, HOUR=1/8766)
  ay <- age * f[age_cod]
  fifelse(!is.na(ay) & ay>=0 & ay<=120, ay, NA_real_)
}]
demo[, age_grp := fcase(age_y<18, "<18", age_y<45, "18-44", age_y<65, "45-64",
                        age_y<75, "65-74", !is.na(age_y), ">=75", default="不明")]
demo[, rep_grp := fcase(occp %chin% c("MD","PH","OT"), "卫生专业人员",
                        occp=="CN", "消费者", occp=="LY", "律师",
                        !is.na(occp) & occp!="", "其他", default="缺失")]
demo[, sex_s := fcase(sex=="M", "男性", sex=="F", "女性", default="其他/缺失")]
demo[, yr := fcase(fda_dt>=20040101 & fda_dt<=20261231,
                   as.integer(fda_dt %/% 10000), default=NA_integer_)]
demo[, era := fcase(yr<=2014, "2004-2014", !is.na(yr), "2015-2026", default="无年份")]
demo[, reg := fcase(country=="US", "美国", country=="JP", "日本", country=="CN", "中国",
                    !is.na(country) & country!="", "其他", default="缺失")]

## 基线特征表: 阿片×CNS5 报告 vs 其余报告
demo[, tgt := primaryid %chin% intersect(exp_ids, ev_cns)]
t1 <- rbindlist(lapply(
  list(sex="sex_s", age_grp="age_grp", rep_grp="rep_grp", era="era", reg="reg"),
  function(v) {
    tv <- demo[, .N, by = c(v, "tgt")]
    dcast(tv, ... ~ tgt, value.var = "N")[
      , .(variable = v, level = get(v), n_target = `TRUE`, n_other = `FALSE`)][
      , `:=`(pct_target = round(100*n_target/sum(n_target),1),
             pct_other  = round(100*n_other/sum(n_other),1))]
  }))
fwrite(t1, file.path(outdir, "faers_table1.csv"))
message("\n== FAERS 基线特征 =="); print(t1)

## 分层 2x2
strat_one <- function(layer, universe, expset, evset, svar) {
  sub <- demo[primaryid %chin% universe]
  sub[, e := primaryid %chin% expset]; sub[, v := primaryid %chin% evset]
  rbindlist(lapply(sort(unique(sub[[svar]])), function(lv) {
    s <- sub[get(svar)==lv]
    a <- sum(s$e & s$v); b <- sum(s$e & !s$v)
    cc <- sum(!s$e & s$v); d <- sum(!s$e & !s$v)
    if (a < 3 || b < 3 || cc < 3 || d < 3)
      return(data.table(layer=layer, stratum=svar, level=lv, a=a,b=b,c=cc,d=d, N=a+b+cc+d,
                        ROR=NA_real_,ROR_lo=NA_real_,ROR_hi=NA_real_,PRR=NA_real_,
                        PRR_lo=NA_real_,PRR_hi=NA_real_,chi2=NA_real_,IC=NA_real_,
                        IC025=NA_real_,sig_ROR=NA,sig_PRR=FALSE,sig_IC=NA))
    cbind(data.table(layer=layer, stratum=svar, level=lv, a=a,b=b,c=cc,d=d), meas(a,b,cc,d))
  }))
}
st <- rbindlist(list(
  strat_one("全库", keep_ids, exp_ids, ev_cns, "rep_grp"),
  strat_one("全库", keep_ids, exp_ids, ev_cns, "sex_s"),
  strat_one("全库", keep_ids, exp_ids, ev_cns, "age_grp"),
  strat_one("全库", keep_ids, exp_ids, ev_cns, "era"),
  strat_one("肝病层", liver_ids, exp_ids, ev_cns, "rep_grp"),
  strat_one("肝病层", liver_ids, exp_ids, ev_cns, "sex_s")))
fwrite(st, file.path(outdir, "faers_stratified.csv"))
message("\n== FAERS 分层 =="); print(st)
}

## ---------- C. JADER ----------
read_j <- function(file, select_cols) {
  cmd <- sprintf("unzip -p %s %s | iconv -f CP932 -t UTF-8", shQuote(jzip), file)
  fread(cmd = cmd, encoding = "UTF-8", select = select_cols, showProgress = FALSE)
}
hdrcmd <- sprintf("unzip -p %s demo202608.csv | head -1 | iconv -f CP932 -t UTF-8", shQuote(jzip))
jcols <- strsplit(system(hdrcmd, intern = TRUE), ",", fixed = TRUE)[[1]]
pick <- function(pat) grep(pat, jcols, value = TRUE)[1]
c_sex <- pick("^性別$"); c_age <- pick("^年齢$"); c_occ <- pick("報告者の資格")
message("JADER demo 列定位: 性別=", c_sex, " 年齢=", c_age, " 資格=", c_occ)
demo_j <- read_j("demo202608.csv", c("識別番号","報告回数", c_sex, c_age, c_occ))
setnames(demo_j, c(c_sex,c_age,c_occ), c("sex","age","occ"))
demo_j[, rk := suppressWarnings(as.integer(報告回数))]
setorder(demo_j, 識別番号, rk); demo_j <- demo_j[, .SD[.N], by = 識別番号]
n_all_j <- nrow(demo_j)

# 病例标记: 与 jader_analysis.R 相同逻辑
cls_map <- list(morphine="モルヒネ", hydromorphone="ヒドロモルフォン",
                fentanyl="(フェンタニル|デュロテップ)",
                oxycodone="(オキシコドン|オキシコンティン)",
                tramadol="(トラマドール|トラマセット)")
any_opi <- paste0("(", paste(cls_map, collapse="|"), ")")
kd <- function(dt) dt[demo_j[, .(識別番号, 報告回数)], on=c("識別番号","報告回数"), nomatch=0]
drug <- kd(read_j("drug202608.csv", c("識別番号","報告回数","医薬品の関与",
        "医薬品（一般名）","医薬品（販売名）")))
sus <- drug[医薬品の関与 == "被疑薬"]
opi_ids <- unique(sus[grepl(any_opi, `医薬品（一般名）`) | grepl(any_opi, `医薬品（販売名）`)]$識別番号)
reac <- kd(read_j("reac202608.csv", c("識別番号","報告回数","有害事象")))
pt_cns <- c("肝性脳症","傾眠","意識障害","錯乱状態","昏睡")
ev_cns_j <- unique(reac[有害事象 %chin% pt_cns]$識別番号)
message("JADER 病例 ", n_all_j, "; 阿片被疑药 ", length(opi_ids), "; CNS5 ", length(ev_cns_j))
message("報告者の資格 原始分布(全部报告): "); print(demo_j[, .N, by = occ][order(-N)])
if (TEST) {
  message("性別: "); print(demo_j[, .N, by = sex])
  message("年齢 summary: "); print(summary(suppressWarnings(as.numeric(demo_j$age))))
} else {
demo_j[, occ_grp := fcase(grepl("医師|薬剤師|看護師|その他の医療|医療関係者|医療従事者|保健", occ), "医疗专业人员",
                          grepl("消費者|消费者|本人|家族", occ), "消费者",
                          grepl("^不明$|^$|不明", occ), "缺失",
                          !is.na(occ) & occ!="", "其他", default="缺失")]
demo_j[, sex_c := fcase(sex=="男性","男性", sex=="女性","女性", default="不明")]
demo_j[, age_grp := fcase(grepl("新生児|乳幼児|小児|10歳未満|トライメスター", age), "<10岁",
                          age=="10歳代", "10-19", age=="20歳代", "20-29",
                          age=="30歳代", "30-39", age=="40歳代", "40-49",
                          age=="50歳代", "50-59", age=="60歳代", "60-69",
                          grepl("^70歳代$|^80歳代$|^90歳代$|^100歳代$", age), "70+",
                          age=="高齢者", "高龄(不明档)", default="不明(含成人)")]
demo_j[, age65g := fcase(grepl("^[6-9]0歳代$|^100歳代$|^高齢者$", age), "≥60代",
                        grepl("歳代|歳未満|小児|新生児|乳幼児|青少年|トライメスター", age), "<60代",
                        default="不明(含成人)")]
demo_j[, tgt := 識別番号 %chin% intersect(opi_ids, ev_cns_j)]
t1j <- rbindlist(lapply(list(sex_c="sex_c", occ_grp="occ_grp", age_grp="age_grp"),
  function(v) {
    tv <- demo_j[, .N, by = c(v, "tgt")]
    dcast(tv, ... ~ tgt, value.var = "N")[
      , .(variable=v, level=get(v), n_target=`TRUE`, n_other=`FALSE`)][
      , `:=`(pct_target=round(100*n_target/sum(n_target),1),
             pct_other=round(100*n_other/sum(n_other),1))]
  }))
fwrite(t1j, file.path(outdir, "jader_table1.csv"))
fwrite(demo_j[, .N, by = occ][order(-N)], file.path(outdir, "jader_reporter_raw.csv"))
strat_j <- function(svar) {
  sub <- copy(demo_j)
  sub[, e := 識別番号 %chin% opi_ids]; sub[, v := 識別番号 %chin% ev_cns_j]
  rbindlist(lapply(sort(unique(sub[[svar]])), function(lv) {
    s <- sub[get(svar)==lv]
    a <- sum(s$e & s$v); b <- sum(s$e & !s$v); cc <- sum(!s$e & s$v); d <- sum(!s$e & !s$v)
    if (a<3 || b<3 || cc<3 || d<3)
      return(data.table(layer="全库", stratum=svar, level=lv, a=a,b=b,c=cc,d=d, N=a+b+cc+d,
                        ROR=NA_real_,ROR_lo=NA_real_,ROR_hi=NA_real_,PRR=NA_real_,
                        PRR_lo=NA_real_,PRR_hi=NA_real_,chi2=NA_real_,IC=NA_real_,
                        IC025=NA_real_,sig_ROR=NA,sig_PRR=FALSE,sig_IC=NA))
    cbind(data.table(layer="全库", stratum=svar, level=lv, a=a,b=b,c=cc,d=d), meas(a,b,cc,d))
  }))
}
stj <- rbindlist(list(strat_j("occ_grp"), strat_j("sex_c"), strat_j("age65g")))
fwrite(stj, file.path(outdir, "jader_stratified.csv"))
message("\n== JADER 基线特征 =="); print(t1j)
message("\n== JADER 分层 =="); print(stj)
}
message("敏感性分析完成 TEST=", TEST)
