# inspire_validation.R —— INSPIRE 1.4.2 验证层分析
# 主线: 手术人群 → 术后首剂阿片(ATC N02A) → 病房GCS意识抑制事件(有分母发生率+时间窗+肝病/MELD分层)
# 输出仅聚合统计(PhysioNet DUA compliance); 患者级数据不出集群
library(data.table)
set.seed(20260930); setDTthreads(8)
TEST <- identical(Sys.getenv("TEST", "0"), "1")
d <- path.expand("~/nurse_inspire")
outdir <- path.expand("~/nurse_faers/analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------- 读取 ----------
ops <- fread(file.path(d, "operations.csv.gz"))
diag <- fread(file.path(d, "diagnosis.csv.gz"))
med <- fread(file.path(d, "medications.csv.gz"),
             select = c("subject_id","chart_time","drug_name","route",
                        "atc_code","atc_code2","atc_code3"))
wv <- fread(file.path(d, "ward_vitals.csv.gz"),
            select = c("subject_id","chart_time","item_name","value"))
labs <- fread(file.path(d, "labs.csv.gz"),
              select = c("subject_id","chart_time","item_name","value"))
for (dt in list(ops, diag, med, wv, labs)) dt[, subject_id := as.character(subject_id)]
message("读入: ops ", nrow(ops), "; diag ", nrow(diag), "; med ", nrow(med),
        "; ward_vitals ", nrow(wv), "; labs ", nrow(labs))

num <- function(x) suppressWarnings(as.numeric(x))
for (cl in c("opstart_time","opend_time","admission_time","discharge_time",
             "icuin_time","inhosp_death_time","age"))
  ops[, (cl) := num(get(cl))]
med[, chart_time := num(chart_time)]
wv[, chart_time := num(chart_time)]
labs[, chart_time := num(chart_time)]

## 时间单位自检(中位手术时长 2-6h: 数值~2-6=小时, ~120-360=分钟)
dur <- num(ops$opend_time - ops$opstart_time); dur <- dur[!is.na(dur)]
unit_h <- if (median(dur) > 20) 1/60 else 1
message("时间单位推断: 中位手术时长=", round(median(dur),1),
        " 原始单位 -> 小时换算系数=", unit_h)

## ---------- 队列: 单次手术成人 ----------
n_ops_subj <- ops[, .N, by = subject_id]
single <- n_ops_subj[N == 1, subject_id]
coh <- ops[subject_id %chin% single & !is.na(age) & age >= 18]
message("队列: ops ", nrow(ops), " -> 单次手术 ", length(single),
        " -> 成人 ", nrow(coh))

## ---------- 肝病标识 (ICD-10-CM: K70-K77 肝病 + C22 恶性) ----------
diag[, icd := toupper(sub("\\..*$", "", icd10_cm))]
liv_rx <- "^(K7[0-7]|C22)$"
liver_subj <- unique(diag[grepl(liv_rx, icd), subject_id])
coh[, liver := subject_id %chin% liver_subj]
message("肝病标识: ", sum(coh$liver), "/", nrow(coh))

## ---------- 阿片暴露 (任一ATC以N02A开头) ----------
atc_all <- paste(med$atc_code, med$atc_code2, med$atc_code3, sep = "|")
opi <- med[grepl("(^|\\|)N02A", atc_all)]
message("ATC N02A 给药行: ", nrow(opi), " / med 总行 ", nrow(med))
## 术后首剂: 在该次手术结束后、出院前
setkey(coh, subject_id)
opi2 <- opi[coh[, .(subject_id, opend = opend_time, disc = discharge_time)],
            on = "subject_id", nomatch = 0]
opi2 <- opi2[!is.na(chart_time) & !is.na(opend) &
             chart_time >= opend & (is.na(disc) | chart_time <= disc)]
first_opi <- opi2[, .(t_opi = min(chart_time)), by = subject_id]
coh <- first_opi[coh, on = "subject_id"]
coh[, exposed := !is.na(t_opi)]
message("术后阿片暴露: ", sum(coh$exposed), "/", nrow(coh))

## ---------- 病房GCS ----------
gcs <- wv[item_name %chin% c("gcs_e","gcs_m","gcs_v")]
gcs_w <- dcast(gcs, subject_id + chart_time ~ item_name, value.var = "value",
               fun.aggregate = function(x) x[1])
gcs_w[, gcs_total := num(gcs_e) + num(gcs_m) + num(gcs_v)]
gcs_w <- gcs_w[!is.na(gcs_total)]
message("完整GCS测量点: ", nrow(gcs_w), " (", uniqueN(gcs_w$subject_id), " 人)")

## 基线(首剂前最后一次GCS)与事件(首剂后首次GCS<=12/<=8)
base_gcs <- gcs_w[coh[, .(subject_id, t_opi, disc = discharge_time)],
                  on = "subject_id", nomatch = 0][chart_time < t_opi,
                  .(g_base = gcs_total[which.max(chart_time)]), by = subject_id]
coh <- base_gcs[coh, on = "subject_id"]
at_risk <- coh[exposed == TRUE & !is.na(g_base) & g_base >= 13]
message("风险队列(暴露且基线GCS>=13): ", nrow(at_risk))
gcs_ev <- gcs_w[at_risk[, .(subject_id, t_opi, disc = discharge_time)], on = "subject_id", nomatch = 0]
gcs_ev <- gcs_ev[chart_time > t_opi & (is.na(disc) | chart_time <= disc)]
ev12 <- gcs_ev[gcs_total <= 12, .(t_ev12 = min(chart_time)), by = subject_id]
ev8  <- gcs_ev[gcs_total <= 8,  .(t_ev8  = min(chart_time)), by = subject_id]
at_risk <- ev8[ev12[at_risk, on = "subject_id"], on = "subject_id"]
message("事件: GCS<=12 ", nrow(ev12), "; GCS<=8 ", nrow(ev8))

## 呼吸抑制次级结局(SpO2<90 或 RR<=8, 首剂后)
resp <- wv[item_name %chin% c("spo2","rr") & subject_id %chin% at_risk$subject_id]
resp <- resp[!is.na(num(value))]
resp_w <- dcast(resp, subject_id + chart_time ~ item_name, value.var = "value",
                fun.aggregate = function(x) x[1])
resp_w[, `:=`(spo2 = num(spo2), rr = num(rr))]
resp_w <- resp_w[at_risk[, .(subject_id, t_opi, disc = discharge_time)], on = "subject_id", nomatch = 0]
resp_w <- resp_w[chart_time > t_opi & (is.na(disc) | chart_time <= disc)]
resp_ev <- resp_w[(spo2 < 90) | (rr <= 8), .(t_resp = min(chart_time)), by = subject_id]
at_risk <- resp_ev[at_risk, on = "subject_id"]
message("呼吸抑制事件: ", nrow(resp_ev))

## ---------- MELD-Na (首剂前最近一次四项) ----------
need <- c("total_bilirubin","creatinine","ptinr","sodium")
l4 <- labs[item_name %chin% need]
l4[, value := num(value)]
lw <- dcast(l4, subject_id + chart_time ~ item_name, value.var = "value",
            fun.aggregate = function(x) x[1])
lw <- lw[at_risk[, .(subject_id, t_opi)], on = "subject_id", nomatch = 0]
lw <- lw[chart_time <= t_opi]
meld_in <- lw[, lapply(.SD, function(x) {
                v <- num(x); v[which.max(chart_time)] }), by = subject_id,
       .SDcols = need]
setnames(meld_in, need, c("bili","cr","inr","na"))
meld_in[, `:=`(bili=pmax(bili,1), cr=pmin(pmax(cr,1),4), inr=pmax(inr,1), na=pmax(na,125))]
meld_in[, meldi := 10*(0.957*log(cr) + 0.378*log(bili) + 1.120*log(inr) + 0.643)]
meld_in[, meld_na := pmin(pmax(meldi + 1.32*(137-na) - 0.033*meldi*(137-na), 6), 40)]
at_risk <- meld_in[, .(subject_id, meld_na)][at_risk, on = "subject_id"]
message("MELD-Na 可算: ", sum(!is.na(at_risk$meld_na)), "/", nrow(at_risk))

## ---------- 分药 ----------
opi_cls <- function(x) {
  x <- tolower(x)
  fcase(grepl("fentanyl", x), "fentanyl",
        grepl("morphine", x), "morphine",
        grepl("oxycodone", x), "oxycodone",
        grepl("tramadol", x), "tramadol",
        grepl("pethidine|meperidine", x), "pethidine",
        grepl("codeine", x), "codeine",
        grepl("hydromorphone", x), "hydromorphone",
        grepl("sufentanil|remifentanil|alfentanil", x), "fentanyl族(其他)",
        default = "其他")
}
first_drug <- opi2[opi2[, .I[chart_time == min(chart_time)], by = subject_id]$V1]
first_drug[, cls := opi_cls(drug_name)]

if (TEST) {
  message("\n== TEST 诊断 ==")
  print(coh[, .N, by = .(liver, exposed)])
  print(at_risk[, .(n=.N, ev12=sum(!is.na(t_ev12)), ev8=sum(!is.na(t_ev8)),
                   resp=sum(!is.na(t_resp)), death=sum(!is.na(inhosp_death_time)),
                   meld_avail=sum(!is.na(meld_na))), by = liver])
  print(at_risk[, .(meld_p25=round(quantile(meld_na,.25,na.rm=TRUE),1),
                    meld_med=round(median(meld_na,na.rm=TRUE),1),
                    meld_p75=round(quantile(meld_na,.75,na.rm=TRUE),1))])
  print(head(first_drug[, .N, by = cls][order(-N)], 10))
  print(at_risk[!is.na(t_ev12), .(hrs=(t_ev12-t_opi)*unit_h)][, summary(hrs)])
  quit(save = "no")
}

## ---------- 汇总输出(全部聚合) ----------
rr_ci <- function(a, n1, c, n0) {
  r <- (a/pmax(n1,1)) / (c/pmax(n0,1))
  se <- sqrt(1/pmax(a, 1e-9) - 1/n1 + 1/pmax(c, 1e-9) - 1/n0)
  out <- sprintf("%.2f (%.2f-%.2f)", r, exp(log(r)-1.96*se), exp(log(r)+1.96*se))
  out[a<3 | c<3] <- sprintf("%.2f (事件过少)", r[a<3 | c<3])
  out
}
tab_inc <- at_risk[, {
  n <- .N
  e1 <- sum(!is.na(t_ev12)); e2 <- sum(!is.na(t_ev8)); e3 <- sum(!is.na(t_resp))
  e4 <- sum(!is.na(inhosp_death_time))
  .(n=n, gcs_le12=e1, pct12=sprintf("%.1f%%",100*e1/n),
    gcs_le8=e2, pct8=sprintf("%.1f%%",100*e2/n),
    resp=e3, pct_resp=sprintf("%.1f%%",100*e3/n),
    death=e4, pct_death=sprintf("%.1f%%",100*e4/n))
}, by = liver]
A <- tab_inc[liver==TRUE]; B <- tab_inc[liver==FALSE]
tab_rr <- data.table(
  outcome = c("GCS<=12","GCS<=8","呼吸抑制","院内死亡"),
  a = c(A$gcs_le12, A$gcs_le8, A$resp, A$death),
  n1 = A$n, c = c(B$gcs_le12, B$gcs_le8, B$resp, B$death), n0 = B$n)
tab_rr[, RR := rr_ci(a, n1, c, n0)]

tim <- at_risk[!is.na(t_ev12), .(liver, hrs = (t_ev12 - t_opi) * unit_h)]
tab_time <- tim[, .(n=.N,
  median_h=round(median(hrs),1),
  p_le24h=sprintf("%.0f%%", 100*mean(hrs<=24)),
  p_le48h=sprintf("%.0f%%", 100*mean(hrs<=48)),
  p_le7d=sprintf("%.0f%%", 100*mean(hrs<=168))), by = liver]
## 事件时间直方(天桶)用于后续作图
tab_hist <- tim[, .(day_bucket = floor(hrs/24)), by = liver][
  , .N, by = .(liver, day_bucket)]

at_risk[, meld_grp := fcase(meld_na < 10, "<10", meld_na < 15, "10-14",
                            meld_na < 20, "15-19", !is.na(meld_na), ">=20",
                            default="缺失")]
tab_meld <- at_risk[, .(n=.N, ev=sum(!is.na(t_ev12)),
  pct=sprintf("%.1f%%",100*mean(!is.na(t_ev12)))), by=.(liver, meld_grp)]

fd <- first_drug[, .(subject_id, cls)][
  at_risk[, .(subject_id, liver, ev = !is.na(t_ev12))], on = "subject_id"]
tab_drug <- fd[, .(n=.N, ev=sum(ev), pct=sprintf("%.1f%%",100*mean(ev))),
               by=.(cls, liver)]

fwrite(tab_inc, file.path(outdir, "insp_incidence.csv"))
fwrite(tab_rr, file.path(outdir, "insp_rr.csv"))
fwrite(tab_time, file.path(outdir, "insp_timing.csv"))
fwrite(tab_hist, file.path(outdir, "insp_timing_hist.csv"))
fwrite(tab_meld, file.path(outdir, "insp_meld.csv"))
fwrite(tab_drug, file.path(outdir, "insp_drug.csv"))
message("\n== 发生率 =="); print(tab_inc)
message("\n== 肝病 vs 非肝病 RR =="); print(tab_rr)
message("\n== 事件时间(首剂后) =="); print(tab_time)
message("\n== MELD 分层 =="); print(tab_meld)
message("\n== 分药 =="); print(tab_drug)
message("INSPIRE 验证分析完成")
