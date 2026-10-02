# sens2_analysis.R —— v2 rerun (2026-10-02): strict liver definition, report-level outcomes, same-day TTO
# 修复: 1)肝层正则剔除CANCER PAIN+三分解 2)结局构成报告级分母 3)Weibull t>=0
# 增补: 4)扩展PT敏感性(CNS9) 5)PS+SS暴露敏感性 6)JADER结局/TTO修复
library(data.table)
source(path.expand("~/nurse_faers/read_drug_fixed.R"))
set.seed(20261002); setDTthreads(8)
TEST <- identical(Sys.getenv("TEST", "0"), "1")
zipdir <- path.expand("~/nurse_faers/ascii_zip")
outdir <- path.expand("~/nurse_faers/analysis")
jzip   <- path.expand("~/nurse_faers/pmda_case_report_202608.zip")
zips <- sort(list.files(zipdir, pattern = "\\.zip$", full.names = TRUE))
if (TEST) zips <- zips[1:3]

norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
read_q <- function(zip, tag, want) {
  lst <- suppressWarnings(system(paste("unzip -Z1", shQuote(zip)), intern = TRUE))
  fn <- lst[grepl(paste0("(^|/)", tag, "[0-9].*\\.txt$"), lst, ignore.case = TRUE)][1]
  if (is.na(fn)) return(NULL)
  cmd <- paste("unzip -p", shQuote(zip), shQuote(fn))
  hdr <- fread(cmd = cmd, sep = "$", header = TRUE, quote = "", nrows = 0, fill = TRUE, showProgress = FALSE)
  nm <- norm(names(hdr))
  idx <- vapply(want, function(p) { i <- grep(paste0("^", p, "$"), nm)[1]
                                    if (is.na(i)) NA_integer_ else i }, integer(1))
  keep <- !is.na(idx); if (!keep[[1]]) return(NULL)
  dt <- tryCatch(fread(cmd = cmd, sep = "$", header = TRUE, quote = "", na.strings = c("", "NA"),
                       fill = Inf, showProgress = FALSE, select = as.integer(idx[keep])),
                 error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  setnames(dt, names(idx[keep])); dt[]
}
ck <- readRDS(file.path(outdir, "ckpt_ids.rds"))
keep_ids <- ck$keep_ids; exp_ids <- ck$exp_ids; ev_cns <- ck$ev_cns
ev_he <- ck$ev_he; opi_rows <- ck$opi_cls; n_all <- ck$n_all

## ---------- 1. INDI re-read + strict liver layer ----------
liver_ind_strict <- "^(LIVER (CANCER|NEOPLASM|CARCINOMA|TUMOUR|TUMOR)|HEPATIC CIRRHOSIS|CIRRHOSIS|HEPATIC CANCER|HEPATOCELLULAR CARCINOMA|HEPATIC FAILURE|LIVER FAILURE)"
liver_pt_rx <- "^(CIRRHOSIS|HEPATIC CIRRHOSIS|LIVER CIRRHOSIS|HEPATOCELLULAR CARCINOMA|LIVER DISORDER|HEPATIC FAILURE|PORTAL HYPERTENSION|OESOPHAGEAL VARICES|ESOPHAGEAL VARICES|ASCITES)"
indi <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "INDI", c(primaryid = "primaryid|isr", indi_pt = "indipt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), indi_pt = toupper(indi_pt))]
}), use.names = TRUE, fill = TRUE)
indi <- indi[primaryid %chin% keep_ids]

## ---------- 2. REAC 重读 (CNS9 + 肝PT + JADER不需要) ----------
pt_cns5 <- c("HEPATIC ENCEPHALOPATHY","SOMNOLENCE","ALTERED STATE OF CONSCIOUSNESS","CONFUSIONAL STATE","COMA")
pt_ext  <- c("DEPRESSED LEVEL OF CONSCIOUSNESS","MENTAL STATUS CHANGES","LETHARGY","OBTUNDATION","STUPOR")
reac <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "REAC", c(primaryid = "primaryid|isr", pt = "pt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), pt_u = toupper(trimws(pt)))]
}), use.names = TRUE, fill = TRUE)
reac <- reac[primaryid %chin% keep_ids]
ev_cns5 <- unique(reac[pt_u %chin% pt_cns5]$primaryid)
ev_cns9 <- unique(reac[pt_u %chin% c(pt_cns5, pt_ext)]$primaryid)
liver_pt_ids <- unique(reac[grepl(liver_pt_rx, pt_u)]$primaryid)
liver_ind_ids <- unique(indi[grepl(liver_ind_strict, indi_pt)]$primaryid)
liverA <- union(as.character(liver_ind_ids), as.character(liver_pt_ids))  # primary definition: strict union
liverB <- as.character(liver_ind_ids)                                     # INDI-only
liverC <- as.character(liver_pt_ids)                                      # PT-only
message("liver layer sizes: primary ", length(liverA), "; INDI ", length(liverB), "; PT ", length(liverC))

mk <- function(ev, ex, n) {
  a <- length(intersect(ev, ex)); cc <- length(setdiff(ev, ex))
  b <- length(setdiff(ex, ev));  d <- n - a - b - cc
  a<-as.numeric(a); b<-as.numeric(b); cc<-as.numeric(cc); d<-as.numeric(d)
  ror <- (a*d)/(b*cc); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(a=a,b=b,c=cc,d=d,ROR=ror,lo=exp(log(ror)-1.96*se),hi=exp(log(ror)+1.96*se))
}
mkx <- function(ev, ex, n) {
  a <- length(intersect(ev,ex))
  if (a < 3) return(data.table(a=a,b=NA_real_,c=NA_real_,d=NA_real_,ROR=NA_real_,lo=NA_real_,hi=NA_real_))
  mk(ev,ex,n)
}
opl <- readRDS(file.path(outdir, "opi_link.rds"))
exp_of <- function(cl) unique(opl[cls == cl]$primaryid)
cls_vec <- c("morphine","hydromorphone","fentanyl","oxycodone","tramadol","hydrocodone")
liv_defs <- list(primary_union=liverA, INDI_only=liverB, PT_only=liverC)
tab_liv <- rbindlist(lapply(names(liv_defs), function(nm2) {
  L <- liv_defs[[nm2]]
  rbindlist(list(
    data.table(def=nm2, drug="全阿片类", mkx(ev_cns5[ev_cns5 %in% L], exp_ids[exp_ids %in% L], length(L))),
    rbindlist(lapply(cls_vec, function(cl) {
      ex <- exp_of(cl)
      data.table(def=nm2, drug=cl, mkx(ev_cns5[ev_cns5 %in% L], ex[ex %in% L], length(L)))
    }))))
}), fill = TRUE)
fwrite(tab_liv, file.path(outdir, "faers2_liver_ror.csv"))

## ---------- 3. CNS9 / PS+SS 敏感性 ----------
tab_ext <- rbindlist(list(
  data.table(set="CNS5", mkx(ev_cns5, exp_ids, n_all)),
  data.table(set="CNS9扩展", mkx(ev_cns9, exp_ids, n_all))))
if (!TEST) {
  drug_ps2 <- rbindlist(lapply(zips, function(z) {
    d <- read_drug(z)
    if (is.null(d)) return(NULL)
    d[!is.na(rolecod) & toupper(rolecod) %chin% c("PS","SS"),
      .(primaryid = as.character(primaryid),
        drugname = toupper(drugname),
        prod_ai = toupper(if ("prodai" %in% names(d)) prodai else NA_character_))]
  }), use.names = TRUE, fill = TRUE)
  drug_ps2 <- drug_ps2[primaryid %chin% keep_ids]
  anyopi <- paste0("(", paste(c("^(MORPHINE|MS ?CONTIN)","^(HYDROMORPHONE|DILAUDID)","^(FENTANYL|SUFENTANIL|REMIFENTANIL|DURAGESIC|SUBSYS|ACTIQ)","^(OXYCODONE|OXYCONTIN)","^(TRAMADOL|ULTRAM)","^HYDROCODONE"), collapse = "|"), ")")
  exp_ps_ss <- unique(drug_ps2[grepl(anyopi, drugname) | grepl(anyopi, prod_ai)]$primaryid)
  tab_ext <- rbind(tab_ext, data.table(set=paste0("PS+SS(", length(exp_ps_ss), "报告)"), mkx(ev_cns5, exp_ps_ss, n_all)))
  fwrite(tab_ext, file.path(outdir, "faers2_cns9_ps_ss.csv"))
}

## ---------- 4. OUTC 报告级分母 ----------
oc_map <- function(x) {
  y <- rep("未编码", length(x))
  y[grepl("^OT$|^OTHER", x)] <- "其他严重"
  y[grepl("^HO$|^HOSP", x)] <- "住院"
  y[grepl("^DE$|DEAT|DIED", x)] <- "死亡"
  y[grepl("^LT$|LIFE", x)] <- "危及生命"
  y[grepl("^DS$|DISA", x)] <- "致残"
  y[grepl("CONGEN|^CA$", x)] <- "先天异常"
  y[grepl("REQUIR|^RI$", x)] <- "需干预"
  y }
outc <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "OUTC", c(primaryid = "primaryid|isr", oc = "outcod|outccod|outcomecode"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), oc = toupper(oc))]
}), use.names = TRUE, fill = TRUE)
outc <- outc[primaryid %chin% keep_ids]
outc[, ocg := oc_map(oc)]
tar <- intersect(exp_ids, ev_cns5); tar_liv <- intersect(tar, liverA)
rep_level <- function(ids) {
  n_rep <- uniqueN(ids)
  outc[primaryid %chin% ids, .(N = uniqueN(primaryid)), by = ocg][
    , pct := round(100*N/n_rep, 1)][order(-pct)]
}
oc_rep <- rep_level(tar); oc_rep_liv <- rep_level(tar_liv)
fwrite(oc_rep, file.path(outdir, "faers2_outcome_reports.csv"))
fwrite(oc_rep_liv, file.path(outdir, "faers2_outcome_liver_reports.csv"))

## ---------- 5. THER 重读 + Weibull t>=0 ----------
ther <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "THER", c(primaryid = "primaryid|isr", dsgdrugseq = "dsgdrugseq|drugseq",
                           start = "startdt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid),
        drug_seq = suppressWarnings(as.integer(dsgdrugseq)),
        st = suppressWarnings(as.integer(start)))]
}), use.names = TRUE, fill = TRUE)
seqs <- opl[, .(primaryid, drug_seq)]
starts <- ther[seqs, on = c("primaryid", "drug_seq"), nomatch = 0][
  !is.na(st) & st >= 19000101 & st <= 21001231]
start_opi <- starts[, .(st = min(st)), by = primaryid]
demo_ev <- ck$ev_dt[primaryid %chin% keep_ids]
tto_dt <- merge(start_opi, demo_ev, by = "primaryid")
tto_dt[, d := as.numeric(as.Date(substr(as.character(event_dt),1,8), format="%Y%m%d") -
                         as.Date(substr(as.character(st),1,8), format="%Y%m%d"))]
tto_dt <- tto_dt[!is.na(d) & d >= 0 & d <= 730]
weib0 <- function(t) {
  n0 <- sum(t >= 0); t <- pmax(t[t >= 0], 0.5)  # 同日事件按0.5天连续性校正入模
  n <- length(t)
  if (n < 20) return(list(n=n, shape=NA, scale=NA, median=NA, p_le7d=NA, p_le30d=NA))
  ll <- function(p) { k <- exp(p[1]); lam <- exp(p[2])
    -sum(log(k) - log(lam) + (k-1)*(log(t)-log(lam)) - (t/lam)^k) }
  o <- tryCatch(optim(log(c(1, max(median(t), 1))), ll), error = function(e) NULL)
  if (is.null(o)) return(list(n=n, shape=NA, scale=NA, median=NA, p_le7d=NA, p_le30d=NA))
  k <- exp(o$par[1]); lam <- exp(o$par[2])
  list(n=n, shape=k, scale=lam, median=qweibull(0.5,k,lam),
       p_le7d=mean(t<=7), p_le30d=mean(t<=30))
}
tto2 <- rbindlist(list(
  data.table(group="FAERS全库·CNS5·t>=0", as.list(weib0(tto_dt[primaryid %chin% intersect(exp_ids, ev_cns5)]$d))),
  data.table(group="FAERS_liver_primary_CNS5_tge0", as.list(weib0(tto_dt[primaryid %chin% intersect(intersect(exp_ids, ev_cns5), liverA)]$d)))))
fwrite(tto2, file.path(outdir, "faers2_tto.csv"))
message("\n== strict liver layer ROR =="); print(tab_liv)
message("== CNS9/PS+SS =="); print(tab_ext)
message("== 报告级结局 =="); print(oc_rep); print(oc_rep_liv)
message("== TTO t>=0 =="); print(tto2)
message("sens2 完成 TEST=", TEST)

## ---------- 6. JADER 修复(t>=0 + 报告级结局) ----------
read_j <- function(file, sel) fread(cmd = sprintf("unzip -p %s %s | iconv -f CP932 -t UTF-8",
        shQuote(jzip), file), encoding = "UTF-8", select = sel, showProgress = FALSE)
demoJ <- read_j("demo202608.csv", c("識別番号","報告回数"))
demoJ[, rk := suppressWarnings(as.integer(報告回数))]
setorder(demoJ, 識別番号, rk); kj <- demoJ[, .SD[.N], by = 識別番号][, .(識別番号, 報告回数)]
kdJ <- function(dt) dt[kj, on = c("識別番号","報告回数"), nomatch = 0]
clsJ <- list(morphine="モルヒネ", hydromorphone="ヒドロモルフォン",
             fentanyl="(フェンタニル|デュロテップ)", oxycodone="(オキシコドン|オキシコンティン)",
             tramadol="(トラマドール|トラマセット)")
anyJ <- paste0("(", paste(clsJ, collapse="|"), ")")
drugJ <- kdJ(read_j("drug202608.csv", c("識別番号","報告回数","医薬品の関与","医薬品（一般名）","医薬品（販売名）","投与開始日")))
susJ <- drugJ[医薬品の関与 == "被疑薬"]
opiJ <- susJ[grepl(anyJ, `医薬品（一般名）`) | grepl(anyJ, `医薬品（販売名）`)]
startJ <- opiJ[!is.na(投与開始日) & 投与開始日 != "",
               .(st = min(suppressWarnings(as.numeric(投与開始日)))), by = 識別番号]
ptJ <- c("肝性脳症","傾眠","意識障害","錯乱状態","昏睡")
reacJ <- kdJ(read_j("reac202608.csv", c("識別番号","報告回数","有害事象","転帰","有害事象の発現日")))
evJ <- unique(reacJ[有害事象 %chin% ptJ]$識別番号)
tarJ <- intersect(unique(opiJ$識別番号), evJ)
edateJ <- reacJ[有害事象 %chin% ptJ & !is.na(有害事象の発現日) & 有害事象の発現日 != "",
                .(ed = min(suppressWarnings(as.numeric(有害事象の発現日)))), by = 識別番号]
ttoJ <- merge(startJ, edateJ, by = "識別番号")
ttoJ[, d := as.numeric(as.Date(substr(as.character(ed),1,8), format="%Y%m%d") -
                       as.Date(substr(as.character(st),1,8), format="%Y%m%d"))]
ttoJ <- ttoJ[!is.na(d) & d >= 0 & d <= 730]
jto2 <- data.table(group="JADER·CNS5·t>=0", as.list(weib0(ttoJ$d)))
fwrite(jto2, file.path(outdir, "jader2_tto.csv"))
ocJ <- reacJ[識別番号 %chin% tarJ, .(id = 識別番号, oc = 転帰)]
ocJ_rep <- ocJ[, .(N = uniqueN(id)), by = oc][, pct := round(100*N/uniqueN(tarJ), 1)][order(-pct)]
fwrite(ocJ_rep, file.path(outdir, "jader2_outcome_reports.csv"))
message("== JADER t>=0 =="); print(jto2); print(ocJ_rep)
