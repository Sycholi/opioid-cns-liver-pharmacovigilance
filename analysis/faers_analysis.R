# faers_analysis.R —— FAERS 正式分析第一轮
# 全库+肝病层: 总体/分药 ROR(CNS5与HE), Weibull 时间起点, 严重结局构成
# 兼容 2004-2011 AERS 与 2012+ FAERS 两代 schema; 拼接行用角色码锚定拆分
library(data.table)
source(path.expand("~/nurse_faers/read_drug_fixed.R"))
set.seed(20260917)
setDTthreads(8)

zipdir <- path.expand("~/nurse_faers/ascii_zip")
outdir <- path.expand("~/nurse_faers/analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
zips <- sort(list.files(zipdir, pattern = "\\.zip$", full.names = TRUE))
message(Sys.time(), " 共 ", length(zips), " 个季度文件")
ckpt_file <- file.path(outdir, "ckpt_ids.rds")
RESUME <- FALSE
if (file.exists(ckpt_file)) {
  ck <- readRDS(ckpt_file)
  if (identical(ck$n_zips, length(zips))) RESUME <- TRUE
}

norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# 角色码锚定的拼接行拆分: 新记录起点 = $主键(7-9位)$药品序号(1-3位)$(PS|SS|C|I)$
JUNCTION <- "| perl -pe 's/\\$(\\d{7,9}\\$\\d{1,3}\\$(PS|SS|C|I)\\$)/\\n$1/g'"

read_q <- function(zip, tag, want, junction = FALSE) {
  lst <- suppressWarnings(system(paste("unzip -Z1", shQuote(zip)), intern = TRUE))
  fn <- lst[grepl(paste0("(^|/)", tag, "[0-9].*\\.txt$"), lst, ignore.case = TRUE)][1]
  if (is.na(fn)) return(NULL)
  cmd <- paste("unzip -p", shQuote(zip), shQuote(fn), if (junction) JUNCTION else "")
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

cls_map <- list(
  morphine = "^(MORPHINE|MS ?CONTIN)",
  hydromorphone = "^(HYDROMORPHONE|DILAUDID)",
  fentanyl = "^(FENTANYL|SUFENTANIL|REMIFENTANIL|DURAGESIC|SUBSYS|ACTIQ)",
  oxycodone = "^(OXYCODONE|OXYCONTIN)",
  tramadol = "^(TRAMADOL|ULTRAM)",
  hydrocodone = "^HYDROCODONE")
any_opioid_rx <- paste0("(", paste(cls_map, collapse = "|"), ")")

if (RESUME) {
  ck <- readRDS(ckpt_file)
  keep_ids <- ck$keep_ids; ev_dt <- ck$ev_dt; n_all <- ck$n_all
  exp_ids <- ck$exp_ids; ev_he <- ck$ev_he; ev_cns <- ck$ev_cns
  liver_ids <- ck$liver_ids
  opi_rows <- tryCatch(readRDS(file.path(outdir, "opi_link.rds")),
                       error = function(e) ck$opi_cls)
  message(Sys.time(), " 从检查点恢复: 季度数一致; opi_rows=",
          nrow(opi_rows), "行, 含drug_seq=", "drug_seq" %chin% names(opi_rows))
}

if (!RESUME) {
# ---- 1. DEMO: 去重 + 事件日期 ----
demo_all <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "DEMO", c(primaryid = "primaryid|isr", caseid = "caseid|case",
                           caseversion = "caseversion", fda_dt = "fdadt",
                           event_dt = "eventdt"))
  if (is.null(d)) return(NULL)
  if (!"caseversion" %in% names(d)) d[, caseversion := 0L]
  d[, .(primaryid, caseid,
        caseversion = suppressWarnings(as.integer(caseversion)),
        fda_dt = suppressWarnings(as.integer(fda_dt)),
        event_dt = suppressWarnings(as.integer(event_dt)))]
}), use.names = TRUE, fill = TRUE)
setorder(demo_all, caseid, caseversion, fda_dt, primaryid)
keepdemo <- demo_all[, .SD[.N], by = caseid]
keep_ids <- as.character(keepdemo$primaryid)
ev_dt <- keepdemo[, .(primaryid = as.character(primaryid), event_dt)]
ev_dt <- ev_dt[event_dt >= 19000101 & event_dt <= 21001231]
n_all <- length(keep_ids)
message(Sys.time(), " 去重后病例 ", n_all)
}

if (!RESUME) {
# ---- 2. DRUG: PS行 + 阿片分类 ----

if (!RESUME) drug_ps <- rbindlist(lapply(zips, function(z) {
  d <- read_drug(z)          # 结构化拼接行修复的读取器, 列名已规范化
  if (is.null(d)) return(NULL)
  d[!is.na(rolecod) & toupper(rolecod) == "PS",
    .(primaryid = as.character(primaryid),
      drug_seq = suppressWarnings(as.integer(drugseq)),
      drugname = toupper(drugname),
      prod_ai = toupper(if ("prodai" %in% names(d)) prodai else NA_character_))]
}), use.names = TRUE, fill = TRUE)
drug_ps <- drug_ps[primaryid %chin% keep_ids]
opi_rows <- drug_ps[grepl(any_opioid_rx, drugname) | grepl(any_opioid_rx, prod_ai)]
opi_rows[, cls := NA_character_]
for (cl in names(cls_map))
  opi_rows[is.na(cls) & (grepl(cls_map[[cl]], drugname) | grepl(cls_map[[cl]], prod_ai)), cls := cl]
exp_ids <- unique(opi_rows$primaryid)
message(Sys.time(), " 阿片主诉报告 ", length(exp_ids),
        " (分类缺失 ", sum(is.na(opi_rows$cls)), " 行)")
}

if (!RESUME) {
# ---- 3. REAC: CNS5 与 HE ----
pt_cns <- c("HEPATIC ENCEPHALOPATHY","SOMNOLENCE",
            "ALTERED STATE OF CONSCIOUSNESS","CONFUSIONAL STATE","COMA")
if (!RESUME) reac <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "REAC", c(primaryid = "primaryid|isr", pt = "pt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), pt_u = toupper(trimws(pt)))]
}), use.names = TRUE, fill = TRUE)
reac <- reac[primaryid %chin% keep_ids]
ev_he  <- unique(reac[pt_u == "HEPATIC ENCEPHALOPATHY"]$primaryid)
ev_cns <- unique(reac[pt_u %chin% pt_cns]$primaryid)
message(Sys.time(), " HE ", length(ev_he), "; CNS5 ", length(ev_cns))
}

if (!RESUME) {
# ---- 4. 肝病人群(id 先转字符再合并, 规避 integer64 union 损坏) ----
liver_ind_rx <- "^(HEPATOCELLULAR|LIVER (CANCER|NEOPLASM|CARCINOMA|TUMOUR|TUMOR)|HEPATIC CIRRHOSIS|CIRRHOSIS|CANCER PAIN|HEPATIC CANCER)"
liver_pt_rx  <- "^(CIRRHOSIS|HEPATIC CIRRHOSIS|LIVER CIRRHOSIS|HEPATOCELLULAR CARCINOMA|LIVER DISORDER|HEPATIC FAILURE|PORTAL HYPERTENSION|OESOPHAGEAL VARICES|ESOPHAGEAL VARICES|ASCITES)"
if (!RESUME) indi <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "INDI", c(primaryid = "primaryid|isr", indi_pt = "indipt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), indi_pt = toupper(indi_pt))]
}), use.names = TRUE, fill = TRUE)
indi <- indi[primaryid %chin% keep_ids]
liver_ind_ids <- unique(indi[grepl(liver_ind_rx, indi_pt)]$primaryid)
liver_pt_ids  <- unique(reac[grepl(liver_pt_rx, pt_u)]$primaryid)
liver_ids <- union(as.character(liver_ind_ids), as.character(liver_pt_ids))
message(Sys.time(), " 肝病人群 ", length(liver_ids))
ckpt <- file.path(outdir, "ckpt_ids.rds")
saveRDS(list(keep_ids=keep_ids, ev_dt=ev_dt, n_all=n_all,
             opi_cls=opi_rows[, .(primaryid, drug_seq, cls)], exp_ids=exp_ids,
             ev_he=ev_he, ev_cns=ev_cns, liver_ids=liver_ids,
             n_zips = length(zips)), ckpt)
}

# ---- 5. ROR: 总体 + 肝病层 + 分药 ----
mk <- function(ev, expv, n) {
  a <- length(intersect(ev, expv)); cc <- length(setdiff(ev, expv))
  b <- length(setdiff(expv, ev));  d <- n - a - b - cc
  if (a < 3) return(data.table(a=a, b=b, c=cc, d=d, ROR=NA, lo=NA, hi=NA))
  a<-as.numeric(a); b<-as.numeric(b); cc<-as.numeric(cc); d<-as.numeric(d)
  ror <- (a*d)/(b*cc); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(a=a, b=b, c=cc, d=d, ROR=ror,
             lo=exp(log(ror)-1.96*se), hi=exp(log(ror)+1.96*se))
}
ror_all <- rbindlist(list(
  data.table(layer="全库", drug="全阿片类", outcome="CNS5", mk(ev_cns, exp_ids, n_all)),
  data.table(layer="全库", drug="全阿片类", outcome="HE",   mk(ev_he,  exp_ids, n_all)),
  data.table(layer="肝病层", drug="全阿片类", outcome="CNS5",
             mk(ev_cns[ev_cns %in% liver_ids], exp_ids[exp_ids %in% liver_ids], length(liver_ids))),
  data.table(layer="肝病层", drug="全阿片类", outcome="HE",
             mk(ev_he[ev_he %in% liver_ids],  exp_ids[exp_ids %in% liver_ids],  length(liver_ids)))))
for (cl in names(cls_map)) {
  ex <- unique(opi_rows[cls == cl]$primaryid)
  ror_all <- rbind(ror_all,
    data.table(layer="全库", drug=cl, outcome="CNS5", mk(ev_cns, ex, n_all)),
    data.table(layer="肝病层", drug=cl, outcome="CNS5",
               mk(ev_cns[ev_cns %in% liver_ids], ex[ex %in% liver_ids], length(liver_ids))))
}

# ---- 6. Weibull 时间起点(事件日-阿片起始日, 0-730天) ----
weibull_fit <- function(t) {
  t <- t[t > 0]; n <- length(t)
  if (n < 20) return(c(n=n, shape=NA, scale=NA, median=NA, q1=NA, q3=NA,
                       p7=NA, p30=NA))
  ll <- function(p) { k <- exp(p[1]); lam <- exp(p[2])
    -sum(log(k) - log(lam) + (k-1)*(log(t)-log(lam)) - (t/lam)^k) }
  o <- optim(log(c(1, max(median(t), 1))), ll)
  k <- exp(o$par[1]); lam <- exp(o$par[2])
  c(n=n, shape=k, scale=lam, median=qweibull(0.5,k,lam),
    q1=qweibull(0.25,k,lam), q3=qweibull(0.75,k,lam),
    p7=mean(t<=7), p30=mean(t<=30))
}
to_date <- function(x) as.Date(substr(sprintf("%08d", as.integer(x)), 1, 8),
                              format = "%Y%m%d", optional = TRUE)
ther <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "THER", c(primaryid = "primaryid|isr", drug_seq = "drugseq|dsgdrugseq",
                           start_dt = "startdt"))
  if (is.null(d) || !"drug_seq" %in% names(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid),
        drug_seq = suppressWarnings(as.integer(drug_seq)),
        start_dt = suppressWarnings(as.integer(start_dt)))]
}), use.names = TRUE, fill = TRUE)
ther <- ther[!is.na(start_dt) & start_dt >= 19000101 & start_dt <= 21001231]
opi_start <- merge(opi_rows[, .(primaryid, drug_seq)], ther,
                   by = c("primaryid", "drug_seq"), all.x = TRUE)
opi_start <- opi_start[!is.na(start_dt)][, .(start = min(start_dt)), by = primaryid]
tto_dt <- merge(ev_dt, opi_start, by = "primaryid")
tto_dt[, tto := as.numeric(to_date(event_dt) - to_date(start))]
tto_dt <- tto_dt[!is.na(tto) & tto >= 0 & tto <= 730]
tto_tab <- rbind(
  as.data.table(c(list(group="全库·阿片×CNS5"), as.list(weibull_fit(tto_dt[primaryid %chin% intersect(exp_ids, ev_cns)]$tto)))),
  as.data.table(c(list(group="肝病·阿片×CNS5"), as.list(weibull_fit(tto_dt[primaryid %chin% intersect(intersect(exp_ids, ev_cns), liver_ids)]$tto)))),
  as.data.table(c(list(group="肝病·阿片×HE"),   as.list(weibull_fit(tto_dt[primaryid %chin% intersect(intersect(exp_ids, ev_he), liver_ids)]$tto)))))
setnames(tto_tab, names(tto_tab), c("group","n","shape","scale","median","q1","q3","p_le7d","p_le30d"))

# ---- 7. OUTC 严重结局构成(阿片×CNS5, 全库+肝病层) ----
outc <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "OUTC", c(primaryid = "primaryid|isr",
                           outc_cod = "outccod|outccode|outcomecode"))
  if (is.null(d) || !"outc_cod" %in% names(d)) return(NULL)
  d[, .(primaryid = as.character(primaryid), oc = toupper(outc_cod))]
}), use.names = TRUE, fill = TRUE)
outc <- outc[primaryid %chin% keep_ids]
oc_map <- function(x) {
  y <- rep("OT", length(x))
  y[grepl("DEATH|^DE$|^D$", x)] <- "死亡"
  y[grepl("LIFE|^LT$", x)] <- "危及生命"
  y[grepl("HOSP|^HO$", x)] <- "住院"
  y[grepl("DISAB|^DS$", x)] <- "致残"
  y[grepl("CONGEN|^CA$", x)] <- "先天异常"
  y[grepl("REQUIR|^RI$", x)] <- "需干预"
  y }
outc[, oc := oc_map(oc)]
tar <- intersect(exp_ids, ev_cns)
out_tab <- outc[primaryid %chin% tar, .N, by = oc][order(-N)]
out_tab[, pct := round(100*N/sum(N), 1)]
out_tab_liv <- outc[primaryid %chin% intersect(tar, liver_ids), .N, by = oc][order(-N)]
out_tab_liv[, pct := round(100*N/sum(N), 1)]

# ---- 8. 输出 ----
fwrite(ror_all, file.path(outdir, "faers_ror.csv"))
fwrite(tto_tab, file.path(outdir, "faers_tto.csv"))
fwrite(out_tab,  file.path(outdir, "faers_outcome_cns5.csv"))
fwrite(out_tab_liv, file.path(outdir, "faers_outcome_cns5_liver.csv"))
saveRDS(list(exp_ids=exp_ids, ev_he=ev_he, ev_cns=ev_cns, liver_ids=liver_ids,
             cls_rows=opi_rows), file.path(outdir, "faers_ids.rds"))
message(Sys.time(), " FAERS 正式分析完成")
print(ror_all); print(tto_tab); print(out_tab)
