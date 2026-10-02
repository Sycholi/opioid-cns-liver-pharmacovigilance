# gate_faers.R —— FAERS/AERS 全量去重与阿片×肝性脑病门控统计 (R 4.4 / data.table)
# 兼容两代schema: 2004-2011 AERS(ISR/CASE,无caseversion) 与 2012+ FAERS(primaryid/caseid/caseversion)
library(data.table)
set.seed(20260914)
setDTthreads(8)

zipdir <- path.expand("~/nurse_faers/ascii_zip")
outdir <- path.expand("~/nurse_faers/gate")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
zips <- sort(list.files(zipdir, pattern = "\\.zip$", full.names = TRUE))
message(Sys.time(), " 共 ", length(zips), " 个季度文件")

norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# 从zip读表; want=命名向量(c(规范名="归一化后表头正则")), 缺失可选列由调用方补NA
read_q <- function(zip, tag, want) {
  lst <- suppressWarnings(system(paste("unzip -Z1", shQuote(zip)), intern = TRUE))
  fn <- lst[grepl(paste0("(^|/)", tag, "[0-9].*\\.txt$"), lst, ignore.case = TRUE)][1]
  if (is.na(fn)) return(NULL)
  # 2015+部分季度存在记录拼接行(>=4个连续$后接数字主键), 先拆行再读, 防止fread提前截断
  cmd <- paste("unzip -p", shQuote(zip), shQuote(fn),
               "| perl -pe 's/\\$\\$\\$\\$+(\\d{5,}\\$)/\\n\\$1/g'")
  hdr <- fread(cmd = cmd, sep = "$", header = TRUE, quote = "", nrows = 0,
               fill = TRUE, showProgress = FALSE)
  nm <- norm(names(hdr))
  idx <- vapply(want, function(p) { i <- grep(paste0("^", p, "$"), nm)[1]
                                    if (is.na(i)) NA_integer_ else i }, integer(1))
  keep <- !is.na(idx)
  if (!keep[[1]]) return(NULL)                 # 第一项必须是主键列
  dt <- tryCatch(fread(cmd = cmd, sep = "$", header = TRUE, quote = "",
                       na.strings = c("", "NA"), fill = TRUE, showProgress = FALSE,
                       select = as.integer(idx[keep])), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  setnames(dt, names(idx[keep]))
  dt[]
}

# ---- 1. DEMO: 每caseid保留最高版本 ----
demo_all <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "DEMO", c(primaryid = "primaryid|isr", caseid = "caseid|case",
                           caseversion = "caseversion", fda_dt = "fdadt"))
  if (is.null(d)) return(NULL)
  if (!"caseversion" %in% names(d)) d[, caseversion := 0L]
  d <- d[, .(primaryid, caseid,
             caseversion = suppressWarnings(as.integer(caseversion)),
             fda_dt = suppressWarnings(as.integer(fda_dt)))]
  d[!is.na(caseid) & !is.na(primaryid)]
}), use.names = TRUE, fill = TRUE)
setorder(demo_all, caseid, caseversion, fda_dt, primaryid)
keep_ids <- demo_all[, .SD[.N], by = caseid]$primaryid
message(Sys.time(), " DEMO行数 ", nrow(demo_all), "; 去重后病例数 ", length(keep_ids))

# ---- 2. DRUG: 主诉(PS)药物行 ----
opioid_rx <- "^(MORPHINE|HYDROMORPHONE|FENTANYL|OXYCODONE|TRAMADOL|HYDROCODONE|SUFENTANIL|REMIFENTANIL|DILAUDID|DURAGESIC|MS ?CONTIN|OXYCONTIN|ULTRAM|SUBSYS|ACTIQ)"
drug_ps <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "DRUG", c(primaryid = "primaryid|isr", drugname = "drugname",
                           prod_ai = "prodai", role_cod = "rolecod"))
  if (is.null(d)) return(NULL)
  d <- d[!is.na(role_cod) & toupper(role_cod) == "PS"]
  if (!"prod_ai" %in% names(d)) d[, prod_ai := NA_character_]
  d[, .(primaryid, drugname = toupper(drugname), prod_ai = toupper(prod_ai))]
}), use.names = TRUE, fill = TRUE)
drug_ps <- drug_ps[primaryid %chin% as.character(keep_ids)]
exp_ids <- drug_ps[grepl(opioid_rx, drugname) | grepl(opioid_rx, prod_ai), unique(primaryid)]
message(Sys.time(), " PS药物行 ", nrow(drug_ps), "; 阿片主诉报告 ", length(exp_ids))

# ---- 3. REAC ----
pt_cns <- c("Hepatic encephalopathy","Somnolence","Altered state of consciousness",
            "Confusional state","Coma")
reac <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "REAC", c(primaryid = "primaryid|isr", pt = "pt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid, pt = trimws(pt))]
}), use.names = TRUE, fill = TRUE)
reac <- reac[primaryid %chin% as.character(keep_ids)]
reac[, pt_u := toupper(pt)]   # 旧版AERS的PT为全大写, 新版首字母大写, 统一大写后匹配
ev_he  <- reac[pt_u == "HEPATIC ENCEPHALOPATHY", unique(primaryid)]
ev_cns <- reac[pt_u %chin% toupper(pt_cns), unique(primaryid)]
message(Sys.time(), " HE报告 ", length(ev_he), "; CNS5报告 ", length(ev_cns))

# ---- 4. INDI + 肝病人群代理 ----
liver_ind_rx <- "^(HEPATOCELLULAR|LIVER (CANCER|NEOPLASM|CARCINOMA|TUMOUR|TUMOR)|HEPATIC CIRRHOSIS|CIRRHOSIS|CANCER PAIN|HEPATIC CANCER)"
liver_pt_rx  <- "^(CIRRHOSIS|HEPATIC CIRRHOSIS|LIVER CIRRHOSIS|HEPATOCELLULAR CARCINOMA|LIVER DISORDER|HEPATIC FAILURE|PORTAL HYPERTENSION|OESOPHAGEAL VARICES|ESOPHAGEAL VARICES|ASCITES)"
indi <- rbindlist(lapply(zips, function(z) {
  d <- read_q(z, "INDI", c(primaryid = "primaryid|isr", indi_pt = "indipt"))
  if (is.null(d)) return(NULL)
  d[, .(primaryid, indi_pt = toupper(indi_pt))]
}), use.names = TRUE, fill = TRUE)
indi <- indi[primaryid %chin% as.character(keep_ids)]
liver_ind_ids <- indi[grepl(liver_ind_rx, indi_pt), unique(primaryid)]
liver_pt_ids  <- reac[grepl(liver_pt_rx, pt_u), unique(primaryid)]
liver_ids     <- union(liver_ind_ids, liver_pt_ids)
message(Sys.time(), " 肝病适应症报告 ", length(liver_ind_ids),
        "; 共报肝病PT报告 ", length(liver_pt_ids), "; 肝病人群合计 ", length(liver_ids))

# ---- 5. 门控计数与粗ROR ----
n_all  <- length(keep_ids); exp_ids <- as.character(exp_ids)
ev_he  <- as.character(ev_he);  ev_cns <- as.character(ev_cns)
liver_ids <- as.character(liver_ids)
mk <- function(ev) {
  a <- length(intersect(exp_ids, ev)); c <- length(setdiff(ev, exp_ids))
  b <- length(setdiff(exp_ids, ev));   d <- n_all - a - b - c
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  ror <- (a*d)/(b*c); se <- sqrt(1/a + 1/b + 1/c + 1/d)
  data.table(a=a, b=b, c=c, d=d, ROR=ror,
             ROR_lo=exp(log(ror) - 1.96*se), ROR_hi=exp(log(ror) + 1.96*se),
             P=fisher.test(matrix(c(a,b,c,d), 2))$p.value)
}
gate_tab <- rbindlist(list(
  data.table(outcome = "Hepatic_encephalopathy", mk(ev_he)),
  data.table(outcome = "CNS5_set",               mk(ev_cns))))
sub_tab <- data.table(
  item = c("全库去重报告","肝病人群报告","肝病×阿片主诉","肝病×阿片×CNS5","肝病×阿片×HE","全库阿片主诉"),
  n = c(n_all, length(liver_ids), length(intersect(exp_ids, liver_ids)),
        length(intersect(intersect(exp_ids, liver_ids), ev_cns)),
        length(intersect(intersect(exp_ids, liver_ids), ev_he)), length(exp_ids)))
# 肝病交叠诊断: 0是不可能的, 查类型与来源
diag <- data.table(
  item = c("liver_ids类", "exp_ids类", "liver∩exp(字符)",
           "liver∩exp(数值)", "IND来源肝病", "PT来源肝病"),
  val = c(class(liver_ids), class(exp_ids),
          length(intersect(liver_ids, exp_ids)),
          length(intersect(as.numeric(liver_ids), as.numeric(exp_ids))),
          length(intersect(as.character(liver_ind_ids), exp_ids)),
          length(intersect(as.character(liver_pt_ids), exp_ids))))
print(diag)
saveRDS(list(keep_ids = keep_ids, exp_ids = exp_ids, ev_he = ev_he,
             ev_cns = ev_cns, liver_ids = liver_ids,
             liver_ind_ids = as.character(liver_ind_ids),
             liver_pt_ids = as.character(liver_pt_ids)),
        file.path(outdir, "gate_intermediates.rds"))
fwrite(gate_tab, file.path(outdir, "gate_ror.csv"))
fwrite(sub_tab,  file.path(outdir, "gate_liver_subset.csv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message(Sys.time(), " 门控完成"); print(gate_tab); print(sub_tab)
