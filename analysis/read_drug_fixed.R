# read_drug_fixed.R —— 结构化拼接行修复的 FAERS/AERS DRUG 读取器
# 原理: 全列读入(fill=Inf, 不会提前截断), 宽于表头的行=记录无缝拼接,
#        按表头列数整块重切; 输出列名为规范化小写(isr统一为primaryid)
library(data.table)
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

read_drug <- function(zip) {
  lst <- suppressWarnings(system(paste("unzip -Z1", shQuote(zip)), intern = TRUE))
  fn <- lst[grepl("(^|/)DRUG[0-9].*\\.txt$", lst, ignore.case = TRUE)][1]
  if (is.na(fn)) return(NULL)
  cmd <- paste("unzip -p", shQuote(zip), shQuote(fn))
  hdr <- fread(cmd = cmd, sep = "$", header = TRUE, quote = "", nrows = 0,
               fill = TRUE, showProgress = FALSE)
  canon <- norm(names(hdr))
  canon[canon == "isr"] <- "primaryid"
  canon[canon == "case"] <- "caseid"
  N0 <- length(canon)

  d <- fread(cmd = cmd, sep = "$", header = TRUE, quote = "",
             na.strings = c("", "NA"), fill = Inf, showProgress = FALSE)
  if (ncol(d) <= N0) { setnames(d, canon); return(d) }

  extra <- (N0 + 1):ncol(d)
  has_extra <- !is.na(d[[extra[1]]])
  for (j in extra[-1]) has_extra <- has_extra | !is.na(d[[j]])
  if (!any(has_extra)) { d <- d[, 1:N0]; setnames(d, canon); return(d) }

  wide   <- d[has_extra]
  normal <- d[!has_extra][, 1:N0]
  setnames(normal, canon)

  rows_all <- list(); i <- 0
  for (r in seq_len(nrow(wide))) {
    v <- as.character(wide[r])            # 保留中间 NA(即空字段)
    last <- max(which(!is.na(v)))
    v <- v[1:last]
    k <- ceiling(last / N0)
    for (b in seq_len(k)) {
      seg <- v[((b - 1) * N0 + 1):min(b * N0, last)]
      row <- as.list(seg)
      if (length(row) < N0)
        for (jj in (length(row) + 1):N0) row[[jj]] <- NA_character_
      i <- i + 1; rows_all[[i]] <- as.data.table(row)
    }
  }
  fixed <- rbindlist(rows_all, fill = TRUE)
  setnames(fixed, canon)
  rbindlist(list(normal, fixed), use.names = TRUE, fill = TRUE)
}
