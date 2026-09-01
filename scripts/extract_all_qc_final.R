#!/usr/bin/env Rscript
# ============================================================================
# 最终版：提取所有60个样本的Q20/Q30质控指标
# ============================================================================

library(tidyverse)
library(jsonlite)
library(rvest)

setwd("/mnt/g/B_rapa_work/B_rapa_Heat_MetaAnalysis")

cat("========================================\n")
cat("提取所有60个样本的质控指标\n")
cat("========================================\n")

samples <- c(
  "SRR12632445", "SRR12632457", "SRR12632456", "SRR12632437", "SRR12632436",
  "SRR12632435", "SRR12632407", "SRR12632406", "SRR12632405", "SRR12632430",
  "SRR12632429", "SRR12632428", "SRR12632446", "SRR12632444", "SRR12632447",
  "SRR12632404", "SRR12632455", "SRR12632454", "SRR12632453", "SRR12632452",
  "SRR12632451", "SRR12632424", "SRR12632422", "SRR12632421", "SRR12632443",
  "SRR12632442", "SRR12632441", "SRR26447755", "SRR26447754", "SRR26447745",
  "SRR26447751", "SRR26447750", "SRR26447749", "SRR26447741", "SRR26447740",
  "SRR26447739", "SRR26447744", "SRR26447743", "SRR26447742", "SRR26447748",
  "SRR26447747", "SRR26447746", "SRR26447738", "SRR26447753", "SRR26447752",
  "SRR21227071", "SRR21227073", "SRR21227074", "SRR21227072", "SRR21227069",
  "SRR21227070", "SRR21227075", "SRR21227076", "SRR21227065",
  "SRR12214241", "SRR12214243", "SRR12214249", "SRR12214245", "SRR12214247",
  "SRR12214251"
)

counts_dir <- "data/counts"
results <- list()

# ============================================================
# 解析 JSON（最可靠）
# ============================================================
parse_json <- function(f) {
  tryCatch({
    d <- fromJSON(f)
    s <- d$summary
    before <- s$before_filtering
    after <- s$after_filtering
    return(list(
      Total_Reads = before$total_reads,
      Pass_Reads = after$total_reads,
      Q20_Rate = after$q20_rate * 100,
      Q30_Rate = after$q30_rate * 100,
      GC_Content = after$gc_content * 100,
      Read_Length = after$read1_mean_length
    ))
  }, error = function(e) return(NULL))
}

# ============================================================
# 解析 LOG
# ============================================================
parse_log <- function(f) {
  tryCatch({
    lines <- read_lines(f)
    res <- list(Total_Reads=NA, Pass_Reads=NA, Q20_Rate=NA, Q30_Rate=NA, GC_Content=NA, Read_Length=NA)
    for (line in lines) {
      if (grepl("total reads:", line)) {
        res$Total_Reads <- as.numeric(str_extract(line, "[0-9,]+") %>% str_replace_all(",", ""))
      }
      if (grepl("passed filters reads:", line)) {
        res$Pass_Reads <- as.numeric(str_extract(line, "[0-9,]+") %>% str_replace_all(",", ""))
      }
      if (grepl("Q20 rate:", line)) {
        res$Q20_Rate <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("Q30 rate:", line)) {
        res$Q30_Rate <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("GC content:", line)) {
        res$GC_Content <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("length of reads:", line)) {
        res$Read_Length <- as.numeric(str_extract(line, "[0-9.]+"))
      }
    }
    return(res)
  }, error = function(e) return(NULL))
}

# ============================================================
# 解析 HTML（从表格提取）
# ============================================================
parse_html <- function(f) {
  tryCatch({
    html <- read_html(f)
    text <- html %>% html_text()
    
    # 找到 After filtering 部分
    after_idx <- str_locate(text, "After filtering")[2]
    if (is.na(after_idx)) after_idx <- 1
    after_text <- substr(text, after_idx, nchar(text))
    
    # 提取数值（含单位转换）
    extract_val <- function(label, txt) {
      pattern <- paste0(label, ":.*?>([0-9.]+)")
      m <- str_extract(txt, pattern)
      if (!is.na(m)) {
        val <- as.numeric(str_extract(m, "[0-9.]+"))
        # 检查是否有 M/G 单位
        if (grepl("M", m)) val <- val * 1000000
        return(val)
      }
      return(NA)
    }
    
    extract_pct <- function(label, txt) {
      pattern <- paste0(label, ":.*?\\(([0-9.]+)%\\)")
      m <- str_extract(txt, pattern)
      if (!is.na(m)) return(as.numeric(str_extract(m, "[0-9.]+")))
      return(NA)
    }
    
    total_reads <- extract_val("total reads", after_text)
    pass_reads <- extract_val("reads passed filters", after_text)
    q20_rate <- extract_pct("Q20 bases", after_text)
    q30_rate <- extract_pct("Q30 bases", after_text)
    gc_content <- extract_val("GC content", after_text)
    
    # 如果GC没提取到，从before部分提取
    if (is.na(gc_content)) {
      before_idx <- str_locate(text, "Before filtering")[2]
      if (!is.na(before_idx)) {
        before_text <- substr(text, 1, before_idx)
        gc_content <- extract_val("GC content", before_text)
      }
    }
    
    return(list(
      Total_Reads = total_reads,
      Pass_Reads = pass_reads,
      Q20_Rate = q20_rate,
      Q30_Rate = q30_rate,
      GC_Content = gc_content,
      Read_Length = NA
    ))
  }, error = function(e) return(NULL))
}

# ============================================================
# 遍历所有样本
# ============================================================
cat("\n处理60个样本...\n")

for (s in samples) {
  html_f <- file.path(counts_dir, paste0(s, "_fastp.html"))
  json_f <- file.path(counts_dir, paste0(s, "_fastp.json"))
  log_f <- file.path(counts_dir, paste0(s, "_fastp.log"))
  
  result <- NULL
  src <- "unknown"
  
  # 优先级：JSON > HTML > LOG
  if (file.exists(json_f)) {
    result <- parse_json(json_f)
    src <- "json"
  } else if (file.exists(html_f)) {
    result <- parse_html(html_f)
    src <- "html"
  } else if (file.exists(log_f)) {
    result <- parse_log(log_f)
    src <- "log"
  }
  
  if (!is.null(result) && !is.na(result$Q20_Rate) && result$Q20_Rate > 0) {
    row <- data.frame(
      Sample = s,
      Total_Reads = ifelse(is.na(result$Total_Reads), 0, result$Total_Reads),
      Pass_Reads = ifelse(is.na(result$Pass_Reads), 0, result$Pass_Reads),
      Q20_Rate = round(result$Q20_Rate, 2),
      Q30_Rate = round(result$Q30_Rate, 2),
      GC_Content = round(ifelse(is.na(result$GC_Content), 0, result$GC_Content), 2),
      Read_Length = round(ifelse(is.na(result$Read_Length), 0, result$Read_Length), 2),
      Source = src,
      stringsAsFactors = FALSE
    )
    results[[s]] <- row
    cat("  ✅", s, "(", src, ")\n")
  } else {
    cat("  ❌", s, "\n")
  }
}

# ============================================================
# 保存
# ============================================================
final_df <- bind_rows(results)

output_dir <- "results/QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(final_df, file.path(output_dir, "qc_metrics_all_60_final.csv"))

cat("\n========================================\n")
cat("✅ 完成！\n")
cat("成功提取:", nrow(final_df), "/", length(samples), "\n")
cat("输出文件: results/QC/qc_metrics_all_60_final.csv\n")
cat("========================================\n")

