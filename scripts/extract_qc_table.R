#!/usr/bin/env Rscript
# ============================================================================
# 从 fastp 报告中提取 Q20、Q30 等质控指标
# 支持 HTML、JSON、LOG 三种格式
# ============================================================================

library(tidyverse)
library(jsonlite)
library(rvest)

setwd("/mnt/g/B_rapa_work/B_rapa_Heat_MetaAnalysis")

cat("========================================\n")
cat("提取 fastp 质控指标\n")
cat("========================================\n")

# 定义60个样本
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
results <- data.frame()

# 解析 JSON 文件
parse_json <- function(file_path) {
  tryCatch({
    data <- fromJSON(file_path)
    summary <- data$summary
    before <- summary$before_filtering
    after <- summary$after_filtering
    
    return(data.frame(
      Total_Reads = before$total_reads,
      Pass_Reads = after$total_reads,
      Q20_Rate = after$q20_rate * 100,
      Q30_Rate = after$q30_rate * 100,
      GC_Content = after$gc_content * 100,
      Read_Length = after$read1_mean_length,
      Source = "json"
    ))
  }, error = function(e) { return(NULL) })
}

# 解析 LOG 文件
parse_log <- function(file_path) {
  tryCatch({
    lines <- read_lines(file_path)
    
    total_reads <- NA
    pass_reads <- NA
    q20_rate <- NA
    q30_rate <- NA
    gc_content <- NA
    read_length <- NA
    
    for (line in lines) {
      if (grepl("total reads:", line)) {
        total_reads <- as.numeric(str_extract(line, "[0-9,]+") %>% str_replace_all(",", ""))
      }
      if (grepl("passed filters reads:", line)) {
        pass_reads <- as.numeric(str_extract(line, "[0-9,]+") %>% str_replace_all(",", ""))
      }
      if (grepl("Q20 rate:", line)) {
        q20_rate <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("Q30 rate:", line)) {
        q30_rate <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("GC content:", line)) {
        gc_content <- as.numeric(str_extract(line, "[0-9.]+"))
      }
      if (grepl("length of reads:", line)) {
        read_length <- as.numeric(str_extract(line, "[0-9.]+"))
      }
    }
    
    return(data.frame(
      Total_Reads = total_reads,
      Pass_Reads = pass_reads,
      Q20_Rate = q20_rate,
      Q30_Rate = q30_rate,
      GC_Content = gc_content,
      Read_Length = read_length,
      Source = "log"
    ))
  }, error = function(e) { return(NULL) })
}

# 解析 HTML 文件（尝试提取嵌入的 JSON）
parse_html <- function(file_path) {
  tryCatch({
    html <- read_html(file_path)
    json_text <- html %>% html_nodes("script") %>% html_text() %>% paste(collapse = "\n")
    json_match <- str_extract(json_text, '\\{"passed_filter".*\\}')
    if (!is.na(json_match)) {
      data <- fromJSON(json_match)
      summary <- data$summary
      before <- summary$before_filtering
      after <- summary$after_filtering
      return(data.frame(
        Total_Reads = before$total_reads,
        Pass_Reads = after$total_reads,
        Q20_Rate = after$q20_rate * 100,
        Q30_Rate = after$q30_rate * 100,
        GC_Content = after$gc_content * 100,
        Read_Length = after$read1_mean_length,
        Source = "html"
      ))
    }
    return(NULL)
  }, error = function(e) { return(NULL) })
}

# 遍历所有样本
cat("\n处理样本...\n")
for (s in samples) {
  html_file <- file.path(counts_dir, paste0(s, "_fastp.html"))
  json_file <- file.path(counts_dir, paste0(s, "_fastp.json"))
  log_file <- file.path(counts_dir, paste0(s, "_fastp.log"))
  
  result <- NULL
  source_type <- "unknown"
  
  if (file.exists(json_file)) {
    result <- parse_json(json_file)
    source_type <- "json"
  } else if (file.exists(html_file)) {
    result <- parse_html(html_file)
    source_type <- "html"
  } else if (file.exists(log_file)) {
    result <- parse_log(log_file)
    source_type <- "log"
  }
  
  if (!is.null(result) && nrow(result) > 0) {
    row <- data.frame(
      Sample = s,
      Total_Reads = result$Total_Reads,
      Pass_Reads = result$Pass_Reads,
      Q20_Rate = round(result$Q20_Rate, 2),
      Q30_Rate = round(result$Q30_Rate, 2),
      GC_Content = round(result$GC_Content, 2),
      Read_Length = round(result$Read_Length, 2),
      Source = source_type,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, row)
    cat("  ✅", s, "(", source_type, ")\n")
  } else {
    cat("  ❌", s, "\n")
  }
}

# 保存表格
output_dir <- "results/QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(results, file.path(output_dir, "qc_metrics_60samples.csv"))

# 生成摘要报告
sink(file.path(output_dir, "qc_summary_report.txt"))
cat("========================================\n")
cat("fastp 质控指标汇总报告\n")
cat("========================================\n\n")
cat("总样本数:", nrow(results), "/", length(samples), "\n")
cat("缺失样本数:", length(samples) - nrow(results), "\n\n")

if (nrow(results) > 0) {
  cat("质控指标统计 (平均值):\n")
  cat("  Total Reads:", format(mean(results$Total_Reads), scientific = FALSE, big.mark = ","), "\n")
  cat("  Pass Reads:", format(mean(results$Pass_Reads), scientific = FALSE, big.mark = ","), "\n")
  cat("  Q20 Rate (%):", round(mean(results$Q20_Rate), 2), "\n")
  cat("  Q30 Rate (%):", round(mean(results$Q30_Rate), 2), "\n")
  cat("  GC Content (%):", round(mean(results$GC_Content), 2), "\n")
  cat("  Read Length:", round(mean(results$Read_Length), 2), "\n")
  
  cat("\n按来源分类:\n")
  table(results$Source) %>% print()
}

sink()

cat("\n✅ 完成！\n")
cat("输出文件: results/QC/qc_metrics_60samples.csv\n")
cat("报告: results/QC/qc_summary_report.txt\n")

