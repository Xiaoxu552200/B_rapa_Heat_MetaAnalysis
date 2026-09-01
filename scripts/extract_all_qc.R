#!/usr/bin/env Rscript
# ============================================================================
# 综合提取 fastp 质控指标（LOG + JSON + HTML 三种格式）
# 输出：60个样本的完整质控表格
# ============================================================================

library(tidyverse)
library(jsonlite)
library(rvest)

setwd("/mnt/g/B_rapa_work/B_rapa_Heat_MetaAnalysis")

cat("========================================\n")
cat("综合提取 fastp 质控指标\n")
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
results <- list()

# ============================================================
# 1. 解析 LOG 文件
# ============================================================
parse_log <- function(file_path) {
  tryCatch({
    lines <- read_lines(file_path)
    res <- list(Total_Reads = NA, Pass_Reads = NA, Q20_Rate = NA, 
                Q30_Rate = NA, GC_Content = NA, Read_Length = NA)
    
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
    # 如果所有指标都是 NA，返回 NULL
    if (all(is.na(unlist(res)))) return(NULL)
    res$Source <- "log"
    return(res)
  }, error = function(e) { return(NULL) })
}

# ============================================================
# 2. 解析 JSON 文件
# ============================================================
parse_json <- function(file_path) {
  tryCatch({
    data <- fromJSON(file_path)
    summary <- data$summary
    before <- summary$before_filtering
    after <- summary$after_filtering
    
    list(
      Total_Reads = before$total_reads,
      Pass_Reads = after$total_reads,
      Q20_Rate = after$q20_rate * 100,
      Q30_Rate = after$q30_rate * 100,
      GC_Content = after$gc_content * 100,
      Read_Length = after$read1_mean_length,
      Source = "json"
    )
  }, error = function(e) { return(NULL) })
}

# ============================================================
# 3. 解析 HTML 文件（从表格中提取）
# ============================================================
parse_html <- function(file_path) {
  tryCatch({
    html <- read_html(file_path)
    
    # 方法1：提取所有表格文本
    tables <- html %>% html_nodes("table")
    
    for (table in tables) {
      table_text <- table %>% html_text()
      
      # 提取 Q20/Q30/GC
      q20_match <- str_extract(table_text, "Q20[^%]+%")
      q30_match <- str_extract(table_text, "Q30[^%]+%")
      gc_match <- str_extract(table_text, "GC content[^%]+%")
      
      if (!is.na(q20_match) && !is.na(q30_match)) {
        q20_val <- as.numeric(str_extract(q20_match, "[0-9.]+"))
        q30_val <- as.numeric(str_extract(q30_match, "[0-9.]+"))
        gc_val <- ifelse(!is.na(gc_match), as.numeric(str_extract(gc_match, "[0-9.]+")), NA)
        
        # 提取 reads 数量
        total_match <- str_extract(table_text, "total reads[:\\s]*[0-9,]+")
        pass_match <- str_extract(table_text, "passed filter[:\\s]*[0-9,]+")
        len_match <- str_extract(table_text, "read length[:\\s]*[0-9]+")
        
        total_reads <- ifelse(!is.na(total_match), 
                              as.numeric(str_extract(total_match, "[0-9,]+") %>% str_replace_all(",", "")), 
                              NA)
        pass_reads <- ifelse(!is.na(pass_match), 
                             as.numeric(str_extract(pass_match, "[0-9,]+") %>% str_replace_all(",", "")), 
                             NA)
        read_length <- ifelse(!is.na(len_match), 
                              as.numeric(str_extract(len_match, "[0-9]+")), 
                              NA)
        
        return(list(
          Total_Reads = total_reads,
          Pass_Reads = pass_reads,
          Q20_Rate = q20_val,
          Q30_Rate = q30_val,
          GC_Content = gc_val,
          Read_Length = read_length,
          Source = "html"
        ))
      }
    }
    return(NULL)
  }, error = function(e) { return(NULL) })
}

# ============================================================
# 4. 遍历所有样本
# ============================================================
cat("\n处理样本...\n")

for (s in samples) {
  html_file <- file.path(counts_dir, paste0(s, "_fastp.html"))
  json_file <- file.path(counts_dir, paste0(s, "_fastp.json"))
  log_file <- file.path(counts_dir, paste0(s, "_fastp.log"))
  
  result <- NULL
  source_type <- "unknown"
  
  # 按优先级：JSON > HTML > LOG
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
  
  if (!is.null(result) && !is.na(result$Q20_Rate) && result$Q20_Rate > 0) {
    results[[s]] <- data.frame(
      Sample = s,
      Total_Reads = ifelse(is.na(result$Total_Reads), 0, result$Total_Reads),
      Pass_Reads = ifelse(is.na(result$Pass_Reads), 0, result$Pass_Reads),
      Q20_Rate = round(result$Q20_Rate, 2),
      Q30_Rate = round(result$Q30_Rate, 2),
      GC_Content = round(ifelse(is.na(result$GC_Content), 0, result$GC_Content), 2),
      Read_Length = round(ifelse(is.na(result$Read_Length), 0, result$Read_Length), 2),
      Source = source_type,
      stringsAsFactors = FALSE
    )
    cat("  ✅", s, "(", source_type, ")\n")
  } else {
    cat("  ❌", s, "\n")
  }
}

# ============================================================
# 5. 合并结果
# ============================================================
if (length(results) > 0) {
  final_df <- bind_rows(results)
} else {
  final_df <- data.frame()
}

# ============================================================
# 6. 保存
# ============================================================
output_dir <- "results/QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(final_df, file.path(output_dir, "qc_metrics_60samples_final.csv"))

# 生成报告
sink(file.path(output_dir, "qc_summary_report_final.txt"))
cat("========================================\n")
cat("fastp 质控指标汇总报告\n")
cat("========================================\n\n")
cat("处理时间:", Sys.time(), "\n")
cat("总样本数:", length(samples), "\n")
cat("成功提取:", nrow(final_df), "\n")
cat("缺失:", length(samples) - nrow(final_df), "\n\n")

if (nrow(final_df) > 0) {
  cat("各指标平均值:\n")
  cat("  Total Reads:", format(mean(final_df$Total_Reads), scientific = FALSE, big.mark = ","), "\n")
  cat("  Pass Reads:", format(mean(final_df$Pass_Reads), scientific = FALSE, big.mark = ","), "\n")
  cat("  Q20 Rate (%):", round(mean(final_df$Q20_Rate), 2), "\n")
  cat("  Q30 Rate (%):", round(mean(final_df$Q30_Rate), 2), "\n")
  cat("  GC Content (%):", round(mean(final_df$GC_Content), 2), "\n")
  cat("  Read Length:", round(mean(final_df$Read_Length), 2), "\n")
  
  cat("\n按来源分类:\n")
  table(final_df$Source) %>% print()
}
sink()

cat("\n========================================\n")
cat("✅ 完成！\n")
cat("成功提取:", nrow(final_df), "/", length(samples), "\n")
cat("输出文件: results/QC/qc_metrics_60samples_final.csv\n")
cat("报告: results/QC/qc_summary_report_final.txt\n")
cat("========================================\n")

