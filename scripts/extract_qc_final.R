#!/usr/bin/env Rscript
# ============================================================================
# 最终版：从 HTML 表格中提取 QC 数据
# ============================================================================

library(tidyverse)
library(rvest)
library(jsonlite)

setwd("/mnt/g/B_rapa_work/B_rapa_Heat_MetaAnalysis")

cat("========================================\n")
cat("提取所有60个样本的质控指标 (最终版)\n")
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
# 解析 HTML（使用 rvest 直接提取表格）
# ============================================================
parse_html <- function(f) {
  tryCatch({
    html <- read_html(f)
    
    # 提取 After filtering 表格中的数据
    # 方法：找到包含 "After filtering" 的 div，然后提取表格
    after_div <- html %>% html_nodes("div#after_filtering_summary")
    if (length(after_div) == 0) {
      # 尝试另一种方式：查找包含 "After filtering" 的节点
      all_divs <- html %>% html_nodes("div")
      for (d in all_divs) {
        text <- d %>% html_text()
        if (grepl("After filtering", text) && grepl("total reads", text)) {
          # 这个 div 包含我们要的数据
          tbl <- d %>% html_node("table")
          if (!is.na(tbl)) {
            rows <- tbl %>% html_nodes("tr")
            total_reads <- NA
            q20 <- NA
            q30 <- NA
            gc <- NA
            pass_reads <- NA
            
            for (row in rows) {
              cells <- row %>% html_nodes("td") %>% html_text(trim = TRUE)
              if (length(cells) >= 2) {
                label <- cells[1]
                value <- cells[2]
                
                if (grepl("total reads", label)) {
                  num <- as.numeric(str_extract(value, "[0-9.]+"))
                  if (grepl("M", value)) num <- num * 1000000
                  if (grepl("G", value)) num <- num * 1000000000
                  total_reads <- num
                }
                if (grepl("Q20 bases", label)) {
                  num <- as.numeric(str_extract(value, "[0-9.]+"))
                  if (grepl("G", value)) num <- num * 1000000000
                  q20 <- num
                }
                if (grepl("Q30 bases", label)) {
                  num <- as.numeric(str_extract(value, "[0-9.]+"))
                  if (grepl("G", value)) num <- num * 1000000000
                  q30 <- num
                }
                if (grepl("GC content", label)) {
                  gc <- as.numeric(str_extract(value, "[0-9.]+"))
                }
                if (grepl("reads passed filters", label)) {
                  num <- as.numeric(str_extract(value, "[0-9.]+"))
                  if (grepl("M", value)) num <- num * 1000000
                  if (grepl("G", value)) num <- num * 1000000000
                  pass_reads <- num
                }
              }
            }
            
            # 提取 Q20 和 Q30 百分比
            q20_pct <- NA
            q30_pct <- NA
            if (!is.na(q20)) {
              # 计算 Q20 百分比：Q20 bases / total bases
              # 或者直接从文本中提取百分比
              q20_text <- html %>% html_text()
              q20_match <- str_extract(q20_text, "Q20 bases:.*?\\(([0-9.]+)%\\)")
              if (!is.na(q20_match)) {
                q20_pct <- as.numeric(str_extract(q20_match, "[0-9.]+"))
              }
              q30_match <- str_extract(q20_text, "Q30 bases:.*?\\(([0-9.]+)%\\)")
              if (!is.na(q30_match)) {
                q30_pct <- as.numeric(str_extract(q30_match, "[0-9.]+"))
              }
            }
            
            return(list(
              Total_Reads = total_reads,
              Pass_Reads = pass_reads,
              Q20_Rate = q20_pct,
              Q30_Rate = q30_pct,
              GC_Content = gc,
              Read_Length = NA
            ))
          }
        }
      }
    }
    return(NULL)
  }, error = function(e) { return(NULL) })
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
# 解析 JSON
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
# 主循环
# ============================================================
cat("\n处理60个样本...\n")

for (s in samples) {
  html_f <- file.path(counts_dir, paste0(s, "_fastp.html"))
  json_f <- file.path(counts_dir, paste0(s, "_fastp.json"))
  log_f <- file.path(counts_dir, paste0(s, "_fastp.log"))
  
  result <- NULL
  src <- "unknown"
  
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
write_csv(final_df, file.path(output_dir, "qc_metrics_60_final.csv"))

cat("\n========================================\n")
cat("✅ 完成！成功提取:", nrow(final_df), "/", length(samples), "\n")
cat("输出文件: results/QC/qc_metrics_60_final.csv\n")
cat("========================================\n")

