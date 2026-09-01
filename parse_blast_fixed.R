#!/usr/bin/env Rscript
# ============================================================================
# 解析BLAST结果并生成拟南芥同源ID
# 修正版：直接读取 dates/PPI 目录下的文件
# ============================================================================

library(tidyverse)
library(ggplot2)

# 设置工作目录到 PPI 文件夹
# setwd 已移除，使用相对路径

cat("========================================\n")
cat("步骤: 解析BLAST结果\n")
cat("========================================\n")

# 输出目录（当前目录）
output_dir <- "results/homology"
figures_dir <- file.path(output_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# 查找BLAST结果文件
blast_file <- file.path("data", "blast_results.txt")

if(!file.exists(blast_file)) {
  cat("错误: BLAST结果文件不存在\n")
  quit()
}

# 读取BLAST结果
blast_results <- read.table(blast_file, 
                            stringsAsFactors = FALSE,
                            sep = "\t",
                            fill = TRUE)

colnames(blast_results) <- c("query_id", "subject_id", "identity", 
                             "alignment_length", "mismatches", "gap_opens",
                             "q_start", "q_end", "s_start", "s_end",
                             "evalue", "bit_score")

cat("BLAST结果行数:", nrow(blast_results), "\n")

# 提取基础基因ID
extract_base_id <- function(query_id) {
  base_id <- str_extract(query_id, "^Bra[^\\s\\.]+")
  if(is.na(base_id)) {
    base_id <- str_extract(query_id, "^Bra[A-Za-z0-9]+g[0-9]+")
  }
  return(base_id)
}

blast_results$query_base <- sapply(blast_results$query_id, extract_base_id)

# 读取方向信息
up_genes <- read_lines("data/key_genes_up.txt")
down_genes <- read_lines("data/key_genes_down.txt")

# 构建映射结果
mapping_results <- data.frame(
  Brapa_Gene = blast_results$query_base,
  Ath_Gene = blast_results$subject_id,
  Identity = blast_results$identity,
  Alignment_Length = blast_results$alignment_length,
  E_value = blast_results$evalue,
  Bit_Score = blast_results$bit_score,
  Direction = ifelse(blast_results$query_base %in% up_genes, "Up", "Down"),
  stringsAsFactors = FALSE
)

# 去重：每个白菜基因保留E-value最低的匹配
unique_genes <- unique(mapping_results$Brapa_Gene)
mapping_unique <- data.frame()

for(gene in unique_genes) {
  sub <- mapping_results[mapping_results$Brapa_Gene == gene, ]
  sub <- sub[order(sub$E_value), ]
  mapping_unique <- rbind(mapping_unique, sub[1, ])
}

mapping_results <- mapping_unique

# 统计结果
total_genes <- length(c(up_genes, down_genes))
mapped_count <- nrow(mapping_results)

cat("\n========== 映射统计 ==========\n")
cat("总关键基因数:", total_genes, "\n")
cat("成功映射:", mapped_count, "\n")
cat("映射率:", round(mapped_count/total_genes*100, 2), "%\n")

# 按方向统计
up_mapped <- mapping_results[mapping_results$Direction == "Up", ]
down_mapped <- mapping_results[mapping_results$Direction == "Down", ]

cat("\n上调基因映射:", nrow(up_mapped), "/", length(up_genes), "\n")
cat("下调基因映射:", nrow(down_mapped), "/", length(down_genes), "\n")

# 提取拟南芥同源基因
ath_up <- unique(mapping_results$Ath_Gene[mapping_results$Direction == "Up"])
ath_down <- unique(mapping_results$Ath_Gene[mapping_results$Direction == "Down"])

cat("\n拟南芥同源基因:\n")
cat("  上调同源基因:", length(ath_up), "个\n")
cat("  下调同源基因:", length(ath_down), "个\n")
cat("  总同源基因:", length(unique(c(ath_up, ath_down))), "个\n")

# 保存结果
write_csv(mapping_results, "homology_mapping_results.csv")
write_lines(ath_up, "ath_homologs_up.txt")
write_lines(ath_down, "ath_homologs_down.txt")
write_lines(c(ath_up, ath_down), "ath_homologs_all.txt")

# 未映射基因
mapped_genes <- unique(mapping_results$Brapa_Gene)
unmapped <- setdiff(c(up_genes, down_genes), mapped_genes)
if(length(unmapped) > 0) {
  write_lines(unmapped, "unmapped_genes.txt")
  cat("\n未映射基因数:", length(unmapped), "\n")
}

# ============================================================================
# 可视化
# ============================================================================
if(nrow(mapping_results) > 0) {
  p1 <- ggplot(mapping_results, aes(x = -log10(E_value), fill = Direction)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    labs(title = "BLAST E-value Distribution",
         x = "-log10(E-value)", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(figures_dir, "evalue_distribution.pdf"), p1, width = 8, height = 6)
  
  p2 <- ggplot(mapping_results, aes(x = Identity, fill = Direction)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    labs(title = "Sequence Identity Distribution",
         x = "Identity (%)", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(figures_dir, "identity_distribution.pdf"), p2, width = 8, height = 6)
  
  p3 <- ggplot(mapping_results, aes(x = Identity, y = -log10(E_value), color = Direction)) +
    geom_point(alpha = 0.6, size = 2) +
    labs(title = "BLAST Results: E-value vs Identity",
         x = "Identity (%)", y = "-log10(E-value)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(figures_dir, "evalue_vs_identity.pdf"), p3, width = 8, height = 6)
}

# ============================================================================
# 生成报告
# ============================================================================
sink("homology_mapping_report.txt")
cat("==================================================================\n")
cat("同源映射报告: 白菜(B. rapa) -> 拟南芥(A. thaliana)\n")
cat("==================================================================\n\n")

cat("【输入基因】\n")
cat("  总关键基因:", total_genes, "\n")
cat("  上调基因:", length(up_genes), "\n")
cat("  下调基因:", length(down_genes), "\n\n")

cat("【映射结果】\n")
cat("  成功映射:", mapped_count, "个 (", round(mapped_count/total_genes*100, 2), "%)\n")
cat("  上调映射:", nrow(up_mapped), "个\n")
cat("  下调映射:", nrow(down_mapped), "个\n\n")

cat("【拟南芥同源基因】\n")
cat("  上调同源物:", length(ath_up), "个\n")
cat("  下调同源物:", length(ath_down), "个\n")
cat("  总同源物:", length(unique(c(ath_up, ath_down))), "个\n\n")

cat("【最佳匹配（E-value最小前20）】\n")
best_for_report <- mapping_results[order(mapping_results$E_value), ]
best_for_report <- head(best_for_report[, c("Brapa_Gene", "Ath_Gene", "Identity", "E_value", "Direction")], 20)
print(best_for_report)

sink()

cat("\n")
cat("==================================================================\n")
cat("✅ 同源映射完成！\n")
cat("==================================================================\n")
cat("\n输出文件:\n")
cat("  homology_mapping_results.csv - 完整映射结果\n")
cat("  ath_homologs_up.txt - 上调拟南芥同源ID\n")
cat("  ath_homologs_down.txt - 下调拟南芥同源ID\n")
cat("  ath_homologs_all.txt - 全部拟南芥同源ID\n")
cat("  homology_mapping_report.txt - 详细报告\n")
cat("  figures/*.pdf - 可视化图表\n")
cat("==================================================================\n")

