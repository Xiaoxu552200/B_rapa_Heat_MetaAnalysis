#!/usr/bin/env Rscript
# ============================================================================
# PCA分析：批次校正（ComBat）前后对比
# 输出：校正前后PCA图（TIFF/PNG/PDF）、校正前后矩阵、R会话信息
# 运行方式：在项目根目录下执行 Rscript scripts/PCA_analysis_clean.R
# ============================================================================

# 设置随机种子（确保可复现）
set.seed(12345)

# 加载必要的包
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

cat("========================================\n")
cat("PCA分析：批次校正前后对比\n")
cat("========================================\n")

# ============================================================================
# 1. 设置路径（全部使用相对路径，基于项目根目录）
# ============================================================================
# 获取项目根目录（假设运行脚本时当前工作目录即为项目根目录）
project_root <- getwd()
cat("项目根目录:", project_root, "\n")

# 定义数据目录和输出目录
data_dir <- file.path(project_root, "data")
counts_dir <- file.path(data_dir, "counts")          # 存放所有 *_counts.txt 文件
output_dir <- file.path(project_root, "results", "PCA")

# 创建输出目录（如果不存在）
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)

# ============================================================================
# 2. 读取样本信息（metadata）
# ============================================================================
cat("\n[1/6] 读取样本信息...\n")

# 请将您的 metadata 文件命名为 "Final_60_Samples_Metadata.csv" 并放在 data/ 目录下
metadata_file <- file.path(data_dir, "Final_60_Samples_Metadata.csv")
if (!file.exists(metadata_file)) {
  stop("错误: 未找到 ", metadata_file, "，请确保 metadata 文件已放入 data/ 目录")
}
csv_keep <- read.csv(metadata_file, stringsAsFactors = FALSE)

cat("样本数:", nrow(csv_keep), "\n")
cat("数据集分布:\n")
print(table(csv_keep$dataset))
cat("\n处理组分布:\n")
print(table(csv_keep$group))

# 转换为 condition（用于分析）
csv_keep$condition <- ifelse(csv_keep$group == "treatment", "Stress", "Normal")

# ============================================================================
# 3. 读取原始 counts 矩阵（从 counts_dir 下的所有 *_counts.txt 文件）
# ============================================================================
cat("\n[2/6] 读取原始counts矩阵...\n")

all_counts <- list.files(path = counts_dir, pattern = "*_counts.txt", full.names = TRUE)
if (length(all_counts) == 0) {
  stop("错误: counts_dir 下没有找到 *_counts.txt 文件，请将 counts 文件放入 data/counts/ 目录")
}

# 确保 counts 文件与 metadata 中的 sample_id 匹配
counts_ids <- gsub("_counts.txt", "", basename(all_counts))
csv_keep <- csv_keep[csv_keep$sample_id %in% counts_ids, ]
cat("匹配到的样本数:", nrow(csv_keep), "\n")

# 读取第一个文件获取基因列表
first <- read.table(all_counts[1], header = TRUE)
raw_count_mat <- matrix(0, nrow = nrow(first), ncol = nrow(csv_keep))
rownames(raw_count_mat) <- first$Geneid
colnames(raw_count_mat) <- csv_keep$sample_id

for (i in 1:nrow(csv_keep)) {
  f <- all_counts[grep(csv_keep$sample_id[i], basename(all_counts))]
  if (length(f) > 0) {
    d <- read.table(f[1], header = TRUE)
    raw_count_mat[, i] <- d[match(rownames(raw_count_mat), d$Geneid), ncol(d)]
    raw_count_mat[is.na(raw_count_mat)] <- 0
  }
}

# 过滤低表达基因（总和为0的基因删除）
raw_count_mat <- raw_count_mat[rowSums(raw_count_mat) > 0, ]
cat("原始counts矩阵维度:", dim(raw_count_mat), "\n")

# 保存原始 counts 矩阵
saveRDS(raw_count_mat, file.path(output_dir, "tables", "raw_count_matrix.rds"))
write.csv(as.data.frame(raw_count_mat), 
          file.path(output_dir, "tables", "raw_count_matrix.csv"), 
          row.names = TRUE)
cat("✓ 原始counts矩阵已保存\n")

# ============================================================================
# 4. 原始数据 PCA（无批次校正）
# ============================================================================
cat("\n[3/6] 原始数据PCA...\n")

# 转换为 log2CPM
logcpm_raw <- cpm(raw_count_mat, log = TRUE, prior.count = 1)
logcpm_raw_t <- t(logcpm_raw)

# 保存校正前的 log2CPM 矩阵
saveRDS(logcpm_raw, file.path(output_dir, "tables", "log2cpm_before_correction.rds"))
write.csv(as.data.frame(logcpm_raw), 
          file.path(output_dir, "tables", "log2cpm_before_correction.csv"), 
          row.names = TRUE)
cat("✓ 校正前log2CPM矩阵已保存\n")

# PCA 分析
pca_raw <- prcomp(logcpm_raw_t, scale. = TRUE)
pca_raw_df <- data.frame(
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  dataset = csv_keep$dataset,
  condition = csv_keep$condition
)

var_pc1_raw <- round(summary(pca_raw)$importance[2, 1] * 100, 1)
var_pc2_raw <- round(summary(pca_raw)$importance[2, 2] * 100, 1)

# 定义形状（根据实际数据集名称调整）
shape_values <- c(
  "PRJNA872510"  = 16,
  "PRJNA1030162" = 17,
  "PRJNA646007"  = 15,
  "PRJNA663233"  = 18
)

p_raw <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, color = condition, shape = dataset)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = c("Normal" = "#2E5A9C", "Stress" = "#C93312")) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("PC1 (", var_pc1_raw, "%)"),
       y = paste0("PC2 (", var_pc2_raw, "%)"),
       title = "Before Batch Correction (Raw Data)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 保存 PCA 结果
saveRDS(pca_raw, file.path(output_dir, "tables", "pca_before_correction.rds"))
write.csv(pca_raw$x, 
          file.path(output_dir, "tables", "pca_before_correction_scores.csv"), 
          row.names = TRUE)

# ============================================================================
# 5. ComBat 批次校正
# ============================================================================
cat("\n[4/6] ComBat批次校正...\n")

logcpm_for_combat <- cpm(raw_count_mat, log = TRUE, prior.count = 1)
batch <- csv_keep$dataset
mod <- model.matrix(~ condition, data = csv_keep)

corrected_mat <- ComBat(dat = logcpm_for_combat, batch = batch, mod = mod, par.prior = TRUE)
cat("校正后矩阵维度:", dim(corrected_mat), "\n")

# 保存校正后的矩阵
saveRDS(corrected_mat, file.path(output_dir, "tables", "combat_corrected_matrix.rds"))
write.csv(as.data.frame(corrected_mat), 
          file.path(output_dir, "tables", "combat_corrected_matrix.csv"), 
          row.names = TRUE)
cat("✓ ComBat校正矩阵已保存\n")

# 转置用于 PCA
corrected_t <- t(corrected_mat)

# ============================================================================
# 6. 校正后 PCA
# ============================================================================
cat("\n[5/6] 校正后PCA...\n")

pca_corrected <- prcomp(corrected_t, scale. = TRUE)
pca_corrected_df <- data.frame(
  PC1 = pca_corrected$x[, 1],
  PC2 = pca_corrected$x[, 2],
  dataset = csv_keep$dataset,
  condition = csv_keep$condition
)

var_pc1_corr <- round(summary(pca_corrected)$importance[2, 1] * 100, 1)
var_pc2_corr <- round(summary(pca_corrected)$importance[2, 2] * 100, 1)

p_corrected <- ggplot(pca_corrected_df, aes(x = PC1, y = PC2, color = condition, shape = dataset)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = c("Normal" = "#2E5A9C", "Stress" = "#C93312")) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("PC1 (", var_pc1_corr, "%)"),
       y = paste0("PC2 (", var_pc2_corr, "%)"),
       title = "After Batch Correction (ComBat)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 保存校正后 PCA 结果
saveRDS(pca_corrected, file.path(output_dir, "tables", "pca_after_correction.rds"))
write.csv(pca_corrected$x, 
          file.path(output_dir, "tables", "pca_after_correction_scores.csv"), 
          row.names = TRUE)

# ============================================================================
# 7. 保存图片（TIFF、PNG、PDF）
# ============================================================================
cat("\n[6/6] 保存图片...\n")

fig_dir <- file.path(output_dir, "figures")

# TIFF 格式（期刊推荐，300 DPI）
tiff(file.path(fig_dir, "PCA_Before_Correction.tiff"), 
     width = 7, height = 6, units = "in", res = 300, compression = "lzw")
print(p_raw)
dev.off()

tiff(file.path(fig_dir, "PCA_After_Correction.tiff"), 
     width = 7, height = 6, units = "in", res = 300, compression = "lzw")
print(p_corrected)
dev.off()

# 组合图
combined_plot <- p_raw + p_corrected +
  plot_annotation(
    title = "PCA Comparison: Before vs After ComBat Batch Correction",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

tiff(file.path(fig_dir, "PCA_Comparison_Combined.tiff"), 
     width = 14, height = 6, units = "in", res = 300, compression = "lzw")
print(combined_plot)
dev.off()

# PNG 格式（预览）
png(file.path(fig_dir, "PCA_Before_Correction.png"), 
    width = 7, height = 6, units = "in", res = 300)
print(p_raw)
dev.off()

png(file.path(fig_dir, "PCA_After_Correction.png"), 
    width = 7, height = 6, units = "in", res = 300)
print(p_corrected)
dev.off()

png(file.path(fig_dir, "PCA_Comparison_Combined.png"), 
    width = 14, height = 6, units = "in", res = 300)
print(combined_plot)
dev.off()

# PDF 矢量图（可编辑）
pdf(file.path(fig_dir, "PCA_Comparison_Combined.pdf"), 
    width = 14, height = 6)
print(combined_plot)
dev.off()

cat("\n✓ 图片已保存（TIFF、PNG、PDF）\n")

# ============================================================================
# 8. 保存样本信息和会话信息
# ============================================================================
write.csv(csv_keep, 
          file.path(output_dir, "tables", "sample_metadata.csv"), 
          row.names = FALSE)

sink(file.path(output_dir, "session_info.txt"))
cat("R Session Info\n")
cat("==============\n\n")
sessionInfo()
cat("\n\nPackages used:\n")
cat("==============\n")
print(loadedNamespaces())
sink()

# ============================================================================
# 9. 输出汇总
# ============================================================================
cat("\n")
cat("================================================================================\n")
cat("✅ PCA分析完成！\n")
cat("================================================================================\n")
cat("输出目录:", output_dir, "\n")
cat("\n矩阵文件 (tables/):\n")
cat("  - raw_count_matrix.rds/csv\n")
cat("  - log2cpm_before_correction.rds/csv\n")
cat("  - combat_corrected_matrix.rds/csv\n")
cat("  - pca_before_correction.rds/csv\n")
cat("  - pca_after_correction.rds/csv\n")
cat("  - sample_metadata.csv\n")
cat("\n图片文件 (figures/):\n")
cat("  - PCA_Before_Correction.tiff/png\n")
cat("  - PCA_After_Correction.tiff/png\n")
cat("  - PCA_Comparison_Combined.tiff/png/pdf\n")
cat("\n其他文件:\n")
cat("  - session_info.txt\n")
cat("================================================================================\n")
