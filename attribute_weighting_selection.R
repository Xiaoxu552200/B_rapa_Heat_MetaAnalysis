#!/usr/bin/env Rscript
# ============================================================================
# 属性加权筛选核心基因
# 功能：从201个高置信度关键基因中，通过四种属性加权方法筛选核心基因
# 方法：信息增益(IG)、增益比率(GR)、基尼指数(Gini)、Relief算法
# 筛选标准：累计权重 >= 3.0
# 输出：33个核心基因列表 + 完整权重得分表
# ============================================================================

library(tidyverse)
library(CORElearn)

cat("========================================\n")
cat("属性加权筛选核心基因\n")
cat("========================================\n")

# ============================================================================
# 1. 设置路径
# ============================================================================
base_dir <- getwd()

# 输出目录
out_dir <- file.path(base_dir, "results", "ML")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n[1/5] 输出目录:", out_dir, "\n")

# ============================================================================
# 2. 读取数据
# ============================================================================
cat("\n[2/5] 读取数据...\n")

expr <- readRDS("data/combat_corrected_matrix.rds")
sample_info <- read.csv("data/metadata.csv")
y <- ifelse(sample_info$group == "treatment", 1, 0)
cat("  样本数:", length(y), "\n")

# 读取201个关键基因
up_genes <- read_lines("data/key_genes_up.txt")
down_genes <- read_lines("data/key_genes_down.txt")
up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]
all_genes <- unique(c(up_genes, down_genes))
all_genes <- intersect(all_genes, rownames(expr))

cat("  上调基因:", length(up_genes), "\n")
cat("  下调基因:", length(down_genes), "\n")
cat("  总关键基因:", length(all_genes), "\n")

# ============================================================================
# 3. 属性加权计算
# ============================================================================
cat("\n[3/5] 计算四种属性权重...\n")

# 提取表达矩阵并标准化
X <- t(expr[all_genes, ])
X_scaled <- scale(X)
X_df <- as.data.frame(X_scaled)
X_df$target <- as.factor(y)

# 四种属性加权方法
cat("  计算信息增益 (Information Gain)...\n")
ig_weights <- attrEval(target ~ ., data = X_df, estimator = "InfGain")

cat("  计算增益比率 (Gain Ratio)...\n")
gr_weights <- attrEval(target ~ ., data = X_df, estimator = "GainRatio")

cat("  计算基尼指数 (Gini Index)...\n")
gini_weights <- attrEval(target ~ ., data = X_df, estimator = "Gini")

cat("  计算Relief算法...\n")
relief_weights <- attrEval(target ~ ., data = X_df, estimator = "ReliefFequalK")

# ============================================================================
# 4. 归一化并计算累计权重
# ============================================================================
cat("\n[4/5] 归一化并计算累计权重...\n")

normalize_01 <- function(x) {
  if (max(x) - min(x) == 0) return(rep(0, length(x)))
  return((x - min(x)) / (max(x) - min(x)))
}

weight_df <- data.frame(
  Gene = all_genes,
  IG = normalize_01(ig_weights[all_genes]),
  GainRatio = normalize_01(gr_weights[all_genes]),
  Gini = normalize_01(gini_weights[all_genes]),
  Relief = normalize_01(relief_weights[all_genes])
)

weight_df$Cumulative_Weight <- weight_df$IG + weight_df$GainRatio + 
                               weight_df$Gini + weight_df$Relief

# 按累计权重排序
weight_df <- weight_df %>% arrange(desc(Cumulative_Weight))

# 筛选累计权重 >= 3.0 的基因
selected_genes <- weight_df %>% filter(Cumulative_Weight >= 3.0) %>% pull(Gene)

cat("  筛选阈值: Cumulative_Weight >= 3.0\n")
cat("  筛选出核心基因数:", length(selected_genes), "\n")

# ============================================================================
# 5. 保存结果
# ============================================================================
cat("\n[5/5] 保存结果...\n")

# 保存完整权重表
write_csv(weight_df, file.path(out_dir, "attribute_weighting_all_201genes.csv"))
cat("  已保存: attribute_weighting_all_201genes.csv\n")

# 保存筛选后的33个核心基因列表
write_lines(selected_genes, file.path(out_dir, "core_33_genes.txt"))
cat("  已保存: core_33_genes.txt\n")

# 保存筛选后的核心基因详细信息（含权重）
core_genes_detail <- weight_df %>% filter(Gene %in% selected_genes)
write_csv(core_genes_detail, file.path(out_dir, "core_33_genes_with_weights.csv"))
cat("  已保存: core_33_genes_with_weights.csv\n")

# 保存筛选报告
report <- c(
  "========================================",
  "属性加权筛选报告",
  "========================================",
  paste("分析时间:", Sys.time()),
  "",
  "输入数据:",
  paste("  总关键基因数:", length(all_genes)),
  paste("  上调基因:", length(up_genes)),
  paste("  下调基因:", length(down_genes)),
  "",
  "属性加权方法:",
  "  1. 信息增益 (Information Gain)",
  "  2. 增益比率 (Gain Ratio)",
  "  3. 基尼指数 (Gini Index)",
  "  4. Relief算法",
  "",
  "筛选标准:",
  "  累计权重 (IG + GainRatio + Gini + Relief) >= 3.0",
  "",
  "筛选结果:",
  paste("  核心基因数:", length(selected_genes)),
  "",
  "核心基因列表:",
  paste("  ", 1:length(selected_genes), ".", selected_genes, sep = ""),
  "",
  "========================================"
)

write_lines(report, file.path(out_dir, "attribute_weighting_report.txt"))
cat("  已保存: attribute_weighting_report.txt\n")

# ============================================================================
# 6. 打印结果摘要
# ============================================================================
cat("\n========================================\n")
cat("筛选结果摘要\n")
cat("========================================\n")
cat("总关键基因:", length(all_genes), "\n")
cat("筛选后核心基因:", length(selected_genes), "\n")
cat("\nTop 10 核心基因:\n")
print(head(core_genes_detail, 10))

cat("\n✅ 完成！结果保存在:", out_dir, "\n")
cat("  - attribute_weighting_all_201genes.csv (全部201个基因的权重)\n")
cat("  - core_33_genes.txt (33个核心基因列表)\n")
cat("  - core_33_genes_with_weights.csv (33个核心基因+权重)\n")
cat("  - attribute_weighting_report.txt (分析报告)\n")
