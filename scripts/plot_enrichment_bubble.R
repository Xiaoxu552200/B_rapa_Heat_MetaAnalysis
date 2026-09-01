#!/usr/bin/env Rscript
# ============================================================================
# 富集分析气泡图绘制脚本
# 功能：绘制上调/下调基因的GO和KEGG富集气泡图
# 输入：enrichment_correct/GO_up.csv, GO_down.csv, KEGG_up_final.csv, KEGG_down_final.csv
# 输出：enrichment_correct/figures/ 下的TIFF图片
# ============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

cat("========================================\n")
cat("富集分析气泡图绘制\n")
cat("========================================\n")

# ============================================================================
# 1. 定义路径
# ============================================================================
project_root <- getwd()
out_dir <- file.path(project_root, "results", "enrichment_correct")
fig_dir <- file.path(out_dir, "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. 读取数据
# ============================================================================
cat("\n[1/3] 读取富集结果...\n")

# 读取GO结果
up_go <- read_csv(file.path(out_dir, "GO_up.csv"), show_col_types = FALSE)
down_go <- read_csv(file.path(out_dir, "GO_down.csv"), show_col_types = FALSE)

# 读取KEGG结果
up_kegg <- read_csv(file.path(out_dir, "KEGG_up_final.csv"), show_col_types = FALSE)
down_kegg <- read_csv(file.path(out_dir, "KEGG_down_final.csv"), show_col_types = FALSE)

cat("上调GO:", nrow(up_go), "条\n")
cat("下调GO:", nrow(down_go), "条\n")
cat("上调KEGG:", nrow(up_kegg), "条\n")
cat("下调KEGG:", nrow(down_kegg), "条\n")

# ============================================================================
# 3. 绘图函数
# ============================================================================
plot_bubble <- function(data, title, max_terms = 20) {
  if (is.null(data) || nrow(data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No significant terms") + 
             theme_void())
  }
  
  plot_data <- data %>%
    arrange(padj) %>%
    head(max_terms) %>%
    mutate(
      log10padj = -log10(padj),
      GeneRatio = Count / TargetSize,
      Name = factor(Name, levels = rev(unique(Name)))
    )
  
  ggplot(plot_data, aes(x = GeneRatio, y = Name)) +
    geom_point(aes(size = Count, color = log10padj), stroke = 0) +
    scale_color_gradient(low = "blue", high = "red", 
                         name = expression(-log[10](padj))) +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
}

# ============================================================================
# 4. 生成组合图
# ============================================================================
cat("\n[2/3] 绘制组合图...\n")

p_up_go <- plot_bubble(up_go, "Up-regulated GO BP")
p_up_kegg <- plot_bubble(up_kegg, "Up-regulated KEGG Pathways")
p_down_go <- plot_bubble(down_go, "Down-regulated GO BP")
p_down_kegg <- plot_bubble(down_kegg, "Down-regulated KEGG Pathways")

# 组合图：GO + KEGG 并排
combined_up <- p_up_go + p_up_kegg + plot_layout(widths = c(2, 1))
combined_down <- p_down_go + p_down_kegg + plot_layout(widths = c(2, 1))

# ============================================================================
# 5. 保存图片
# ============================================================================
cat("\n[3/3] 保存图片...\n")

height_up <- max(8, nrow(up_go) * 0.35 + 2)
tiff(file.path(fig_dir, "Figure_Up_GO_KEGG_bubble.tiff"), 
     width = 14, height = height_up, units = "in", res = 300, compression = "lzw")
print(combined_up)
dev.off()
cat("已保存: Figure_Up_GO_KEGG_bubble.tiff\n")

height_down <- max(8, nrow(down_go) * 0.35 + 2)
tiff(file.path(fig_dir, "Figure_Down_GO_KEGG_bubble.tiff"), 
     width = 14, height = height_down, units = "in", res = 300, compression = "lzw")
print(combined_down)
dev.off()
cat("已保存: Figure_Down_GO_KEGG_bubble.tiff\n")

# ============================================================================
# 6. 也保存PDF版本
# ============================================================================
pdf(file.path(fig_dir, "Figure_Up_GO_KEGG_bubble.pdf"), width = 14, height = height_up)
print(combined_up)
dev.off()
cat("已保存: Figure_Up_GO_KEGG_bubble.pdf\n")

pdf(file.path(fig_dir, "Figure_Down_GO_KEGG_bubble.pdf"), width = 14, height = height_down)
print(combined_down)
dev.off()
cat("已保存: Figure_Down_GO_KEGG_bubble.pdf\n")

cat("\n========================================\n")
cat("✅ 全部完成！\n")
cat("输出目录:", fig_dir, "\n")
cat("========================================\n")
