#!/usr/bin/env Rscript
# ============================================================================
# 联合分析：limma + WGCNA（直接从 combat_corrected_matrix.rds 读取）
# 参数：2+2模块，MM>0.85，GS>0.3
# 输入：results/PCA/tables/combat_corrected_matrix.rds, data/metadata.csv
# 输出：results/joint_analysis_from_matrix/
# ============================================================================

set.seed(12345)

library(limma)
library(WGCNA)
library(tidyverse)

enableWGCNAThreads()

cat("========================================\n")
cat("联合分析：limma + WGCNA（从矩阵读取）\n")
cat("========================================\n")

# ============================================================================
# 1. 定义路径
# ============================================================================
project_root <- getwd()
out_dir <- file.path(getwd(), "results", "joint_analysis_from_matrix")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. 读取数据（直接读取已校正的矩阵）
# ============================================================================
cat("\n[1/5] 读取数据...\n")

expr <- readRDS("data/combat_corrected_matrix.rds")
sample_info <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)

# 确保 group 是因子
y <- factor(sample_info$group, levels = c("control", "treatment"))

cat("基因数:", nrow(expr), "\n")
cat("样本数:", ncol(expr), "\n")
cat("对照组:", sum(y == "control"), "个\n")
cat("处理组:", sum(y == "treatment"), "个\n")

# ============================================================================
# 3. limma 差异分析
# ============================================================================
cat("\n[2/5] limma 差异分析...\n")

design <- model.matrix(~ y)
colnames(design) <- c("Intercept", "treatment")

fit <- lmFit(expr, design)
fit <- eBayes(fit)
limma_results <- topTable(fit, coef = "treatment", number = Inf, adjust.method = "BH")
limma_results$Gene <- rownames(limma_results)

# 筛选差异基因
deg_up <- limma_results %>% filter(logFC > 1, adj.P.Val < 1e-6) %>% pull(Gene)
deg_down <- limma_results %>% filter(logFC < -1, adj.P.Val < 1e-6) %>% pull(Gene)

cat(sprintf("上调DEGs: %d\n", length(deg_up)))
cat(sprintf("下调DEGs: %d\n", length(deg_down)))

write_csv(limma_results, file.path(tab_dir, "limma_all_results.csv"))

# ============================================================================
# 4. WGCNA 分析
# ============================================================================
cat("\n[3/5] WGCNA 分析（前25%% 高变异基因）...\n")

datExpr <- t(expr)

# 选择前25%高变异基因
gene_vars <- apply(datExpr, 2, var, na.rm = TRUE)
n_top <- round(length(gene_vars) * 0.25)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_top]
datExpr_filtered <- datExpr[, top_genes]
cat(sprintf("WGCNA 输入基因数: %d\n", ncol(datExpr_filtered)))

# 软阈值选择
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr_filtered, powerVector = powers, verbose = 0)

best_power <- sft$powerEstimate
if (is.na(best_power)) best_power <- 12
cat("选用的软阈值: β =", best_power, "\n")

# 构建网络
adjacency <- adjacency(datExpr_filtered, power = best_power)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 动态剪枝
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

# 合并相似模块
merge <- mergeCloseModules(datExpr_filtered, dynamicColors,
                           cutHeight = 0.25, verbose = 0)
moduleColors <- merge$colors
MEs <- merge$newMEs
moduleColors <- as.character(moduleColors)
names(moduleColors) <- colnames(datExpr_filtered)

cat(sprintf("模块数量: %d\n", length(unique(moduleColors))))

# ============================================================================
# 5. 模块-性状关联
# ============================================================================
cat("\n[4/5] 模块-性状关联...\n")

trait <- data.frame(Stress = ifelse(y == "treatment", 1, 0))
rownames(trait) <- rownames(datExpr_filtered)

geneModuleMembership <- as.data.frame(cor(datExpr_filtered, MEs, use = "p"))
geneTraitSignificance <- as.data.frame(cor(datExpr_filtered, trait$Stress, use = "p"))
colnames(geneTraitSignificance) <- "GS"

moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))

# 按相关性排序
mod_df <- data.frame(
  Module = gsub("ME", "", rownames(moduleTraitCor)),
  Correlation = moduleTraitCor[,1],
  Pvalue = moduleTraitPvalue[,1]
) %>% arrange(desc(Correlation))

write_csv(mod_df, file.path(tab_dir, "module_trait_correlations.csv"))

cat("\n模块相关性排序:\n")
for (i in 1:min(10, nrow(mod_df))) {
  cat(sprintf("  %s: r = %.4f (p = %.4f)\n",
              mod_df$Module[i], mod_df$Correlation[i], mod_df$Pvalue[i]))
}

# ============================================================================
# 6. 选择模块（正相关最强2个 + 负相关最强2个）
# ============================================================================
cat("\n[5/5] 选择模块并筛选Hub基因...\n")

up_modules <- mod_df %>% filter(Correlation > 0) %>% head(2)
down_modules <- mod_df %>% filter(Correlation < 0) %>% tail(2)

cat(sprintf("正相关最强2个模块: %s\n", paste(up_modules$Module, collapse=", ")))
cat(sprintf("负相关最强2个模块: %s\n", paste(down_modules$Module, collapse=", ")))

# ============================================================================
# 7. Hub 基因筛选（MM > 0.85, GS > 0.3）
# ============================================================================
hub_up <- c()
hub_down <- c()

for (mod_name in up_modules$Module) {
  mod_genes <- names(moduleColors)[moduleColors == mod_name]
  me_name <- paste0("ME", mod_name)
  mm <- abs(geneModuleMembership[mod_genes, me_name])
  gs <- abs(geneTraitSignificance[mod_genes, "GS"])
  hub <- mod_genes[mm > 0.85 & gs > 0.3]
  hub_up <- c(hub_up, hub)
  cat(sprintf("  上调模块 %s: %d 个hub (MM>0.85, GS>0.3)\n", mod_name, length(hub)))
}

for (mod_name in down_modules$Module) {
  mod_genes <- names(moduleColors)[moduleColors == mod_name]
  me_name <- paste0("ME", mod_name)
  mm <- abs(geneModuleMembership[mod_genes, me_name])
  gs <- abs(geneTraitSignificance[mod_genes, "GS"])
  hub <- mod_genes[mm > 0.85 & gs > 0.3]
  hub_down <- c(hub_down, hub)
  cat(sprintf("  下调模块 %s: %d 个hub (MM>0.85, GS>0.3)\n", mod_name, length(hub)))
}

cat(sprintf("\n总上调hub: %d\n", length(hub_up)))
cat(sprintf("总下调hub: %d\n", length(hub_down)))

# ============================================================================
# 8. 与limma取交集
# ============================================================================
up_intersect <- intersect(deg_up, hub_up)
down_intersect <- intersect(deg_down, hub_down)

cat("\n========================================\n")
cat("联合分析结果\n")
cat("========================================\n")
cat(sprintf("上调交集: %d\n", length(up_intersect)))
cat(sprintf("下调交集: %d\n", length(down_intersect)))
cat(sprintf("总计: %d\n", length(up_intersect) + length(down_intersect)))

# 保存结果
write_lines(up_intersect, file.path(tab_dir, "key_genes_up.txt"))
write_lines(down_intersect, file.path(tab_dir, "key_genes_down.txt"))
write_lines(c(up_intersect, down_intersect), file.path(tab_dir, "key_genes_all.txt"))

# 保存详细信息
if (length(up_intersect) > 0) {
  up_detail <- data.frame(
    Gene = up_intersect,
    logFC = limma_results$logFC[match(up_intersect, limma_results$Gene)],
    FDR = limma_results$adj.P.Val[match(up_intersect, limma_results$Gene)]
  )
  write_csv(up_detail, file.path(tab_dir, "key_genes_up_detail.csv"))
}

if (length(down_intersect) > 0) {
  down_detail <- data.frame(
    Gene = down_intersect,
    logFC = limma_results$logFC[match(down_intersect, limma_results$Gene)],
    FDR = limma_results$adj.P.Val[match(down_intersect, limma_results$Gene)]
  )
  write_csv(down_detail, file.path(tab_dir, "key_genes_down_detail.csv"))
}

# 保存hub基因完整列表
hub_all <- data.frame(
  Gene = c(hub_up, hub_down),
  Direction = c(rep("Up", length(hub_up)), rep("Down", length(hub_down))),
  Module = c(
    rep(up_modules$Module[1], length(hub_up[hub_up %in% names(moduleColors)[moduleColors == up_modules$Module[1]]])),
    # 简化版本，只记录方向
  )
)
write_csv(data.frame(
  Gene = c(hub_up, hub_down),
  Direction = c(rep("Up", length(hub_up)), rep("Down", length(hub_down)))
), file.path(tab_dir, "wgcna_hub_genes.csv"))

cat("\n✅ 结果保存在:", out_dir, "\n")
