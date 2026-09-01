#!/usr/bin/env Rscript
# ============================================================================
# 标准GO/KEGG富集分析（修正版 v3）
# 修复：KEGG和GO注释的列名问题
# ============================================================================

library(tidyverse)

# setwd 已移除，使用相对路径
cat("========================================\n")
cat("GO/KEGG 富集分析（修正版 v3）\n")
cat("========================================\n")

# ============================================================================
# 1. 读取数据
# ============================================================================
cat("\n[1/6] 读取基因列表...\n")

up <- read_lines("data/key_genes_up.txt")
down <- read_lines("data/key_genes_down.txt")
cat("  上调基因:", length(up), "\n")
cat("  下调基因:", length(down), "\n")

cat("\n[2/6] 读取背景基因...\n")
combat <- readRDS("data/combat_corrected_matrix.rds")
background <- rownames(combat)
cat("  背景基因数:", length(background), "\n")

cat("\n[3/6] 读取eggNOG注释...\n")
eggnog <- read_tsv(file.path("data", "out.emapper.annotations"),
                   comment = "##", skip = 4,
                   col_names = c("query","seed_ortholog","evalue","score",
                                 "eggNOG_OGs","max_annot_lvl","COG_category",
                                 "Description","Preferred_name","GOs","EC",
                                 "KEGG_ko","KEGG_Pathway","KEGG_Module",
                                 "KEGG_Reaction","KEGG_rclass","BRITE",
                                 "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"),
                   na = c("-","","NA"), show_col_types = FALSE)

# 清理基因ID
clean_id <- function(x) {
  x <- str_remove(x, "\\.[0-9]+\\.[0-9]+C\\.[0-9]+$")
  x <- str_remove(x, "\\.[0-9]+\\.[0-9]+C$")
  x <- str_remove(x, "\\.[0-9]+$")
  x
}
eggnog$gene <- clean_id(eggnog$query)

# 只保留背景中的基因
eggnog <- eggnog %>% filter(gene %in% background)
cat("  背景中有注释的基因数:", nrow(eggnog), "\n")

# ============================================================================
# 4. 读取GO/KEGG名称映射
# ============================================================================
cat("\n[4/6] 读取GO/KEGG名称映射...\n")

go_names <- read_tsv("data/GO_term2name.tsv",
                     col_names = c("term", "Name"), show_col_types = FALSE)
cat("  GO term名称:", nrow(go_names), "\n")

kegg_names <- read_tsv("data/KEGG_term2name.tsv",
                       col_names = c("term", "Name"), show_col_types = FALSE)
cat("  KEGG通路名称:", nrow(kegg_names), "\n")

# ============================================================================
# 5. 解析GO注释（按BP/CC/MF分类）
# ============================================================================
cat("\n[5/6] 解析GO注释...\n")

go_annot <- eggnog %>%
  filter(!is.na(GOs), GOs != "", GOs != "-") %>%
  select(gene, GOs) %>%
  separate_rows(GOs, sep = ",") %>%
  filter(grepl("^GO:", GOs)) %>%
  distinct() %>%
  rename(term = GOs)

# 根据GO号判断类别
go_annot <- go_annot %>%
  mutate(
    go_num = as.numeric(str_extract(term, "[0-9]+")),
    GO_type = case_when(
      go_num >= 4000000 ~ "BP",
      go_num >= 2000000 & go_num < 4000000 ~ "CC",
      go_num < 2000000 ~ "MF",
      TRUE ~ "unknown"
    )
  ) %>%
  filter(GO_type != "unknown") %>%
  select(-go_num)

cat("  GO注释完成\n")

# ============================================================================
# 6. 解析KEGG注释（修复：明确列名为term）
# ============================================================================
cat("  解析KEGG注释...\n")

kegg_annot <- eggnog %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
  select(gene, KEGG_Pathway) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  filter(KEGG_Pathway != "") %>%
  mutate(
    # 提取ko号
    term = str_extract(KEGG_Pathway, "^ko[0-9]+")
  ) %>%
  filter(!is.na(term), term != "") %>%
  select(gene, term) %>%  # 明确保留 gene 和 term 列
  distinct()

cat("  KEGG注释基因数:", n_distinct(kegg_annot$gene), "\n")

# ============================================================================
# 7. 富集函数（修复：确保列名正确）
# ============================================================================
run_enrich <- function(test_genes, annot_df, background_all, label) {
  if(length(test_genes) < 3) {
    cat("  ", label, ": 测试基因太少 (<3)\n")
    return(NULL)
  }
  
  # 确保 annot_df 有 gene 和 term 列
  if(!all(c("gene", "term") %in% colnames(annot_df))) {
    cat("  ", label, ": 注释数据缺少 gene 或 term 列\n")
    return(NULL)
  }
  
  test_annot <- annot_df %>% filter(gene %in% test_genes)
  if(nrow(test_annot) == 0) {
    cat("  ", label, ": 测试基因无注释\n")
    return(NULL)
  }
  
  # 每个term的背景计数
  bg_counts <- annot_df %>%
    group_by(term) %>%
    summarise(Bg = n_distinct(gene), .groups = "drop")
  
  # 测试基因每个term的计数
  test_counts <- test_annot %>%
    group_by(term) %>%
    summarise(Count = n_distinct(gene), .groups = "drop")
  
  # 合并并计算超几何检验
  res <- test_counts %>%
    left_join(bg_counts, by = "term") %>%
    mutate(
      TotalBg = length(background_all),
      TestGenes = length(test_genes),
      pvalue = phyper(Count - 1, Bg, TotalBg - Bg, TestGenes, lower.tail = FALSE)
    ) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    arrange(padj)
  
  return(res)
}

# ============================================================================
# 8. 执行富集分析
# ============================================================================
cat("\n[6/6] 执行富集分析...\n")

# KEGG
cat("  KEGG 上调...\n")
kegg_up <- run_enrich(up, kegg_annot, background, "KEGG_up")
cat("  KEGG 下调...\n")
kegg_down <- run_enrich(down, kegg_annot, background, "KEGG_down")

# GO BP/CC/MF 分别分析
cat("  GO_BP 上调...\n")
go_up_bp <- run_enrich(up, go_annot %>% filter(GO_type == "BP"), background, "GO_up_BP")
cat("  GO_BP 下调...\n")
go_down_bp <- run_enrich(down, go_annot %>% filter(GO_type == "BP"), background, "GO_down_BP")

cat("  GO_CC 上调...\n")
go_up_cc <- run_enrich(up, go_annot %>% filter(GO_type == "CC"), background, "GO_up_CC")
cat("  GO_CC 下调...\n")
go_down_cc <- run_enrich(down, go_annot %>% filter(GO_type == "CC"), background, "GO_down_CC")

cat("  GO_MF 上调...\n")
go_up_mf <- run_enrich(up, go_annot %>% filter(GO_type == "MF"), background, "GO_up_MF")
cat("  GO_MF 下调...\n")
go_down_mf <- run_enrich(down, go_annot %>% filter(GO_type == "MF"), background, "GO_down_MF")

# ============================================================================
# 9. 添加名称并保存
# ============================================================================
out_dir <- "results/enrichment_standard_fixed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 添加KEGG名称
add_kegg_names <- function(df) {
  if(is.null(df) || nrow(df) == 0) return(df)
  df %>% left_join(kegg_names, by = "term")
}

# 添加GO名称
add_go_names <- function(df) {
  if(is.null(df) || nrow(df) == 0) return(df)
  df %>% left_join(go_names, by = "term")
}

# 保存结果
save_results <- function(df, prefix, suffix, label) {
  if(is.null(df) || nrow(df) == 0) {
    cat("  ", label, ": 无显著结果 (padj < 0.05)\n")
    return()
  }
  # 过滤显著结果
  df_sig <- df %>% filter(padj < 0.05)
  if(nrow(df_sig) == 0) {
    cat("  ", label, ": 无显著结果 (padj < 0.05)\n")
    return()
  }
  file_name <- file.path(out_dir, paste0(prefix, "_", suffix, ".csv"))
  write_csv(df_sig, file_name)
  cat("  ", label, ": 保存", nrow(df_sig), "条显著通路到", basename(file_name), "\n")
}

cat("\n保存结果...\n")

# KEGG
save_results(add_kegg_names(kegg_up), "KEGG", "up", "KEGG上调")
save_results(add_kegg_names(kegg_down), "KEGG", "down", "KEGG下调")

# GO BP
save_results(add_go_names(go_up_bp), "GO_BP", "up", "GO_BP上调")
save_results(add_go_names(go_down_bp), "GO_BP", "down", "GO_BP下调")

# GO CC
save_results(add_go_names(go_up_cc), "GO_CC", "up", "GO_CC上调")
save_results(add_go_names(go_down_cc), "GO_CC", "down", "GO_CC下调")

# GO MF
save_results(add_go_names(go_up_mf), "GO_MF", "up", "GO_MF上调")
save_results(add_go_names(go_down_mf), "GO_MF", "down", "GO_MF下调")

# ============================================================================
# 10. 生成报告
# ============================================================================
report <- c(
  "========================================",
  "GO/KEGG 富集分析报告（修正版 v3）",
  "========================================",
  paste("分析时间:", Sys.time()),
  "",
  "数据统计:",
  paste("  背景基因数:", length(background)),
  paste("  上调基因数:", length(up)),
  paste("  下调基因数:", length(down)),
  "",
  "KEGG显著通路 (padj < 0.05):",
  paste("  上调:", ifelse(is.null(kegg_up), 0, nrow(kegg_up %>% filter(padj < 0.05)))),
  paste("  下调:", ifelse(is.null(kegg_down), 0, nrow(kegg_down %>% filter(padj < 0.05)))),
  "",
  "GO_BP显著通路 (padj < 0.05):",
  paste("  上调:", ifelse(is.null(go_up_bp), 0, nrow(go_up_bp %>% filter(padj < 0.05)))),
  paste("  下调:", ifelse(is.null(go_down_bp), 0, nrow(go_down_bp %>% filter(padj < 0.05)))),
  "",
  "GO_CC显著通路 (padj < 0.05):",
  paste("  上调:", ifelse(is.null(go_up_cc), 0, nrow(go_up_cc %>% filter(padj < 0.05)))),
  paste("  下调:", ifelse(is.null(go_down_cc), 0, nrow(go_down_cc %>% filter(padj < 0.05)))),
  "",
  "GO_MF显著通路 (padj < 0.05):",
  paste("  上调:", ifelse(is.null(go_up_mf), 0, nrow(go_up_mf %>% filter(padj < 0.05)))),
  paste("  下调:", ifelse(is.null(go_down_mf), 0, nrow(go_down_mf %>% filter(padj < 0.05)))),
  "",
  "========================================"
)

write_lines(report, file.path(out_dir, "analysis_report.txt"))
cat("\n", paste(report, collapse = "\n"), "\n")

cat("\n✅ 完成！结果保存在:", out_dir, "\n")
