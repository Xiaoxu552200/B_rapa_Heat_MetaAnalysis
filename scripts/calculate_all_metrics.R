#!/usr/bin/env Rscript
# ============================================================================
# 计算三个模型的所有评估指标
# 基于5折交叉验证结果
# XGBoost, GBM, SVM
# ============================================================================

library(tidyverse)
library(caret)
library(pROC)
library(xgboost)
library(e1071)
library(gbm)
library(CORElearn)

cat("====================================================\n")
cat("三个模型 - 完整评估指标计算\n")
cat("====================================================\n")

set.seed(42)

# ============================================================================
# 1. 设置路径
# ============================================================================
base_dir <- getwd()

out_dir <- file.path(base_dir, "results", "ML")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n[1/6] 输出目录:", out_dir, "\n")

# ============================================================================
# 2. 加载33个核心基因
# ============================================================================
cat("\n[2/6] 加载33个核心基因...\n")

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

X <- t(expr[all_genes, ])
X_scaled <- scale(X)

X_df <- as.data.frame(X_scaled)
X_df$target <- as.factor(y)

# 属性加权筛选33个基因
ig_weights <- attrEval(target ~ ., data = X_df, estimator = "InfGain")
gr_weights <- attrEval(target ~ ., data = X_df, estimator = "GainRatio")
gini_weights <- attrEval(target ~ ., data = X_df, estimator = "Gini")
relief_weights <- attrEval(target ~ ., data = X_df, estimator = "ReliefFequalK")

normalize_01 <- function(x) {
  if (max(x) - min(x) == 0) return(rep(0, length(x)))
  return((x - min(x)) / (max(x) - min(x)))
}

weight_df <- data.frame(
  Gene = all_genes,
  IG = normalize_01(ig_weights[all_genes]),
  GR = normalize_01(gr_weights[all_genes]),
  Gini = normalize_01(gini_weights[all_genes]),
  Relief = normalize_01(relief_weights[all_genes])
)

weight_df$Cumulative_Weight <- weight_df$IG + weight_df$GR + weight_df$Gini + weight_df$Relief
weight_df <- weight_df %>% arrange(desc(Cumulative_Weight))

selected_genes <- weight_df %>% filter(Cumulative_Weight >= 3.0) %>% pull(Gene)
cat(sprintf("  筛选出 %d 个核心基因\n", length(selected_genes)))

# 保存33个基因列表到ML目录
write_lines(selected_genes, file.path(out_dir, "core_33_genes.txt"))

X_sel <- X_scaled[, selected_genes]

# ============================================================================
# 3. 三个模型参数
# ============================================================================
cat("\n[3/6] 设置模型参数...\n")

# XGBoost
xgb_params <- list(
  objective = "binary:logistic",
  max_depth = 3,
  eta = 0.03,
  subsample = 0.5,
  colsample_bytree = 0.5,
  min_child_weight = 2,
  lambda = 0.5,
  alpha = 0.5
)
xgb_nrounds <- 15

# GBM
gbm_n_trees <- 20
gbm_depth <- 2
gbm_shrinkage <- 0.005

# SVM
svm_degree <- 3
svm_cost <- 0.03
svm_gamma <- 0.001

# ============================================================================
# 4. 5折交叉验证预测函数
# ============================================================================
get_cv_predictions <- function(model_type) {
  set.seed(42)
  folds <- createFolds(y, k = 5, list = TRUE)
  all_preds <- c()
  all_true <- c()
  
  for (fold in 1:5) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(1:length(y), test_idx)
    
    X_train <- X_sel[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X_sel[test_idx, ]
    y_test <- y[test_idx]
    
    if (model_type == "xgb") {
      dtrain <- xgb.DMatrix(X_train, label = y_train)
      dtest <- xgb.DMatrix(X_test)
      model <- xgb.train(params = xgb_params, data = dtrain, nrounds = xgb_nrounds, verbose = 0)
      pred <- predict(model, dtest)
      
    } else if (model_type == "gbm") {
      df <- as.data.frame(X_train)
      df$y <- y_train
      model <- gbm(y ~ ., data = df, distribution = "bernoulli",
                   n.trees = gbm_n_trees,
                   interaction.depth = gbm_depth,
                   shrinkage = gbm_shrinkage,
                   verbose = FALSE)
      pred <- predict(model, as.data.frame(X_test), n.trees = gbm_n_trees, type = "response")
      
    } else if (model_type == "svm") {
      svm_model <- svm(X_train, as.factor(y_train), 
                       kernel = "polynomial",
                       degree = svm_degree,
                       cost = svm_cost,
                       gamma = svm_gamma,
                       probability = TRUE)
      pred <- attr(predict(svm_model, X_test, probability = TRUE), "probabilities")[, "1"]
    }
    
    all_preds <- c(all_preds, pred)
    all_true <- c(all_true, y_test)
  }
  
  return(list(preds = all_preds, true = all_true))
}

cat("\n[4/6] 5折交叉验证...\n")

model_names <- c("XGBoost", "GBM", "SVM")
model_types <- c("xgb", "gbm", "svm")

results <- list()

for (i in 1:length(model_types)) {
  cat(sprintf("  运行: %s...\n", model_names[i]))
  res <- get_cv_predictions(model_types[i])
  results[[i]] <- res
}
names(results) <- model_types

# ============================================================================
# 5. 计算评估指标
# ============================================================================
cat("\n[5/6] 计算评估指标...\n")

calculate_metrics <- function(true_labels, pred_probs, threshold = 0.5) {
  pred_class <- ifelse(pred_probs >= threshold, 1, 0)
  
  cm <- table(Predicted = pred_class, Actual = true_labels)
  
  TP <- ifelse(is.na(cm[2, 2]), 0, cm[2, 2])
  FN <- ifelse(is.na(cm[1, 2]), 0, cm[1, 2])
  FP <- ifelse(is.na(cm[2, 1]), 0, cm[2, 1])
  TN <- ifelse(is.na(cm[1, 1]), 0, cm[1, 1])
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  npv <- ifelse((TN + FN) > 0, TN / (TN + FN), 0)
  f1 <- ifelse((precision + sensitivity) > 0, 
               2 * precision * sensitivity / (precision + sensitivity), 0)
  
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE)
  auc_val <- roc_obj$auc
  
  return(list(
    AUC = auc_val,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    NPV = npv,
    F1 = f1,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  ))
}

all_metrics <- data.frame()

for (i in 1:length(model_types)) {
  res <- results[[i]]
  metrics <- calculate_metrics(res$true, res$preds)
  
  all_metrics <- rbind(all_metrics, data.frame(
    Model = model_names[i],
    AUC = round(metrics$AUC, 4),
    Accuracy = round(metrics$Accuracy, 4),
    Sensitivity = round(metrics$Sensitivity, 4),
    Specificity = round(metrics$Specificity, 4),
    Precision = round(metrics$Precision, 4),
    NPV = round(metrics$NPV, 4),
    F1 = round(metrics$F1, 4)
  ))
}

# ============================================================================
# 6. 保存结果
# ============================================================================
cat("\n[6/6] 保存结果...\n")

cat("\n====================================================\n")
cat("模型评估指标汇总\n")
cat("====================================================\n\n")
print(all_metrics, row.names = FALSE)

# 保存CSV
write.csv(all_metrics, file.path(out_dir, "Model_Performance_Metrics.csv"), row.names = FALSE)
cat(sprintf("\n✅ CSV已保存: %s\n", file.path(out_dir, "Model_Performance_Metrics.csv")))

# 保存详细的预测结果
for (i in 1:length(model_types)) {
  res <- results[[i]]
  pred_df <- data.frame(
    Sample = paste0("Sample_", 1:length(res$true)),
    True = res$true,
    Pred_Prob = res$preds,
    Pred_Class = ifelse(res$preds >= 0.5, 1, 0)
  )
  write.csv(pred_df, 
            file.path(out_dir, paste0("Predictions_", model_names[i], ".csv")), 
            row.names = FALSE)
  cat(sprintf("  已保存: Predictions_%s.csv\n", model_names[i]))
}

cat("\n✅ 全部完成！\n")
cat("结果保存在:", out_dir, "\n")

