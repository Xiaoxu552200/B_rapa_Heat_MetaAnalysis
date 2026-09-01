# 白菜高温胁迫转录组学分析与关键基因挖掘

## 项目简介
本研究整合多个公共RNA-seq数据集，通过limma差异表达分析、WGCNA共表达网络构建、属性加权特征筛选和机器学习分类建模（XGBoost/GBM/SVM），系统鉴定了白菜热胁迫响应的核心基因。

## 数据来源
- 原始RNA-seq数据：NCBI SRA (PRJNA872510, PRJNA1030162, PRJNA646007, PRJNA663233)
- 参考基因组：Brassica rapa Chiifu v3.5 (BRAD数据库)

## 运行环境
- R >= 4.0
- Python >= 3.8
- 必需R包：tidyverse, limma, WGCNA, ggplot2, pheatmap, caret, xgboost, e1071, gbm, CORElearn, pROC, here
- 必需Python包：pandas, numpy, matplotlib, seaborn, shap, xgboost, scikit-learn

## 快速开始
1. 双击 `B_rapa_Heat_MetaAnalysis.Rproj` 打开RStudio
2. 按顺序运行 `scripts/` 目录下的脚本

## 目录结构
- `scripts/`: 所有分析脚本
- `data/`: 输入数据（需自行下载）
- `results/`: 输出结果

## 联系方式
徐彦泙 - 22316191@zju.edu.cn
### 数据获取

大文件已上传至 Zenodo：
> https://zenodo.org/record/8411138](https://zenodo.org/records/22236209
