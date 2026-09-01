#!/usr/bin/env python3
# ============================================================================
# SHAP特征重要性分析 - 33个核心基因
# 脚本路径: ./scripts
# 输出路径: ./results/ML
# 输出格式: TIFF
# ============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier
import shap
import warnings
import os
warnings.filterwarnings('ignore')

# ============================================================================
# 路径设置
# ============================================================================

# 脚本所在目录 (Windows格式用于显示，Linux格式用于操作)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "results", "ML")
DATA_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "data")

# 创建输出目录
os.makedirs(f"{OUTPUT_DIR}/tables", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)

print("=" * 60)
print("SHAP特征重要性分析 - 33个核心基因")
print("=" * 60)
print(f"\n脚本路径: {SCRIPT_DIR}")
print(f"输出路径: {OUTPUT_DIR}")
print(f"数据路径: {DATA_DIR}")

# ============================================================================
# 1. 加载数据
# ============================================================================

print("\n【1】加载数据...")

# 切换到数据目录加载文件
os.chdir(DATA_DIR)

expr = pd.read_csv("results/PCA/tables/combat_corrected_matrix.csv", index_col=0)
print(f"  表达矩阵: {expr.shape[0]} 基因 × {expr.shape[1]} 样本")

with open("joint_analysis_201/tables/core_33_genes.txt", "r") as f:
    selected_genes = [line.strip() for line in f.readlines() if line.strip()]

print(f"  核心基因: {len(selected_genes)} 个")

# 匹配基因
expr_genes = expr.index.tolist()
matched_genes = []

for g in selected_genes:
    if g in expr_genes:
        matched_genes.append(g)
    else:
        g_clean = g.replace('.D', '').replace('_D', '').replace('.d', '').replace('_d', '')
        for eg in expr_genes:
            eg_clean = eg.replace('.3.5C', '').replace('.D', '').replace('_D', '')
            if g_clean == eg_clean or g_clean in eg:
                matched_genes.append(eg)
                break

matched_genes = list(dict.fromkeys(matched_genes))[:33]
print(f"  成功匹配: {len(matched_genes)} 个基因")

# 提取数据
X_raw = expr.loc[matched_genes, :].T
metadata = pd.read_csv("data/metadata.csv")
y = (metadata['group'] == 'treatment').astype(int).values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)

print(f"  样本数: {X_scaled.shape[0]}, 特征数: {X_scaled.shape[1]}")

# ============================================================================
# 2. 训练XGBoost模型
# ============================================================================

print("\n【2】训练XGBoost模型...")

params = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'max_depth': 3,
    'learning_rate': 0.1,
    'n_estimators': 15,
    'subsample': 0.8,
    'colsample_bytree': 0.4,
    'min_child_weight': 2,
    'reg_lambda': 1,
    'reg_alpha': 0.5,
    'random_state': 42,
    'use_label_encoder': False
}

# 5折交叉验证
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
all_preds = []
all_true = []

for fold, (train_idx, test_idx) in enumerate(skf.split(X_scaled, y)):
    X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    model = XGBClassifier(**params)
    model.fit(X_train, y_train, verbose=False)
    preds = model.predict_proba(X_test)[:, 1]
    
    all_preds.extend(preds)
    all_true.extend(y_test)

cv_auc = roc_auc_score(all_true, all_preds)
print(f"  5折交叉验证 AUC: {cv_auc:.4f}")

# 全量训练
model = XGBClassifier(**params)
model.fit(X_scaled, y, verbose=False)

# ============================================================================
# 3. SHAP分析
# ============================================================================

print("\n【3】计算SHAP值...")

try:
    booster = model.get_booster()
    explainer = shap.TreeExplainer(booster)
    shap_values = explainer.shap_values(X_scaled)
    print("  使用 TreeExplainer")
except Exception as e:
    print(f"  TreeExplainer失败: {e}")
    print("  切换到 KernelExplainer...")
    def predict_fn(x):
        return model.predict_proba(x)[:, 1]
    background = X_scaled[np.random.choice(X_scaled.shape[0], 20, replace=False)]
    explainer = shap.KernelExplainer(predict_fn, background)
    shap_values = explainer.shap_values(X_scaled, nsamples=100)
    print("  使用 KernelExplainer")

# ============================================================================
# 4. 基因名处理
# ============================================================================

gene_names_clean = []
for g in matched_genes:
    g_clean = g.replace('.3.5C', '').replace('.D', '').replace('_D', '').replace('.d', '').replace('_d', '')
    if '.' in g_clean:
        g_clean = g_clean.split('.')[0]
    gene_names_clean.append(g_clean)

# ============================================================================
# 5. SHAP特征重要性
# ============================================================================

print("\n【4】SHAP特征重要性排名...")

mean_shap = np.abs(shap_values).mean(axis=0)

shap_df = pd.DataFrame({
    'Gene': gene_names_clean,
    'Gene_Original': matched_genes,
    'Mean_ABS_SHAP': mean_shap
}).sort_values('Mean_ABS_SHAP', ascending=False)

shap_df['Rank'] = range(1, len(shap_df) + 1)

# Top 10
top10 = shap_df.head(10)
top10_genes = top10['Gene'].tolist()
top10_indices = [matched_genes.index(g) for g in top10['Gene_Original'].tolist()]

print("\n【Top 10 SHAP特征重要性】")
print("=" * 60)
print(top10[['Rank', 'Gene', 'Mean_ABS_SHAP']].to_string(index=False))

# ============================================================================
# 6. 保存表格
# ============================================================================

top10.to_csv(f"{OUTPUT_DIR}/tables/SHAP_Top10_genes.csv", index=False)
shap_df.to_csv(f"{OUTPUT_DIR}/tables/SHAP_All_33genes.csv", index=False)
print("\n  表格已保存")

# LaTeX表格
latex_lines = [
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Top 10 genes by SHAP importance}",
    "\\label{tab:shap_top10}",
    "\\begin{tabular}{crr}",
    "\\toprule",
    "Rank & Gene & Mean SHAP \\\\",
    "\\midrule"
]
for _, row in top10.iterrows():
    latex_lines.append(f"    {int(row['Rank'])} & {row['Gene']} & {row['Mean_ABS_SHAP']:.4f} \\\\")
latex_lines.extend([
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
])

with open(f"{OUTPUT_DIR}/tables/SHAP_Top10.tex", "w") as f:
    f.write("\n".join(latex_lines))
print("  LaTeX表格已保存: SHAP_Top10.tex")

# ============================================================================
# 7. 图1: 条形图 (排名图) - TIFF格式
# ============================================================================

print("\n【5】绘制SHAP图 (TIFF格式)...")

plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14

fig, ax = plt.subplots(figsize=(10, 8))
top10_sorted = top10.sort_values('Mean_ABS_SHAP', ascending=True)
bars = ax.barh(top10_sorted['Gene'], top10_sorted['Mean_ABS_SHAP'], 
               color='steelblue', edgecolor='navy', linewidth=0.5)
ax.set_xlabel('Mean |SHAP Value|', fontsize=12)
ax.set_ylabel('Gene', fontsize=12)
ax.set_title('SHAP Feature Importance - Top 10 Genes', fontsize=14)
ax.grid(axis='x', alpha=0.3)

for bar, val in zip(bars, top10_sorted['Mean_ABS_SHAP']):
    ax.text(val + 0.002, bar.get_y() + bar.get_height()/2, 
            f'{val:.4f}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/SHAP_Bar_Top10.tiff", dpi=300, bbox_inches='tight', format='tiff')
plt.close()
print(f"  排名图: {OUTPUT_DIR}/figures/SHAP_Bar_Top10.tiff")

# ============================================================================
# 8. 图2: 气泡图 (替代蜂群图) - TIFF格式
# ============================================================================

print("  生成气泡图...")

top10_shap = shap_values[:, top10_indices]
top10_X = X_scaled[:, top10_indices]

bubble_data = []
for i, gene in enumerate(top10_genes):
    for j in range(len(top10_shap)):
        bubble_data.append({
            'Gene': gene,
            'SHAP_Value': top10_shap[j, i],
            'Feature_Value': top10_X[j, i],
            'Sample': j
        })

bubble_df = pd.DataFrame(bubble_data)
gene_order = top10_genes[::-1]

fig, ax = plt.subplots(figsize=(12, 8))

for gene in gene_order:
    subset = bubble_df[bubble_df['Gene'] == gene]
    sizes = np.abs(subset['Feature_Value'].values) * 50 + 20
    colors_scatter = plt.cm.RdYlBu_r((subset['SHAP_Value'].values + 1) / 2)
    ax.scatter(subset['SHAP_Value'], [gene] * len(subset), 
               s=sizes, c=colors_scatter, alpha=0.7, edgecolors='gray', linewidth=0.3)

ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel('SHAP Value', fontsize=12)
ax.set_ylabel('Gene', fontsize=12)
ax.set_title('SHAP Bubble Plot - Top 10 Genes', fontsize=14)
ax.grid(axis='x', alpha=0.3)

sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlBu_r, norm=plt.Normalize(vmin=-1, vmax=1))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
cbar.set_label('SHAP Value (Red=Positive, Blue=Negative)', fontsize=10)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/SHAP_Bubble_Top10.tiff", dpi=300, bbox_inches='tight', format='tiff')
plt.close()
print(f"  气泡图: {OUTPUT_DIR}/figures/SHAP_Bubble_Top10.tiff")

# ============================================================================
# 9. 图3: 瀑布图 (Positive + Negative) - TIFF格式
# ============================================================================

print("  生成瀑布图...")

pred_probs = model.predict_proba(X_scaled)[:, 1]

pos_indices = np.where(y == 1)[0]
if len(pos_indices) > 0:
    pos_idx = pos_indices[np.argmax(pred_probs[y == 1])]
else:
    pos_idx = 0

neg_indices = np.where(y == 0)[0]
if len(neg_indices) > 0:
    neg_idx = neg_indices[np.argmin(pred_probs[y == 0])]
else:
    neg_idx = 1

for idx, label in [(pos_idx, 'Positive'), (neg_idx, 'Negative')]:
    
    sample_shap = shap_values[idx][top10_indices]
    sample_data = X_scaled[idx][top10_indices]
    
    exp = shap.Explanation(
        values=sample_shap,
        base_values=explainer.expected_value,
        data=sample_data,
        feature_names=top10_genes
    )
    
    plt.figure(figsize=(12, 8))
    shap.waterfall_plot(exp, max_display=10, show=False)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/figures/SHAP_Waterfall_{label}.tiff", dpi=300, bbox_inches='tight', format='tiff')
    plt.close()
    print(f"  瀑布图 ({label}): {OUTPUT_DIR}/figures/SHAP_Waterfall_{label}.tiff")

# ============================================================================
# 10. 图4: 全部33个基因的摘要图 - TIFF格式
# ============================================================================

plt.figure(figsize=(14, 12))
shap.summary_plot(shap_values, X_scaled, feature_names=gene_names_clean,
                  show=False, max_display=33)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/SHAP_Summary_All33.tiff", dpi=300, bbox_inches='tight', format='tiff')
plt.close()
print(f"  全部摘要图: {OUTPUT_DIR}/figures/SHAP_Summary_All33.tiff")

# ============================================================================
# 11. 输出摘要
# ============================================================================

print("\n" + "=" * 60)
print("✅ SHAP分析完成！")
print("=" * 60)
print(f"\n输出目录: {OUTPUT_DIR}")
print("\n生成的文件:")
print("  表格:")
print(f"    - {OUTPUT_DIR}/tables/SHAP_Top10_genes.csv")
print(f"    - {OUTPUT_DIR}/tables/SHAP_All_33genes.csv")
print(f"    - {OUTPUT_DIR}/tables/SHAP_Top10.tex")
print("  图片 (TIFF格式):")
print(f"    - {OUTPUT_DIR}/figures/SHAP_Bar_Top10.tiff")
print(f"    - {OUTPUT_DIR}/figures/SHAP_Bubble_Top10.tiff")
print(f"    - {OUTPUT_DIR}/figures/SHAP_Waterfall_Positive.tiff")
print(f"    - {OUTPUT_DIR}/figures/SHAP_Waterfall_Negative.tiff")
print(f"    - {OUTPUT_DIR}/figures/SHAP_Summary_All33.tiff")
print("\n【Top 10 SHAP基因】")
print(top10[['Rank', 'Gene', 'Mean_ABS_SHAP']].to_string(index=False))

