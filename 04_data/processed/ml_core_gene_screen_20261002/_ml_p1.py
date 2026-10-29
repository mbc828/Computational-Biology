# -*- coding: utf-8 -*-
"""Part 1: 数据准备 + 三方法全数据特征选择, 结果存盘"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

BASE = Path(r'<LOCAL_DRIVE>/kimi\kimi\workspace\SHXD-PJI项目')
OUT = BASE / 'ml_core_gene_screen_20261002'
OUT.mkdir(exist_ok=True)
TOP_VAR = 3000
RFE_N = 30
SEED = 20261002

counts = pd.read_csv(BASE / 'GSE255786_DE' / 'GSE255786_gene_counts.txt', sep='\t')
counts = counts.groupby('GeneId', as_index=False).sum(numeric_only=True)
labels = pd.read_csv(BASE / 'GSE255786_DE' / 'sample_groups.csv')
labels['sample'] = labels['sample'].astype(str)
samples = [c for c in counts.columns if c in set(labels['sample'])]
y = labels.set_index('sample').loc[samples, 'group'].map({'PJI': 1, 'AF': 0}).values
X_raw = counts.set_index('GeneId')[samples].astype(float)
cpm = X_raw / X_raw.sum(axis=0) * 1e6
keep = (cpm > 1).sum(axis=1) >= int(0.2 * len(samples))
Xdf = np.log2(cpm[keep] + 0.5).T
genes = np.array(Xdf.columns)
X = Xdf.values.astype(np.float32)
np.save(OUT / 'X.npy', X); np.save(OUT / 'y.npy', y); np.save(OUT / 'genes.npy', genes)
print(f'样本 {X.shape[0]} (PJI={y.sum()}, AF={(y==0).sum()}), 基因 {X.shape[1]}', flush=True)

vidx = np.argsort(X.var(axis=0))[::-1][:TOP_VAR]
sc = StandardScaler().fit(X[:, vidx])
Xs = np.zeros_like(X); Xs[:, vidx] = sc.transform(X[:, vidx])

# LASSO 选C (10个候选 × 5折, 单轮)
Cs = np.logspace(-2.5, -0.3, 10)
inner = StratifiedKFold(5, shuffle=True, random_state=SEED)
best_c, best_auc = None, -1
for c in Cs:
    auc = cross_val_score(LogisticRegression(penalty='l1', solver='liblinear', C=c, max_iter=3000, random_state=SEED),
                          Xs, y, cv=inner, scoring='roc_auc', n_jobs=-1).mean()
    print(f'  C={c:.4f} AUC={auc:.3f}', flush=True)
    if auc > best_auc: best_auc, best_c = auc, c
lasso = LogisticRegression(penalty='l1', solver='liblinear', C=best_c, max_iter=5000, random_state=SEED).fit(Xs, y)
lasso_genes = genes[lasso.coef_[0] != 0]
print(f'LASSO: C={best_c}, 内层AUC={best_auc:.3f}, 选中 {len(lasso_genes)}', flush=True)

rfe = RFE(SVC(kernel='linear', C=1.0, random_state=SEED), n_features_to_select=RFE_N, step=300).fit(Xs, y)
svm_genes = genes[rfe.support_]
print(f'SVM-RFE: 选中 {len(svm_genes)}', flush=True)

rf = RandomForestClassifier(n_estimators=1000, random_state=SEED, n_jobs=-1).fit(Xs, y)
rf_genes = genes[np.argsort(rf.feature_importances_)[::-1][:RFE_N]]
print(f'RF: 前 {len(rf_genes)}', flush=True)

core = sorted(set(lasso_genes) & set(svm_genes) & set(rf_genes))
print('三方法交集:', core, flush=True)

imp = pd.DataFrame({'gene': genes, 'rf_importance': rf.feature_importances_, 'lasso_coef': lasso.coef_[0]})
imp.to_csv(OUT / 'gene_importance_full.csv', index=False)
pd.DataFrame({'LASSO': pd.Series(lasso_genes), 'SVM_RFE': pd.Series(svm_genes), 'RF_top30': pd.Series(rf_genes)}).to_csv(OUT / 'selected_genes_by_method.csv', index=False)
json.dump({'best_c': float(best_c), 'lasso_genes': list(lasso_genes), 'svm_genes': list(svm_genes),
           'rf_genes': list(rf_genes), 'core': core, 'top_var': TOP_VAR},
          open(OUT / 'selection.json', 'w'), ensure_ascii=False, indent=1)
print('Part1 完成', flush=True)
