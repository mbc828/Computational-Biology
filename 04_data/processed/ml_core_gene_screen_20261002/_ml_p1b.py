# -*- coding: utf-8 -*-
"""Part 1b: 无方差预过滤变体 —— LASSO 与 RF 在全部 18192 基因上选"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

OUT = Path(r'<LOCAL_DRIVE>/kimi\kimi\workspace\SHXD-PJI项目\ml_core_gene_screen_20261002')
X = np.load(OUT / 'X.npy'); y = np.load(OUT / 'y.npy'); genes = np.load(OUT / 'genes.npy', allow_pickle=True)
SEED = 20261002

sc = StandardScaler().fit(X)
Xs = sc.transform(X)

# LASSO 全特征 (C 沿用 part1 网格最优附近, 直接复用 best_c)
sel = json.load(open(OUT / 'selection.json'))
best_c = sel['best_c']
lasso_all = LogisticRegression(penalty='l1', solver='liblinear', C=best_c, max_iter=8000, random_state=SEED).fit(Xs, y)
lasso_all_genes = genes[lasso_all.coef_[0] != 0]
print(f'LASSO(全特征): 选中 {len(lasso_all_genes)}', flush=True)
print(list(lasso_all_genes), flush=True)

rf_all = RandomForestClassifier(n_estimators=1000, random_state=SEED, n_jobs=-1).fit(Xs, y)
rf_all_genes = genes[np.argsort(rf_all.feature_importances_)[::-1][:30]]
print('RF(全特征) Top30:', list(rf_all_genes), flush=True)

# 主线基因在全特征RF中的排名
study = ['TLR2', 'LY96', 'MYD88', 'IL1B', 'IL6', 'TNF', 'CASP1', 'NFKBIA', 'TLR4', 'NLRP3', 'GSDMD', 'IKBKB', 'PYCARD']
rank = pd.Series(rf_all.feature_importances_, index=genes).rank(ascending=False)
imp = pd.Series(rf_all.feature_importances_, index=genes)
print('\n主线基因(全特征RF):')
for g in study:
    print(f'  {g:8s} 排名 {int(rank[g]):6d} / {len(genes)}, 重要性 {imp[g]:.5f}', flush=True)

core_b = sorted(set(lasso_all_genes) & set(rf_all_genes))
print('LASSO∩RF(全特征):', core_b, flush=True)

json.dump({'lasso_all_genes': list(lasso_all_genes), 'rf_all_top30': list(rf_all_genes),
           'lasso_rf_overlap': core_b,
           'study_genes_rf_all': {g: {'rank': int(rank[g]), 'importance': float(imp[g])} for g in study}},
          open(OUT / 'selection_fullfeature.json', 'w'), ensure_ascii=False, indent=1)
imp.rename('rf_importance_full').to_csv(OUT / 'rf_importance_fullfeature.csv')
print('Part1b 完成', flush=True)
