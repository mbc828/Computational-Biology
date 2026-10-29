# -*- coding: utf-8 -*-
"""Part 2: 外层重复CV评估(选择+建模全在折内) + 出图 + 汇总"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC']
plt.rcParams['axes.unicode_minus'] = False

OUT = Path(r'<LOCAL_DRIVE>/kimi\kimi\workspace\SHXD-PJI项目\ml_core_gene_screen_20261002')
X = np.load(OUT / 'X.npy'); y = np.load(OUT / 'y.npy'); genes = np.load(OUT / 'genes.npy', allow_pickle=True)
sel = json.load(open(OUT / 'selection.json'))
sel_b = json.load(open(OUT / 'selection_fullfeature.json'))
SEED = 20261002
TOP_VAR = 3000

def fold_eval(fit_predict, label, repeats=3):
    outer = RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=SEED)
    aucs = []
    for tr, te in outer.split(X, y):
        p = fit_predict(X[tr], y[tr], X[te])
        aucs.append(roc_auc_score(y[te], p))
    aucs = np.array(aucs)
    print(f'{label}: AUC {aucs.mean():.3f} ± {aucs.std():.3f} ({len(aucs)}折)', flush=True)
    return aucs

def vidx_of(Xtr):
    return np.sort(np.argsort(Xtr.var(axis=0))[::-1][:TOP_VAR])

def fp_lasso(Xtr, ytr, Xte):
    vi = vidx_of(Xtr)
    sc = StandardScaler().fit(Xtr[:, vi])
    m = LogisticRegression(penalty='l1', solver='liblinear', C=sel['best_c'], max_iter=3000, random_state=SEED)
    m.fit(sc.transform(Xtr[:, vi]), ytr)
    return m.predict_proba(sc.transform(Xte[:, vi]))[:, 1]

def fp_svm(Xtr, ytr, Xte):
    vi = vidx_of(Xtr)
    sc = StandardScaler().fit(Xtr[:, vi])
    sel_rfe = RFE(SVC(kernel='linear', C=1.0, probability=True, random_state=SEED), n_features_to_select=30, step=500)
    sel_rfe.fit(sc.transform(Xtr[:, vi]), ytr)
    return sel_rfe.predict_proba(sc.transform(Xte[:, vi]))[:, 1]

def fp_rf(Xtr, ytr, Xte):
    sc = StandardScaler().fit(Xtr)
    m = RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1)
    m.fit(sc.transform(Xtr), ytr)
    return m.predict_proba(sc.transform(Xte))[:, 1]

# 两票基因(至少两个方法选中): LASSO∪变体 与 SVM/RF 的交集
votes = {}
for g in sel['lasso_genes'] + sel['svm_genes'] + sel['rf_genes'] + sel_b['lasso_all_genes'] + sel_b['rf_all_top30']:
    votes[g] = votes.get(g, 0) + 1
two_plus = sorted([g for g, v in votes.items() if v >= 2])
print('≥2票基因:', two_plus, flush=True)
gidx = np.array([list(genes).index(g) for g in two_plus])

def fp_consensus(Xtr, ytr, Xte):
    sc = StandardScaler().fit(Xtr[:, gidx])
    m = LogisticRegression(max_iter=3000, random_state=SEED)
    m.fit(sc.transform(Xtr[:, gidx]), ytr)
    return m.predict_proba(sc.transform(Xte[:, gidx]))[:, 1]

res = {
    'LASSO': fold_eval(fp_lasso, 'LASSO'),
    'SVM-RFE': fold_eval(fp_svm, 'SVM-RFE'),
    'RandomForest': fold_eval(fp_rf, 'RandomForest(全特征)'),
    '共识基因+LR': fold_eval(fp_consensus, f'共识{len(two_plus)}基因+LR'),
}

# 图
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
ax = axes[0]
names = list(res)
means = [res[k].mean() for k in names]; sds = [res[k].std() for k in names]
ax.bar(names, means, yerr=sds, capsize=5, color=['#4C72B0', '#55A868', '#C44E52', '#8172B2'], alpha=0.85)
ax.axhline(0.5, ls='--', c='grey'); ax.set_ylim(0.4, 1.02)
ax.set_ylabel('AUC (均值±SD, 重复5折CV)')
ax.set_title('GSE255786 PJI(n=53) vs AF(n=40): 内部验证')
for i, m_ in enumerate(means):
    ax.text(i, m_ + 0.03, f'{m_:.3f}', ha='center', fontsize=10)
plt.setp(ax.get_xticklabels(), rotation=12)

ax = axes[1]
imp = pd.Series(np.load(OUT / 'X.npy')[0] * 0, index=genes)  # placeholder
rfimp = pd.read_csv(OUT / 'rf_importance_fullfeature.csv', index_col=0).iloc[:, 0]
top = rfimp.sort_values(ascending=False).head(15).iloc[::-1]
colors = ['#C44E52' if g in two_plus else '#4C72B0' for g in top.index]
ax.barh(top.index, top.values, color=colors, alpha=0.85)
ax.set_title('RF(全特征)重要性Top15 (红=≥2方法共识)')
ax.set_xlabel('importance')
plt.tight_layout()
plt.savefig(OUT / 'ml_core_gene_screen.png', dpi=200)

summary = {
    'cohort': 'GSE255786 纠正标签, PJI=53 vs AF=40, 假体超声处理液',
    'genes_after_filter': int(X.shape[1]),
    'methods': 'LASSO / SVM-RFE(30) / RandomForest; 变体: LASSO+RF全特征',
    'lasso_genes': sel['lasso_genes'],
    'consensus_2of3': two_plus,
    'cv_auc': {k: {'mean': float(v.mean()), 'sd': float(v.std()), 'folds': len(v)} for k, v in res.items()},
    'pathway_genes_not_selected': True,
    'study_genes_rf_rank': sel_b['study_genes_rf_all'],
    'caveat': '单队列内部重复CV, 偏乐观; 无外部验证; 不作诊断宣称; ML判别信号未收敛到NLRP3/TLR4/IKKβ轴基因',
}
(OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')
print('图与汇总已存:', OUT, flush=True)
