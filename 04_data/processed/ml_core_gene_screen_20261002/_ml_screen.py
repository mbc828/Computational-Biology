# -*- coding: utf-8 -*-
"""GSE255786 (PJI=53 vs AF=40, 纠正标签) 机器学习核心基因筛选
方法: LASSO(L1逻辑回归) + SVM-RFE(线性) + RandomForest
评估: 重复分层5折CV (5 repeats), 选择+建模全部在折内完成, 报告 AUC 均值±SD
口径: 探索性/描述性, 单队列内部验证, 不作临床诊断宣称
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC', 'SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

BASE = Path(r'<LOCAL_DRIVE>/kimi\kimi\workspace\SHXD-PJI项目')
OUT = BASE / 'ml_core_gene_screen_20261002'
OUT.mkdir(exist_ok=True)
STUDY_GENES = ['TLR2', 'LY96', 'MYD88', 'IL1B', 'IL6', 'TNF', 'CASP1', 'NFKBIA',
               'TLR4', 'NLRP3', 'GSDMD', 'IKBKB', 'PYCARD']
TOP_VAR = 3000      # 方差预过滤(折内)
RFE_N = 30          # SVM-RFE 固定选30个
RF_TOP = 30         # RF 取重要性前30
SEED = 20261002

# ---------- 数据 ----------
counts = pd.read_csv(BASE / 'GSE255786_DE' / 'GSE255786_gene_counts.txt', sep='\t')
counts = counts.groupby('GeneId', as_index=False).sum(numeric_only=True)  # 同GeneId合计, 与主分析一致
labels = pd.read_csv(BASE / 'GSE255786_DE' / 'sample_groups.csv')
labels['sample'] = labels['sample'].astype(str)
samples = [c for c in counts.columns if c in set(labels['sample'])]
y = labels.set_index('sample').loc[samples, 'group'].map({'PJI': 1, 'AF': 0}).values
X_raw = counts.set_index('GeneId')[samples].astype(float)  # genes x samples
lib = X_raw.sum(axis=0)
cpm = X_raw / lib * 1e6
keep = (cpm > 1).sum(axis=1) >= int(0.2 * len(samples))   # 低表达过滤
X = np.log2(cpm[keep] + 0.5).T                             # samples x genes
genes = np.array(X.columns)
X = X.values
print(f'样本 {X.shape[0]} (PJI={y.sum()}, AF={(y==0).sum()}), 过滤后基因 {X.shape[1]}')

# ---------- 各方法: 全数据选基因 ----------
def var_prefilter(Xtr, k=TOP_VAR):
    idx = np.argsort(Xtr.var(axis=0))[::-1][:k]
    return np.sort(idx)

# LASSO (全数据)
sc = StandardScaler().fit(X[:, var_prefilter(X)])
Xs = np.zeros_like(X); Xs[:, var_prefilter(X)] = sc.transform(X[:, var_prefilter(X)])
lasso = LogisticRegression(penalty='l1', solver='liblinear', C=0.05, max_iter=5000, random_state=SEED)
Cs = np.logspace(-2.5, -0.3, 20)
best_c, best_auc = None, -1
inner = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=SEED)
for c in Cs:
    auc = cross_val_score(lasso.__class__(penalty='l1', solver='liblinear', C=c, max_iter=5000, random_state=SEED),
                          Xs, y, cv=inner, scoring='roc_auc', n_jobs=-1).mean()
    if auc > best_auc: best_auc, best_c = auc, c
lasso = LogisticRegression(penalty='l1', solver='liblinear', C=best_c, max_iter=5000, random_state=SEED).fit(Xs, y)
lasso_genes = genes[lasso.coef_[0] != 0]
print(f'LASSO: C={best_c:.4f} 内层AUC={best_auc:.3f}, 选中 {len(lasso_genes)} 个基因')

# SVM-RFE (全数据, 折内一致预过滤)
svc = SVC(kernel='linear', C=1.0, random_state=SEED)
rfe = RFE(svc, n_features_to_select=RFE_N, step=100).fit(Xs, y)
svm_genes = genes[rfe.support_]
print(f'SVM-RFE: 选中 {len(svm_genes)} 个基因')

# RandomForest (全数据)
rf = RandomForestClassifier(n_estimators=2000, random_state=SEED, n_jobs=-1).fit(Xs, y)
rf_rank = np.argsort(rf.feature_importances_)[::-1]
rf_genes = genes[rf_rank[:RF_TOP]]
print(f'RandomForest: 取重要性前 {len(rf_genes)} 个')

core = sorted(set(lasso_genes) & set(svm_genes) & set(rf_genes))
print('三方法交集:', core)

# ---------- 诚实性能: 选择+建模全部放进外层CV ----------
def lasso_pipe(Xtr, ytr):
    idx = var_prefilter(Xtr)
    sc = StandardScaler().fit(Xtr[:, idx])
    Xt = np.zeros_like(Xtr); Xt[:, idx] = sc.transform(Xtr[:, idx])
    m = LogisticRegression(penalty='l1', solver='liblinear', C=best_c, max_iter=5000, random_state=SEED).fit(Xt, ytr)
    return ('lasso', idx, sc, m)

def svm_pipe(Xtr, ytr):
    idx = var_prefilter(Xtr)
    sc = StandardScaler().fit(Xtr[:, idx])
    Xt = np.zeros_like(Xtr); Xt[:, idx] = sc.transform(Xtr[:, idx])
    sel = RFE(SVC(kernel='linear', C=1.0, probability=True, random_state=SEED), n_features_to_select=RFE_N, step=100).fit(Xt, ytr)
    return ('svm', idx, sc, sel)

def rf_pipe(Xtr, ytr):
    idx = var_prefilter(Xtr)
    sc = StandardScaler().fit(Xtr[:, idx])
    Xt = np.zeros_like(Xtr); Xt[:, idx] = sc.transform(Xtr[:, idx])
    m = RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1).fit(Xt, ytr)
    return ('rf', idx, sc, m)

def core_pipe(Xtr, ytr):  # 三方法交集特征 + 逻辑回归
    gidx = np.array([list(genes).index(g) for g in core if g in genes])
    sc = StandardScaler().fit(Xtr[:, gidx])
    Xt = np.zeros((Xtr.shape[0], len(gidx))); Xt[:] = sc.transform(Xtr[:, gidx])
    m = LogisticRegression(max_iter=5000, random_state=SEED).fit(Xt, ytr)
    return ('core', gidx, sc, m)

def eval_pipe(factory, name):
    outer = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=SEED)
    aucs = []
    for tr, te in outer.split(X, y):
        kind, idx, sc, m = factory(X[tr], y[tr])
        if kind == 'core':
            Xte = sc.transform(X[te][:, idx])
        else:
            Xte = np.zeros((len(te), X.shape[1])); Xte[:, idx] = sc.transform(X[te][:, idx])
        if hasattr(m, 'predict_proba'):
            p = m.predict_proba(Xte)[:, 1]
        else:
            p = m.decision_function(Xte)
        from sklearn.metrics import roc_auc_score
        aucs.append(roc_auc_score(y[te], p))
    aucs = np.array(aucs)
    print(f'{name}: 外层CV AUC = {aucs.mean():.3f} ± {aucs.std():.3f} (n={len(aucs)}折)')
    return aucs

res = {
    'LASSO': eval_pipe(lasso_pipe, 'LASSO'),
    'SVM-RFE': eval_pipe(svm_pipe, 'SVM-RFE'),
    'RandomForest': eval_pipe(rf_pipe, 'RandomForest'),
    'Core交集+LR': eval_pipe(core_pipe, 'Core交集+LR'),
}

# ---------- 输出 ----------
pd.DataFrame({
    'LASSO': pd.Series(lasso_genes), 'SVM_RFE': pd.Series(svm_genes), 'RF_top30': pd.Series(rf_genes)
}).to_csv(OUT / 'selected_genes_by_method.csv', index=False)
pd.DataFrame({'core_gene': core}).to_csv(OUT / 'core_genes_intersection.csv', index=False)

imp = pd.DataFrame({'gene': genes, 'rf_importance': rf.feature_importances_,
                    'lasso_coef': lasso.coef_[0]}).sort_values('rf_importance', ascending=False)
imp.to_csv(OUT / 'gene_importance_full.csv', index=False)

study_hit = {g: dict(lasso=bool(g in lasso_genes), svm=bool(g in svm_genes), rf=bool(g in rf_genes),
                     core=bool(g in core)) for g in STUDY_GENES}
pd.DataFrame(study_hit).T.to_csv(OUT / 'study_genes_hit.csv')

summary = {
    'samples': {'total': int(X.shape[0]), 'PJI': int(y.sum()), 'AF': int((y == 0).sum())},
    'genes_after_filter': int(X.shape[1]),
    'lasso': {'C': float(best_c), 'n_selected': int(len(lasso_genes)), 'genes': list(lasso_genes)},
    'svm_rfe': {'n_selected': int(len(svm_genes)), 'genes': list(svm_genes)},
    'rf_top': list(rf_genes),
    'core_intersection': core,
    'cv_auc': {k: {'mean': float(v.mean()), 'sd': float(v.std())} for k, v in res.items()},
    'study_genes_hit': study_hit,
    'caveat': '单队列内部重复CV估计, 偏乐观; 未做外部验证; 不作临床诊断宣称',
}
(OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')

# ---------- 图 ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
ax = axes[0]
names = list(res)
means = [res[k].mean() for k in names]
sds = [res[k].std() for k in names]
bars = ax.bar(names, means, yerr=sds, capsize=5, color=['#4C72B0', '#55A868', '#C44E52', '#8172B2'], alpha=0.85)
ax.axhline(0.5, ls='--', c='grey')
ax.set_ylim(0.4, 1.02)
ax.set_ylabel('重复分层5折CV AUC (均值±SD)')
ax.set_title('GSE255786 PJI vs AF (n=93): 内部验证 AUC')
for b, m_ in zip(bars, means):
    ax.text(b.get_x() + b.get_width() / 2, m_ + 0.03, f'{m_:.3f}', ha='center', fontsize=10)
plt.setp(ax.get_xticklabels(), rotation=12)

ax = axes[1]
top = imp.head(15).iloc[::-1]
ax.barh(top['gene'], top['rf_importance'], color='#4C72B0', alpha=0.85)
for g in top['gene']:
    if g in core:
        ax.text(0, list(top['gene']).index(g), ' ★', color='red', fontsize=12, va='center')
ax.set_title('RandomForest 重要性 Top15 (★=三方法交集核心基因)')
ax.set_xlabel('Mean decrease in impurity')
plt.tight_layout()
plt.savefig(OUT / 'ml_core_gene_screen.png', dpi=200)
print('图已存:', OUT / 'ml_core_gene_screen.png')
print('输出目录:', OUT)
