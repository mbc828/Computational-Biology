# ML 核心基因筛选结果摘要（GSE255786，PJI vs AF）

日期：2026-10-02 ｜ 数据：GSE255786 纠正标签（PJI=53，AF=40，假体超声处理液，18,192 个过滤后基因）
目的：借鉴 Phytomedicine 2026;160:158683 的"ML 收敛核心基因"思路，检验本文预设通路节点是否也被数据驱动的方法选中。

## 方法

- 预处理：同 GeneId 合计 → CPM>1 且覆盖 ≥20% 样本过滤 → log2-CPM → z 标准化
- 三方法：LASSO（L1 逻辑回归，C 网格内层 5 折选定）、SVM-RFE（线性核选 30 个）、RandomForest（1000 树取前 30）；另跑无方差预过滤的全特征变体（LASSO+RF）验证稳健性
- 性能评估：重复分层 5 折 CV ×3（15 折），**特征选择+建模全部在折内完成**，报告 AUC 均值±SD

## 结果

**判别性能（诚实口径，折内选择）**：

| 模型 | 重复 CV AUC |
|---|---|
| LASSO | 0.991 ± 0.019 |
| SVM-RFE | 0.984 ± 0.032 |
| RandomForest（全特征） | 0.989 ± 0.022 |

PJI 与 AF 在转录层高度可分，三种方法一致。

**核心基因**：
- LASSO 仅选 4 个：EFNA2、NKX3-1、RPS10P7、S1PR2（全特征变体结果完全相同，稳健）
- ≥2 方法共识基因 19 个：EFNA2、S1PR2、NKX3-1、RPS10P7、HAS2、SLC6A14、SEMA3A、PVR、RFPL1-AS1、GPM6A、PAN3、C10orf32、IBA57、RAB2B、FKBP9L、CLEC4GP1、ARSI、LOC100129858、LOC100131208
- 三方法完全交集为空；LASSO∩RF（全特征）= EFNA2、NKX3-1、S1PR2
- 共识 19 基因 + 逻辑回归 AUC 0.997 ± 0.009 —— **注意：该数值的特征是在全体数据上预先选出的，存在选择泄漏，仅示意，不作报告值**

**关键阴性结果**：预设通路基因无一入选。全特征 RF 排名（/18,192）：IL1B 168、TNF 754；TLR2、LY96、MYD88、IL6、CASP1、NFKBIA、TLR4、NLRP3、GSDMD、IKBKB、PYCARD 均约 9,953 名（重要性≈0）。NLRP3/TLR4/IKKβ 轴不构成 ML 判别信号的主要来源。

## 对论文的意义

1. **可以借鉴的叙事成立一半**：ML 确实能给出高 AUC 的判别签名和可聚焦的候选基因集（EFNA2、S1PR2 等），"单队列内部验证 AUC≈0.99"可写入结果，但**必须标注**单队列、内部 CV、偏乐观、无外部验证。
2. **不能学对方的"ML 收敛到机制核心基因"**：本数据里 ML 没有收敛到 NLRP3/TLR4/IKKβ 轴——这与稿件既有口径反而自洽（bulk 差异含细胞组成/临床协变量影响，§四.4 已声明）。若强行把 ML 和计算靶点缝合，会制造新的不一致。
3. **建议写法**：新增一小节"机器学习判别签名"作为独立、描述性结果（判别能力 + 共识基因清单），并明确写"判别签名不覆盖预设通路节点，提示 PJI–AF 可分性来自更广泛的转录重塑而非单一炎症轴"——把阴性结果写成对整体谨慎口径的支持。
4. 若要追求对方那种"ML→机制"闭环，路径是 WB 阶段 B 同批样本的 RNA-seq（药物干预侧签名）+ LINCS 回投，而不是在疾病侧数据上硬凑。

## 文件清单

- `summary.json`、`selection.json`、`selection_fullfeature.json`：全部数值结果
- `selected_genes_by_method.csv`、`core_genes_intersection.csv`（经 selection.json）、`gene_importance_full.csv`、`rf_importance_fullfeature.csv`
- `ml_core_gene_screen.png`：AUC 条形图 + RF 重要性 Top15
- `_ml_p1.py / _ml_p1b.py / _ml_p2.py`：完整可复跑脚本（`_ml_screen.py` 为初版弃稿）
