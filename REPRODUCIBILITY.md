# 复现与数据说明

## 公开可直接核对的内容

- `02_figures/figure_manifest.csv`：24 幅图与 v95 内部媒体、像素尺寸、图注和 SHA-256 的对应关系。
- `04_data/processed/tables_data/`：转录组、外部队列、结构指标、能量与残基分析的处理后表格。
- `04_data/processed/docking/extended_docking_parameters_v95.csv`：IL-17A 与 CASP1 的 PDB 编号、受体链、共晶配体、网格、搜索参数及 Rank-1 重对接 RMSD；各项均从原始 PDB、Vina 日志和 RMSD 记录核对。
- `04_data/processed/wetlab/`：抗菌统计工作簿、相关汇总表及图 6E 检测限报告表 `antibacterial_gradient_reporting_v94.csv`。
- `04_data/processed/ml_core_gene_screen_20261002/`：机器学习结果摘要与图件。
- `04_data/processed/admet_xval/`：计算毒理交叉核对结果。
- `04_data/md_parameters/`：项目中收集到的参数、拓扑文件和逐体系参数映射，按原始相对路径保存。
- `manuscript/SHXD_PJI_v95_text.md`：当前文稿的可检索文本与图件链接。

## 尚需配置的内容

公开代码不是无配置的一键管线。需要配置本地路径、软件环境和原始输入。公共原始数据按论文中的 GEO、PDB、ProteomeXchange 编号获取；完整 MD 轨迹按数据可得性声明获取。

## 已完成与尚待补齐的记录

逐体系 MD 参数映射已完成：85 条正式轨迹、29 个体系已建立“体系—重复—参数文件”对照。WB 文件级全膜—灰度对应已完成，21 张主用全膜及相关材料已整理为本地 Source Data 投稿附件包；原始全膜和该附件包不公开上传主仓库。商业 MRSA 参考菌株 BNCC330041 的身份与相关伦理表述已在 v93 统一，不再按临床分离株追溯伦理编号。v95 已补齐五靶点重对接结果、IL-17A／CASP1 原始对接参数、附录 H4 交叉引用和图 6E 检测限表述，并将 IL-17A 的 PDB 编号统一纠正为 5HI4；单浓度抗菌臂黄芩素 24 μM 及图 6D 的 2／2／3 个视野已按分析记录补入。

网络药理学逐边来源／数据库版本、湿实验批次归属、试剂和抗体货号，以及 WB 原始 ImageJ ROI／背景扣除／曝光记录仍待作者依据原始记录补齐。抗菌实验的培养基消歧、暴露时间、vehicle、冻干粉质量口径、复溶除菌与复板属性亦须依原执行记录确认。公开包不对这些缺失信息作推测性填充。
