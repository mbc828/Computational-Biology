# 公开代码使用说明

本目录保留论文相关的转录组、结构分析、统计、制图和质量控制脚本快照。

- 原始本地绝对路径已替换为 `<PROJECT_ROOT>`、`<MD_DATA_ROOT>`、`<WETLAB_DATA_ROOT>`、`<USER_HOME>` 等占位符。
- 运行前必须按自己的目录结构修改路径，不能直接以占位符运行。
- 代码按原项目相对路径保存，便于追溯；部分脚本属于不同分析阶段，不是一个无需配置的一键管线。
- 主要依赖包括 Python、NumPy、pandas、Matplotlib、Pillow、R／Bioconductor、GROMACS、gmx_MMPBSA、AutoDock Vina 及相关结构分析工具。具体版本以论文方法部分和各运行记录为准。
- 完整 MD 轨迹、原始测序文件和私有湿实验原始记录未包含在公开仓库中。

`MANIFEST_SHA256.json` 用于核对本次公开文件是否在传输后保持不变。
