在进行母婴线粒体基因组学研究时，**APOGEE2** 和 **MitoTIP** 是目前国际公认的线粒体致病性评估核心评分系统。
它们分别针对线粒体基因组的不同功能区域，是量化“突变负荷”时进行加权分析的必要工具。

以下是这两个分数的详细定义、判定标准及获取渠道

### 1. APOGEE2 (Automated Pathogenicity Given Evolutionary Evidence 2)

* **定义：** 专门针对线粒体 **13 个蛋白质编码基因**中的**错义变异（Missense variants）**开发的集成机器学习预测模型。
* **核心机制：** 采用 KNN RusSmote 机器学习算法，通过整合演化保守性、蛋白质结构改变、氨基酸物理化学性质等 13 个独立预测因子及 6 个元预测因子进行评估。
* **判定标准（原始分数范围： 0 至 1）：**
* **Likely Pathogenic (可能致病)：** 原始分数 > 0.75。
* **Possibly Pathogenic (可疑致病)：** 分数区间为 (0.50, 0.75]。
* **Neutral/Likely Benign (中性/可能良性)：** 分数 < 0.25。

* **数据获取渠道：**
* **官方网站：** [MitImpact](https://mitimpact.css-mendel.it/)。
* **下载方式：** 您可以从 MitImpact 顶部的“Download”菜单批量获取全基因组蛋白质编码区变异的预计算分数。目前已更新至 2024.0.1 版本，支持 （针对 tRNA 的变体）。


### 2. MitoTIP (Mitochondrial tRNA Informatics Predictor)

* **定义：** 专门针对线粒体 **22 个 tRNA 基因**开发的致病性预测算法。
* **核心机制：** 基于 tRNA 的二级结构（如茎环结构）、物种间序列保守性、变异发生的结构位点重要性以及历史致病位点分布进行加权评分。
* **判定标准（分数范围：-5.9 至 21.8）：**
* **Likely Pathogenic (可能致病)：** 分数 > 16.25（位于前 25% 分位数 Q1）。
* **Possibly Pathogenic (可疑致病)：** 分数区间为 [12.66, 16.25] (Q2)。
* **Possibly Benign (可能良性)：** 分数区间为 $。
* **Likely Benign (极可能良性)：** 分数 < 8.44 (Q4)。

* **数据获取渠道：**
* **官方数据库：**([https://www.mitomap.org/](https://www.mitomap.org/))。
* **直接下载链接：**([https://mitomap.org/downloads/mitotip_scores.txt](https://mitomap.org/downloads/mitotip_scores.txt))。
* **集成查询：** 也可以通过 MITOMAP 的 **Mitomaster** 接口进行在线变异注释。


### 注释

在研究中，建议建立如下的功能注释矩阵（Functional Annotation Matrix）：

| 变异类型 | 推荐工具 | 获取来源 | 顶刊引用建议 |
| --- | --- | --- | --- |
| 蛋白质编码区 (Protein-coding) | APOGEE2 | MitImpact | 引用 *Nature Communications (2023)* |
| 转运 RNA 区 (tRNA) | MitoTIP | MITOMAP | 引用 *Nucleic Acids Research (2020)* |
| 总体变异频率 (Population AF) | NMVR / Helix / gnomAD | 对应官方资源 | 重点对比 NyuWa 中国人群数据库 |



