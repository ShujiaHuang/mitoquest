# MitoQuest Changelog

## [TBD] - TBD

## [1.5.4] - 2025-05-13

- 改进 VCF 处理模块
- 添加按照实际的样本，对 REF 和 ATL 中多余的基因型进行清理的功能
- INFO 中 `PT` 的标签内容改名为：`Hom`, `Het` 和 `Mixed`, 并增加一个 `Ref` 用于标记纯粹 REF 的记录
- 修复已知 bug

## [1.5.3] - 2025-05-09

- 添加一系列处理 VCF 的模块
- 添加 `subsam` 功能，用于从全体 VCF 中提取目标样本，并根据目标样本更新 VCF INFO 信息

## [1.5.2] - 2025-04-10

- 添加 `-Q` 用于指定过滤测序质量值
- INFO 添加 `HOM_AF`, `HET_AF` 和 `SUM_AF` 记录携带非 REF 的同质性、异质性、同质性+异质性的人群频率
- INFO 添加 `Plasmicity_type: 缩写PT` 用来记录该位点在群体中是 `Hom_only/Het_only/Both`

## [1.5.0] - 2025-03-25

- 修复已知 bug
- 不再移动 Indel 端点，解决 Indel 位点深度覆盖的计算和异质性估算的 bug 问题
- Indel 优先原则：输出时，如果 break point 的 leftmost 是 SNV，那么强制用 Indel 代替
- 添加 `LHF` 记录异质性分值的对数几率模型(对数优势比)

## [1.4.1] - 2025-03-14

- 添加 `HQ` 记录异质性碱基可信值
- 使用 `Fisher's exact test` 计算 `phred-score` 描述异质性碱基的可靠性
- 修正 SB 的计算，并非都是依据 `REF`，而是要改为 `major allele`
- 重写并完善了 `iobgzf.h` 中的 `BGZFile class`，实现对 `htslib/bgzf.h` 的封装
- 使用 `BGZFile class` 不再直接调用 BGZF，这样更安全、更方便，也符合 RAII 原则
- 优化 `fisher_exact_test` 等多个统计学函数

## [1.3.0] - 2025-03-10

- 创建 Github Actions workflows 实现跨平台编译
- 发布 Release 时自动附带编译好的二进制可执行文件
