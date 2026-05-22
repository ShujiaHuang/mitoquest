# MitoQuest Changelog

## [1.7.0] - 2026-05-22

- 将原本独立的 `mtcopynum` 工具重写并融合进主程序，作为 `mitoquest copynum` 子命令调用
- 复用 `ngslib::Bam`/`ngslib::Fasta`/`ThreadPool`/`mean`/`standard_error` 等已有模块，移除重复代码
- 在 `mt_utils.h` 中新增 `is_mitochondrial()` 工具函数
- 重构 `MtCopyNumber` 类：拆分出可独立测试的 `compute_gc_content` / `compute_statistics` / `compute_normalized_ratios` 静态方法
- 改进 PE 配对去重逻辑：read1 + 孤立 mate 通过 qname 去重，避免漏计 orphan reads
- 在 `tests/mt_copynum_test.cpp` 中新增 20 个 GoogleTest 用例覆盖工具函数、统计计算、归一化、端到端 CRAM 输入
- 修复 macOS 下 `-Wl,-no_compact_unwind` 全局链接选项导致 gtest `EXPECT_THROW` 在跨翻译单元抛出异常时调用 `std::terminate` 的问题：将该选项收敛到生产可执行文件 target，不再传染到测试二进制
- 输入 BAM 不存在或 fragment 计数为 0 时改为快速失败并给出清晰错误信息

## [1.6.4] - 2026-05-21

- CMake / GitHub Actions 构建脚本若干修复

## [1.5.5] - 2025-08-19

- 修复已知问题
- 调整多线程模块，使用 `std::bind` 解决成员函数多线程调用的问题

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
