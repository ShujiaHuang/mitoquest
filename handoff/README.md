# MitoQuest 开发历史存档

本目录按时间顺序记录了 MitoQuest 项目从 v1.3.0 到 v1.10.0 的重要版本更新、技术决策和实现细节。

> Promot: 我要求你自动生成一份包含当前进度、决策和下一步计划的 HANDOFF.md 存档文件，文件名最好包括有用的信息和时间信息，然后将它存到 ​handoff​ 目录下。

```txt
每份文档均包含：
- 阶段总览（版本/日期/变更摘要表）
- 核心功能与技术决策（问题背景、解决方案、设计理由）
- 已知限制（遗留至下一阶段的问题）
- 关键文件索引
- 下一阶段入口（互相链接形成完整时间线）
```

---

## 存档索引

| 文档 | 版本范围 | 时间 | 主题 |
|---|---|---|---|
| [HANDOFF_v1.3.0-v1.5.5_early_development.md](HANDOFF_v1.3.0-v1.5.5_early_development.md) | v1.3.0 ~ v1.5.5 | 2025-03 ~ 2025-08 | 初期变异调用器：Phred p-value、FS/SB 修复、Indel 检测 |
| [HANDOFF_v1.6.0-v1.6.4_build_infrastructure.md](HANDOFF_v1.6.0-v1.6.4_build_infrastructure.md) | v1.6.0 ~ v1.6.4 | 2026-01 ~ 2026-05 | 构建系统稳定化：CMake 修复、CI 搭建、静态二进制 |
| [HANDOFF_v1.7.0-v1.7.2_copynum_subcommand.md](HANDOFF_v1.7.0-v1.7.2_copynum_subcommand.md) | v1.7.0 ~ v1.7.2 | 2026-05-22 | mtcopynum 集成为 `copynum` 子命令；NUMT 区域排除 |
| [HANDOFF_v1.8.0-v1.8.6_ne_estimation_pipeline.md](HANDOFF_v1.8.0-v1.8.6_ne_estimation_pipeline.md) | v1.8.0 ~ v1.8.6 | 2026-05-28 ~ 05-30 | Ne 估计流水线：trans-prep/ne-estimate、连续模型、bin-simulation、MMLE 术语 |
| [HANDOFF_v1.9.0-v1.9.2_variant_qc_kimura_ux.md](HANDOFF_v1.9.0-v1.9.2_variant_qc_kimura_ux.md) | v1.9.0 ~ v1.9.2 | 2026-06-01 ~ 06-02 | variant-qc C++ 重写；Kimura trim/outlier UX 修复 |
| [HANDOFF_v1.10.0_20250528.md](HANDOFF_v1.10.0_20250528.md) | v1.10.0 | 2026-06-02 | 三代 G-M-C trio 边际似然；`--gm-fam` 参数；闭式 Beta 函数比 |
| [HANDOFF_per_family_ne_estimation.md](HANDOFF_per_family_ne_estimation.md) | 设计文档 | 2026-06-11 | Per-Family $N_e$ 估计：MMLE + 校正 Kimura 矩估计 + 两步残差法 (vGWAS) + 三篇文献审查 + 数学查漏补缺 + 详细代码升级方案 |
| [HANDOFF_ne_estimate_methods.md](HANDOFF_ne_estimate_methods.md) | 方法学参考 | 2026-06-11 | `ne-estimate` 完整方法学手册：公式推导、代码实现、计算实例 |
| [HANDOFF_ne_estimate_vs_deCODE_methodology.md](HANDOFF_ne_estimate_vs_deCODE_methodology.md) | 方法学对比 | 2026-06-13 | mitoquest `ne-estimate` vs. deCODE (Helgason et al., Cell 2024) 方法论全面对比 |

---

## 版本时间线

```txt
2025-03-10  v1.3.0   项目初始化
2025-03-14  v1.4.0   Phred p-value 异质性检验
2025-03-14  v1.4.1   FS/SB 计算修复
2025-03-25  v1.5.0   Log-odds ratio + Indel 修复
2025-04-10  v1.5.2   Linux 静态二进制
2025-08-19  v1.5.5   稳定性修复
─────────────────────────────────────────────
2026-01-28  v1.6.0   Het/Hom 状态 + 变异类型
2026-03-30  v1.6.1   双平台二进制发布
2026-05-22  v1.6.4   CMake 重大修复 + CI
2026-05-22  v1.7.0   copynum 子命令集成
2026-05-22  v1.7.1   -L/--regions NUMT 排除
─────────────────────────────────────────────
2026-05-28  v1.8.0   trans-prep + ne-estimate
2026-05-28  v1.8.1   Kimura 分母 bug 修复 + bootstrap CI
2026-05-28  v1.8.2   连续 Beta-diffusion 模型（替代离散 MLE）
2026-05-28  v1.8.3   实值 Ne 优化 + README 解释
2026-05-28  v1.8.4   per-bin drift summary TSV 输出
2026-05-30  v1.8.5   Ne-profile 诊断
2026-05-30  v1.8.6   MLE → MMLE 术语修正
─────────────────────────────────────────────
2026-06-01  v1.9.0   variant-qc 子命令（C++ 重写）
2026-06-01  v1.9.1   variant-qc 文档 + bug 修复
2026-06-02  v1.9.2   Kimura trim/outlier UX 修复
2026-06-02  v1.10.0  三代 G-M-C trio 边际似然
2026-06-11  设计    Per-Family Ne 估计方法学文档（文献审查 + 数学查漏补缺 + 详细代码升级方案）
```
