# Handoff — v1.8.0 ~ v1.8.6：Ne 估计流水线（trans-prep / ne-estimate）

> 时间跨度：2026-05-28 ~ 2026-05-30（三日内密集迭代 7 个版本）  
> 核心工作：实现 mtDNA 瓶颈大小（$N_e$）估计完整流水线；从离散 MLE 演进至连续 Beta-diffusion 模型；术语从 MLE 修正为 MMLE  
> 对应 tag：`v1.8.0` → `v1.8.1` → `v1.8.2` → `v1.8.3` → `v1.8.4` → `v1.8.5` → `v1.8.6`

---

## 1. 阶段总览

| 版本 | 日期 | 主要变更 |
|---|---|---|
| v1.8.0 | 2026-05-28 12:10 | **新增** `trans-prep` / `ne-estimate` 子命令；修复 VCF FORMAT cardinality bug |
| v1.8.1 | 2026-05-28 14:06 | **修复** Kimura 采样噪声校正分母交换 bug；新增 bootstrap CI |
| v1.8.2 | 2026-05-28 15:33 | **重大** 替换离散 MLE 为连续 Beta-diffusion 模型（修复 Ne 系统性偏高） |
| v1.8.3 | 2026-05-28 16:54 | 实值 $N_e$ 优化（golden-section search）；README 解释 Ne_MLE < Ne_Kimura 原因 |
| v1.8.4 | 2026-05-28 23:25 | 新增 per-bin drift summary TSV 输出 + 绘图工具 |
| v1.8.5 | 2026-05-30 12:11 | Ne-profile 诊断功能（MLE vs Kimura 差异分析） |
| v1.8.6 | 2026-05-30 22:43 | **术语修正** MLE → MMLE（Maximum Marginal Likelihood Estimator） |

---

## 2. v1.8.0：trans-prep / ne-estimate 流水线

### 2.1 背景与目标

**问题**：MitoQuest 能调用变异，但无法估计线粒体瓶颈大小 $N_e$。用户需要：
1. 从多样本 VCF + PLINK FAM 提取母-子等位基因传递对
2. 对这些传递对执行统计推断，估计 $N_e$ 及 95% 置信区间

### 2.2 两个新子命令

```bash
# Step 1: trans-prep — VCF+FAM → TSV
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f trios.fam \
    -o mc_pairs.tsv

# Step 2: ne-estimate — TSV → Ne 估计
mitoquest ne-estimate \
    -i mc_pairs.tsv \
    --cross-check kimura \
    -o ne_result.json
```

### 2.3 trans-prep 设计

| 输入 | 说明 |
|---|---|
| VCF（多样本） | 由 `mitoquest caller` 产生 |
| FAM（PLINK 格式） | 定义家系关系 |

**输出 TSV 格式**（v1.8.0 为 16 列）：

```
CHROM POS REF ALT FAM_ID MOTHER_ID CHILD_ID
MOTHER_DP MOTHER_AD_REF MOTHER_AD_ALT MOTHER_VAF
CHILD_DP CHILD_AD_REF CHILD_AD_ALT CHILD_VAF
QC
```

### 2.4 ne-estimate 统计模型（v1.8.0 离散版）

**离散 Beta-Binomial MLE**：对每个 mother-child pair，在离散网格 $k \in \{0, 1, \ldots, N_e\}$ 上计算条件似然，然后对所有 pair 取复合对数似然和，最大化得到 $\hat{N}_e$。

**问题发现**：离散网格在深度较高时产生系统性向上偏差（Ne_MLE=23 vs Ne_Kimura=2.78）。

### 2.5 关键文件新增

| 文件 | 用途 |
|---|---|
| `src/trans_prep.h` / `src/trans_prep.cpp` | trans-prep 子命令实现 |
| `src/ne_estimate.h` / `src/ne_estimate.cpp` | ne-estimate 子命令实现 |
| `src/log_factorial.h` / `src/log_factorial.cpp` | 从 ne_estimate 中提取的通用对数阶乘工具 |
| `tests/test_trans_prep.cpp` | trans-prep 单元测试 |
| `tests/test_ne_estimate.cpp` | ne-estimate 单元测试 |

---

## 3. v1.8.1：Kimura 采样噪声校正 Bug 修复

### 3.1 问题诊断

**现象**：真实数据上 Ne_MLE=23 vs Ne_Kimura=2.78，差距极大。

**根因**：Wonnapinij 采样噪声校正公式中分母交换：

```cpp
// BUG（错误）
p_m*(1-p_m)/c_dp + p_c*(1-p_c)/m_dp   // m_dp 和 c_dp 位置反了

// 正确
p_m*(1-p_m)/m_dp + p_c*(1-p_c)/c_dp
```

### 3.2 修复内容

- 修正 `ne_estimate.cpp` 中 `compute_kimura_check()` 的分母
- 新增 `--kimura-bootstrap` / `--kimura-seed` CLI 标志，启用非参数 bootstrap 95% CI
- JSON 输出新增 `b_CI_95_Low`、`b_CI_95_High`、`Ne_Kimura_CI_95_Low`、`Ne_Kimura_CI_95_High`

### 3.3 重要发现

修复后 Ne/Kimura 差距仍然存在。进一步分析表明这是**预期行为**：
- 离散 BB MLE 受网格分辨率限制，在高深度时向上偏
- Kimura 矩估计对高漂移离群值（NUMTs/污染）敏感，向下偏
- 两者差异反映数据异质性，不是 bug

---

## 4. v1.8.2：连续 Beta-diffusion 模型（重大重构）

### 4.1 问题

离散 MLE 在真实数据上系统性高估 $N_e$（如 Ne_MLE=23 vs Ne_Kimura=2.78），根因是离散网格强制 $p \in \{0, 1/N_e, 2/N_e, \ldots, 1\}$，在高深度下产生分辨率偏差。

### 4.2 解决方案

**连续 Beta-diffusion 模型**（对齐 deCODE Cell 论文的 Kimura 估计器）：

$$
p_{\text{child}} \mid p_{\text{mother}} \sim \text{Beta}\bigl(p_m(N_e-1),\; (1-p_m)(N_e-1)\bigr)
$$

**核心改进**：
- $p_{\text{child}}$ 是连续变量，无网格约束
- 在合成数据上验证：Ne_continuous=32 vs Ne_Kimura=31.7（高度一致）

### 4.3 CLI 变更

```bash
--model continuous   # 默认（v1.8.2+）
--model discrete     # 保留旧模型用于对比
```

### 4.4 关键决策

| 决策 | 理由 |
|---|---|
| 连续模型作为默认 | 离散模型在高深度下系统性偏高；连续模型与 Kimura 一致 |
| 保留离散模型选项 | 供用户对比验证；某些低深度场景离散模型仍有诊断价值 |
| 使用 plug-in $p_m = k_m/d_m$ | 避免对母亲等位基因频率积分，简化计算；与 pair 模型保持一致 |

---

## 5. v1.8.3：实值优化与解释文档

### 5.1 实值 $N_e$ 优化

**变更**：将 $N_e$ 搜索从整数网格改为实值 golden-section search。

**结果**：
- JSON 输出 `Ne_MLE` 现为浮点数（如 `Ne_MLE: 31.72`）
- `Result` 结构体 `ne_mle` 字段从 `int` 改为 `double`

### 5.2 README 新增解释段落

**问题**：用户常问"为什么 Ne_MLE 常小于 Ne_Kimura？"

**新增解释**（基于统计学原理）：
1. **Jensen 不等式偏差**：$E[\log f] < \log E[f]$，对数似然凹性导致 MLE 向下偏
2. **全分布建模**：MMLE 使用 Beta 分布建模全部 VAF，Kimura 仅用方差矩
3. **Fisher 最优加权**：MMLE 自动给予低 VAF（高信息量）对更大权重
4. **Plug-in 噪声**：Kimura 使用 $\hat{p}_m$ 引入额外方差

---

## 6. v1.8.4：per-bin drift summary TSV 输出

### 6.1 动机

用户需要一种可视化诊断工具，将观测到的母子漂变按母亲 VAF 分组，并与 Kimura 解析预测 $p(1-p)/N_e$ 进行对比，以直观评估 $N_e$ 拟合质量。

### 6.2 实现

**新增 CLI 标志**：
```bash
--bin-simulation           # Enable per-bin drift summary output
--bin-simulation-bins 20   # bin count (default 20)
```

**TSV 输出格式**：
- `#` 元数据行：拟合的 $N_e$、CI 边界、解析预测公式
- 数据行：每个 maternal-VAF bin 的观测 drift 均值、解析 Kimura 预测 $p(1-p)/N_e$ 曲线值

**绘图工具**：`tools/plot_bottleneck_simulation.py`（matplotlib 双面板布局）

---

## 7. v1.8.5：Ne-profile 诊断

### 7.1 功能

新增 `--ne-profile` 标志，输出 $N_e$ 扫描范围内每个候选值的对数似然，用于：
- 可视化似然曲面形态
- 诊断多峰问题
- 比较 MLE 与 Kimura 估计位置

### 7.2 配套测试

新增 `all-homoplasmic` 回归测试夹具：验证纯质数据下 Ne-profile 行为正确。

---

## 8. v1.8.6：术语修正 MLE → MMLE

### 8.1 动机

统计学严格性：主估计器对潜变量 $p_m$ 执行**解析边缘化**，且对多个 pair 使用**复合似然**，应称为 MMLE（Maximum Marginal Likelihood Estimator）而非 MLE。

### 8.2 变更范围

| 位置 | 变更 |
|---|---|
| JSON 输出字段 | `Max_LogLik` → `Max_Marginal_LogLik` |
| JSON `Estimator` 字段 | 新增 `"MMLE (composite marginal likelihood)"` |
| 内部变量 | `mle_log_lik` → `mmle_log_lik` |
| CLI help 文本 | 所有 "MLE" 替换为 "MMLE" |
| README | 术语解释段落 |

### 8.3 版本决策

**为何选 patch 版本（1.8.6）而非 minor（1.9.0）**：
- JSON 字段重命名是 API 破坏性变更
- 但变更发生在同一开发日内（v1.8.0 → v1.8.6），用户尚未大规模采用
- 遵循"快速修正优于版本膨胀"原则，使用 patch 版本

---

## 9. 关键技术决策汇总

| # | 决策 | 版本 | 理由 |
|---|---|---|---|
| D1 | 离散 Beta-Binomial MLE | v1.8.0 | 初版实现；精确但受网格分辨率限制 |
| D2 | Wonnapinij 采样噪声校正 | v1.8.0 | 使 Kimura 估计在低深度下无偏 |
| D3 | 修复分母交换 bug | v1.8.1 | 实现错误，非设计决策 |
| D4 | 连续 Beta-diffusion 替代离散 MLE | v1.8.2 | 消除网格偏差；与 deCODE 论文对齐 |
| D5 | Plug-in $p_m = k_m/d_m$ | v1.8.2 | 简化计算；与 pair 模型一致 |
| D6 | 实值 golden-section search | v1.8.3 | 消除整数约束；提高精度 |
| D7 | per-bin drift summary TSV | v1.8.4 | 可视化诊断：观测漂变 vs 解析 Kimura 预测 |
| D8 | MLE → MMLE 术语 | v1.8.6 | 统计学准确性；反映边缘化+复合似然 |

---

## 10. 已知限制（遗留至 v1.9.0）

| 问题 | 后续处理 |
|---|---|
| 变异调用器无贝叶斯质量控制 | v1.9.0 新增 `variant-qc` 子命令 |
| Kimura trim 和 outlier 诊断不可见 | v1.9.2 修复 stderr 输出 |
| 仅支持母-子对（2代） | v1.10.0 扩展为三代 G-M-C trio |

---

## 11. 关键文件

| 文件 | 用途 |
|---|---|
| `src/trans_prep.h` / `src/trans_prep.cpp` | trans-prep：VCF+FAM → TSV |
| `src/ne_estimate.h` / `src/ne_estimate.cpp` | ne-estimate：Beta-diffusion MMLE + Kimura cross-check |
| `src/log_factorial.h` / `src/log_factorial.cpp` | 通用对数阶乘工具 |
| `tools/plot_bottleneck_simulation.py` | bin-simulation 绘图 |
| `tests/test_trans_prep.cpp` | trans-prep 测试（20+ 用例） |
| `tests/test_ne_estimate.cpp` | ne-estimate 测试（30+ 用例） |

---

## 12. 下一阶段入口

→ [HANDOFF_v1.9.0-v1.9.2_variant_qc_kimura_ux.md](HANDOFF_v1.9.0-v1.9.2_variant_qc_kimura_ux.md)
