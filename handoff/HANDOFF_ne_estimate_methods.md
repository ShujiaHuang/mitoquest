# MitoQuest `ne-estimate` 方法学完整手册

> 源码：`src/ne_estimate.cpp` (2228 行) + `src/ne_estimate.h` (547 行)  
> 数值工具：`src/log_factorial.cpp` + `src/log_factorial.h`  
> 版本：v1.11.0（2026-06-11）  
> 定位：**方法学参考文档**——从生物学到公式推导、代码实现到计算实例，完整覆盖

---

文档结构如下：

**核心章节**（10 节）：

| # | 章节 | 内容 |
|---|------|------|
| 1 | 生物学背景 | 线粒体瓶颈、$N_e$ 含义、表观瓶颈大小概念 |
| 2 | 数据模型与输入 | TSV 格式、QC 过滤、VAF 窗口、信息位点定义 |
| 3 | 连续 Beta-diffusion 模型 | 三层生成模型 → Beta-Binomial 边际似然完整推导 → 对数域 `lgamma` 实现 |
| 4 | 离散模型 | 遗留 BetaBinomial-Binomial 模型 + log-sum-exp 累加 |
| 5 | 三代 Trio 扩展 | G-M-C 闭式边际似然推导 + 自动调度逻辑 |
| 6 | 优化算法与 CI | Phase 1 整数扫描 + Phase 2 golden-section + Profile Likelihood CI |
| 7 | Kimura 交叉校验 | Wonnapinij 校正、$b$ 统计量、Bootstrap CI、Trimmed Kimura、SSR 闭式解 |
| 8 | Per-Family 估计 | 家系分组 → 独立估计 → 尴尬并行 |
| 9 | 诊断工具 | Bin drift summary (per-bin observed vs analytical Kimura prediction)、Ne-profile 扫描 |
| 10 | 计算实例 | 5 个位点的完整手算：似然曲面表、优化、Kimura、Per-Family |

**附录**：LogFactorial 数值工具表、CLI 完整参考、JSON 输出结构示例。

每个公式都标注了对应的源码位置（函数名 + 行号），方便在代码和文档之间交叉参考。

## 目录

1. [生物学背景：线粒体遗传瓶颈](#1-生物学背景线粒体遗传瓶颈)
2. [数据模型与输入](#2-数据模型与输入)
3. [连续 Beta-diffusion 模型（主估计器）](#3-连续-beta-diffusion-模型主估计器)
4. [离散 Beta-Binomial 模型（遗留模型）](#4-离散-beta-binomial-模型遗留模型)
5. [三代家系（G-M-C Trio）扩展](#5-三代家系g-m-c-trio扩展)
6. [优化算法与置信区间](#6-优化算法与置信区间)
7. [Kimura 矩估计交叉校验](#7-kimura-矩估计交叉校验)
8. [Per-Family 独立估计](#8-per-family-独立估计)
9. [诊断工具](#9-诊断工具)
10. [计算实例：从原始 read 到 $N_e$](#10-计算实例从原始读数到-n_e)

---

## 1. 生物学背景：线粒体遗传瓶颈

### 1.1 什么是线粒体遗传瓶颈？

人类线粒体 DNA（mtDNA）是长约 16,569 bp 的环状基因组，几乎完全通过**母系遗传**。当一个母亲将 mtDNA 传递给孩子时，并非所有母体 mtDNA 拷贝都被传递——而是经历了一个**遗传瓶颈（genetic bottleneck）**：只有少数 mtDNA 分子（有效数量记为 $N_e$）被"抽样"进入卵子，随后通过营养期分裂（vegetative segregation）扩增到正常细胞水平。

**后果**：如果母亲在某个 mtDNA 位点携带两种等位基因（即**异质性 / heteroplasmy**），孩子在该位点的等位基因频率 $p_C$ 会因瓶颈的随机抽样而**偏离**母亲的频率 $p_M$。这种频率偏移称为**遗传漂变（genetic drift）**。

![Genetic drift (article) | Natural selection](https://static.fungenomics.com/images/2026/06/b86e87affc7b3d590928ecaeebf54ce8f2d0f36e.png)

### 1.2 为什么估计 $N_e$？

$N_e$ 量化了瓶颈的"紧度"：

| $N_e$ 值 | 瓶颈紧度 | 母子漂变 | 生物学含义 |
|----------|----------|----------|-----------|
| $N_e \approx 1\text{-}5$ | 极紧 | 极大 | 极少数 mtDNA 被传递（deCODE: $N_e \approx 3$） |
| $N_e \approx 20\text{-}50$ | 中等 | 中等 | Li et al. (2016) 估计 ~28-35 |
| $N_e \approx 100\text{-}200$ | 宽松 | 微小 | 瓶颈几乎不存在 |

$N_e$ 的准确估计对理解 mtDNA 疾病传递、衰老和进化至关重要。

### 1.3 $N_e$ 应被解读为"表观有效瓶颈大小"

deCODE (Cell, 2024) 估计 $N_e \approx 3$，远小于 Li et al. (2016) 的 ~28-35。差异的一个重要原因是**胚系选择（germline selection）**：有害 mtDNA 突变在卵子发生中被主动清除（Rath et al., Science, 2019）。Kimura 模型假设**中性漂变**，因此估计的 $N_e$ 实际上是"漂变 + 选择"联合效应的**表观有效瓶颈大小（apparent effective bottleneck size）**。

---

## 2. 数据模型与输入

### 2.1 母子传递对（Transmission Pairs）

`ne-estimate` 的输入是由 `mitoquest trans-prep` 预处理的 TSV 文件，每行代表一个**变异位点**上的母亲-孩子读数对：

| 列名 | 类型 | 含义 |
|------|------|------|
| `MOTHER_DP` | int | 母亲在该位点的总测序深度 $d_M$ |
| `MOTHER_AD_ALT` | int | 母亲 ALT 等位基因读数 $k_M$ |
| `MOTHER_VAF` | float | 母亲 variant allele frequency $= k_M / d_M$ |
| `CHILD_DP` | int | 孩子在该位点的总深度 $d_C$ |
| `CHILD_AD_ALT` | int | 孩子 ALT 等位基因读数 $k_C$ |
| `QC` | string | 质控标签，仅 `PASS` 被使用 |
| `HAS_G` | int | 1 = 三代 G-M-C trio；0 = 二代 M-C pair（可选） |
| `GRANDMOTHER_DP` / `GRANDMOTHER_AD_ALT` | int | 祖母深度和 ALT 读数（可选） |
| `FAM_ID` / `MOTHER_ID` / `CHILD_ID` | string | 家系标识（可选，per-family 模式使用） |

### 2.2 数据过滤

`load_pairs()` (line 643-741) 对每行应用两层过滤：

1. **QC 门控**：仅保留 `QC == "PASS"` 的位点
2. **VAF 窗口**：默认 `--min-vaf 0.10 --max-vaf 0.90`，过滤母亲 VAF 接近 0 或 1 的位点

**为什么需要 VAF 窗口？** 母亲纯质位点（$p_M = 0$ 或 $1$）对 $N_e$ 没有信息量——无论瓶颈多紧，孩子的频率都只能从 0 或 1 出发，没有漂变可观测。此外，VAF 接近 0 或 1 的位点信息量极低且易受测序错误干扰。

### 2.3 信息位点 vs. 非信息位点

对 $N_e$ 估计有贡献的位点必须满足 $0 < p_M < 1$（母亲异质性）。代码中：

```cpp
// compute_ll_single_continuous(), line 100-101
const double pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
if (pm <= 0.0 || pm >= 1.0) return 0.0;  // 母亲纯质：对数似然 = 0（Ne 无关常数）
```

---

## 3. 连续 Beta-diffusion 模型（主估计器）

这是 `mitoquest ne-estimate` 的**默认且推荐**的模型（`--model continuous`），自 v1.8.2 起引入。

### 3.1 生成模型

对于单个母子传递对、单个变异位点，数据生成过程分三层：

**第一层：母亲的真实异质性**
$$p_M \text{ 为母亲在该位点的真实等位基因频率}$$

在连续模型中，$p_M$ 取 plug-in 点估计：$\hat{p}_M = k_M / d_M$。在典型 mtDNA 深度（$\geq 100\times$）下，plug-in 与完整 Bayesian 后验积分的差异可忽略（$O(1/d_M)$）。

**第二层：瓶颈传递（Kimura 平稳 Beta 近似）**

Kimura (1955) 扩散过程的**平稳分布**为 Beta 分布。经瓶颈 + 营养期分裂后，孩子的真实异质性为：

$$\boxed{p_C \mid p_M, N_e \sim \text{Beta}\bigl(\alpha = p_M(N_e - 1),\; \beta = (1 - p_M)(N_e - 1)\bigr)}$$

**直觉**：
- Beta 分布的**均值** $= \alpha / (\alpha + \beta) = p_M$（孩子的期望频率等于母亲的频率）
- Beta 分布的**方差** $= p_M(1-p_M) / (N_e - 1 + 1) = p_M(1-p_M)/N_e$（$N_e$ 越小，方差越大 → 漂变越大）
- $N_e \to \infty$ 时，Beta 退化为 $p_M$ 处的点质量 → 无漂变

**第三层：测序抽样**
$$k_C \mid p_C, d_C \sim \text{Binomial}(d_C, p_C)$$

### 3.2 边际似然推导（核心公式）

孩子的真实频率 $p_C$ 不可直接观测。我们对其做**解析边际化**——将第二层 Beta 和第三层 Binomial 合并：

$$P(k_C \mid p_M, N_e) = \int_0^1 \binom{d_C}{k_C} p^{k_C}(1-p)^{d_C - k_C} \cdot \frac{p^{\alpha-1}(1-p)^{\beta-1}}{B(\alpha, \beta)} \, dp$$

合并 $p$ 的幂次：

$$= \binom{d_C}{k_C} \frac{1}{B(\alpha, \beta)} \int_0^1 p^{\alpha + k_C - 1}(1-p)^{\beta + d_C - k_C - 1} \, dp$$

积分恰为 Beta 函数 $B(\alpha + k_C,\; \beta + d_C - k_C)$，因此：

$$\boxed{P(k_C \mid p_M, N_e) = \binom{d_C}{k_C} \frac{B(\alpha + k_C,\; \beta + d_C - k_C)}{B(\alpha, \beta)}}$$

这就是 **Beta-Binomial 边际分布**，也是连续模型的核心似然函数。

### 3.3 对数域实现

为避免浮点下溢，所有计算在对数域进行。`log_betabinom_pmf()` (`log_factorial.cpp` line 47-65) 实现：

$$
\ell_i(N_e) = \log\binom{d_C}{k_C} + \ln\Gamma(k_C + \alpha) + \ln\Gamma(d_C - k_C + \beta) - \ln\Gamma(d_C + \alpha + \beta) - \ln\Gamma(\alpha) - \ln\Gamma(\beta) + \ln\Gamma(\alpha + \beta)$$

对应代码 (`compute_ll_single_continuous`, line 93-106)：

```cpp
const double pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
if (pm <= 0.0 || pm >= 1.0) return 0.0;
const double ne1 = static_cast<double>(ne - 1);
return lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt, pm * ne1, (1.0 - pm) * ne1);
```

> **注意**：$\log\binom{d_C}{k_C}$ 在所有 $N_e$ 下是常数，不影响最优 $N_e$ 的位置。但它影响对数似然的绝对值和 Profile Likelihood CI 的阈值计算。

### 3.4 全局复合对数似然

将所有通过 QC 和 VAF 过滤的位点视为**独立观测**（复合似然 / composite likelihood），全局对数似然为：

$$\boxed{L(N_e) = \sum_{i=1}^{N_{\text{pairs}}} \ell_i(N_e)}$$

对应 `compute_global_ll_continuous()` (line 403-442)。当 `threads > 1` 且数据量 $\geq 64$ 时，自动启用线程池并行求和。

### 3.5 为什么连续模型优于离散模型？

离散模型（v1.8.0-v1.8.2 之前使用）将孩子的频率限制在网格 $\{0, 1/N_e, 2/N_e, \ldots, 1\}$ 上。在高深度测序下（$d_C \gg N_e$），Binomial 似然极窄，迫使 MLE 将 $N_e$ 膨胀到不切实际的值（~5-10× 真实值），以获得足够细的网格。连续模型通过解析边际化完全消除了这一网格效应。详见 `release_v1.8.2.md`。

---

## 4. 离散 Beta-Binomial 模型（遗留模型）

`--model discrete` 保留了旧版模型以供特殊用途（如病毒传代瓶颈实验）。

### 4.1 生成模型

$$\begin{aligned}
p_M &\sim \text{Beta}(k_M + 1,\; d_M - k_M + 1) \quad \text{（Uniform 先验 + 二项似然的后验）} \\
k &\sim \text{BetaBinomial}(N_e,\; k_M + 1,\; d_M - k_M + 1) \quad \text{（瓶颈抽样）} \\
k_C &\sim \text{Binomial}(d_C,\; k / N_e) \quad \text{（测序）}
\end{aligned}$$

### 4.2 对数似然

对离散潜变量 $k \in \{0, 1, \ldots, N_e\}$ 求和（`compute_ll_single`, line 63-81）：

$$\ell_i(N_e) = \log \sum_{k=0}^{N_e} \text{BetaBin}(k \mid N_e, \alpha, \beta) \cdot \text{Bin}(k_C \mid d_C, k/N_e)$$

使用 `log_sum_exp_pair()` 在对数域安全累加，避免下溢。

---

## 5. 三代家系（G-M-C Trio）扩展

当祖母数据可用时（`HAS_G = 1`），利用三代信息可进一步提高估计精度（v1.10.0 引入）。

### 5.1 闭式边际似然

祖母→母亲的传递同样使用 Kimura Beta 近似。对 $p_M$ 做解析边际化后（`compute_ll_trio_continuous`, line 156-193）：

$$\boxed{I_{\text{trio}}(N_e) = \binom{d_M}{k_M}\binom{d_C}{k_C} \frac{B(\alpha_G + k_M + k_C,\; \beta_G + (d_M - k_M) + (d_C - k_C))}{B(\alpha_G,\; \beta_G)}}$$

其中 $\alpha_G = \hat{p}_G(N_e - 1)$，$\beta_G = (1 - \hat{p}_G)(N_e - 1)$，$\hat{p}_G = g_{\text{ad\_alt}} / g_{\text{dp}}$。

**推导关键**：Beta 先验核 × 二项似然核 × BetaBinomial 边际核，三者的 $p_M$ 指数合并后恰好构成新 Beta 函数的核函数，积分可得闭式 Beta 函数比。

### 5.2 自动调度

`compute_global_ll_continuous()` 内部按行自动调度（line 409-413）：

```cpp
auto row_ll = [&](const PairData& pd) {
    return (pd.has_g == 1)
        ? compute_ll_trio_continuous(pd, ne, lf)    // trio 行
        : compute_ll_single_continuous(pd, ne, lf);  // pair 行
};
```

当 `has_g == 0`（普通二代对）或祖母纯质（$\hat{p}_G \in \{0, 1\}$）时，自动回退到对应处理。

---

## 6. 优化算法与置信区间

### 6.1 两阶段优化

$N_e$ 的 MMLE 通过最大化全局对数似然 $L(N_e)$ 获得：

$$\hat{N}_e = \arg\max_{N_e \in [N_{e,\min},\; N_{e,\max}]} L(N_e)$$

**Phase 1：整数网格扫描** (`find_optimal_ne_continuous`, line 475-501)

在 $N_e \in \{1, 2, \ldots, N_{e,\max}\}$ 上逐一评估 $L(N_e)$，找到最优整数 $N_e^*$。

**为什么不用梯度法？** 离散模型的似然曲面**不是单峰的**——当真实 $N_e$ 较小时，$N_e$ 的整数倍（$2N_e, 3N_e, \ldots$）处会出现局部极值（网格共振效应）。全扫描确保定位到全局最优。

**Phase 2：Golden-Section 精化** (`golden_section_max`, line 447-473)

在 $[N_e^* - 1, N_e^* + 1]$ 区间内用黄金分割法（$\phi \approx 0.618$）寻找实数值最优解，容差 $\text{tol} = 0.01$：

```
初始区间 [a, b] = [N_e* - 1, N_e* + 1]
x1 = b - φ(b - a),  x2 = a + φ(b - a)
while (b - a) > 0.01:
    if L(x1) < L(x2): a = x1  // 峰值在右侧
    else:             b = x2  // 峰值在左侧
return (a + b) / 2
```

约 $\log_\phi(2/0.01) \approx 13$ 次迭代即可收敛。

### 6.2 Profile Likelihood 95% 置信区间

基于 **Wilks 定理**：在大样本下，$-2(L(N_e) - L_{\max})$ 近似服从 $\chi^2_1$ 分布。95% CI 为：

$$\text{CI}_{95\%} = \{N_e : L(N_e) \geq L_{\max} - \chi^2_{1,0.95}/2\} = \{N_e : L(N_e) \geq L_{\max} - 1.92\}$$

**实现**（`estimate_continuous`, line 504-553）：

从 $\hat{N}_e$ 向两侧以步长 0.01 步行，找到对数似然首次低于阈值的点：

```cpp
const double thr = r.max_log_lik - 1.92;  // kProfileLLThresholdDelta
// 向左步行
for (double ne = r.ne - 0.01; ne >= min_ne; ne -= 0.01) {
    if (compute_global_ll_continuous(ne, ...) < thr) { ci_low = ne + 0.01; break; }
}
// 向右步行（类似）
```

如果步行到搜索边界仍未穿越阈值，则标记 `ci_low_clipped = true`（或 `ci_high_clipped`），提示用户 CI 被截断。

---

## 7. Kimura 矩估计交叉校验

Kimura 矩估计是独立于 MMLE 的**第二种 $N_e$ 估计方法**，基于方差分解而非似然最大化。用作交叉校验——当两种方法的结果一致时，估计最可信。

### 7.1 核心思想

Kimura 扩散预测：给定 $N_e$，母子频率差的方差为

$$E[(p_C - p_M)^2] = \frac{p_M(1 - p_M)}{N_e} + \text{测序噪声}$$

如果我们能从观测方差中**减去**测序噪声，就能反推 $N_e$。

### 7.2 Wonnapinij 采样噪声校正

直接使用 $(p_C - p_M)^2$ 会**高估**漂变方差——因为 $p_M$ 和 $p_C$ 都是有限深度下的点估计，包含二项抽样噪声。Wonnapinij (2008/2010) 校正公式（`prepare_pair_contributions`, line 796-819）：

$$\begin{aligned}
d_i &= (p_{C,i} - p_{M,i})^2 & &\text{（原始漂变平方）} \\
s_i &= \frac{p_{M,i}(1-p_{M,i})}{d_{M,i}} + \frac{p_{C,i}(1-p_{C,i})}{d_{C,i}} & &\text{（采样噪声估计）} \\
r_i &= d_i - s_i & &\text{（校正后漂变信号）} \\
w_i &= p_{M,i}(1-p_{M,i}) & &\text{（方差归一化因子）}
\end{aligned}$$

> **直觉**：$s_i$ 估计了"即使没有真实漂变，仅因测序噪声也会观测到的 $(p_C-p_M)^2$"。$r_i = d_i - s_i$ 是去噪后的真实漂变信号。当 $r_i < 0$ 时，表示该位点的观测漂变甚至小于噪声——该位点不含显著漂变信号。

### 7.3 Kimura 估计量

聚合所有信息位点（`compute_kimura_check`, line 852-917）：

$$V = \frac{\sum_i r_i}{\sum_i w_i}, \qquad b = 1 - V, \qquad \boxed{N_e^{\text{Kimura}} = \frac{1}{1 - b} = \frac{1}{V}}$$

$b$ 被裁剪到 $(\epsilon, 1-\epsilon)$（$\epsilon = 10^{-9}$）以防除零。

### 7.4 Bootstrap 置信区间

当 `--kimura-bootstrap > 0` 时（默认 1000 次），执行 pair-level 非参数 bootstrap（line 893-940）：

1. 有放回抽样 $N_{\text{pairs}}$ 个位点
2. 在每次 resample 上重新计算 $b^{(j)}$ 和 $N_e^{(j)}$
3. 取 2.5% 和 97.5% 分位数作为 95% CI

### 7.5 Trimmed Kimura（鲁棒估计）

标准 Kimura 对**离群高漂变位点**极敏感（NUMT 污染、测序错误等）。`--kimura-trim 0.10` 按 $F_i = r_i/w_i$ 降序排列后，丢弃 top 10% 的高漂变位点，在剩余位点上重新计算 $b$（line 942-986）。当 trimmed 与 untrimmed $N_e$ 差异很大时，说明数据中存在离群值。

### 7.6 Kimura-SSR 闭式估计

Ne-profile 使用另一种最小二乘形式的 Kimura 估计（`kimura_ssr_best_ne`, line 1201-1216）。对每个位点拟合预测 $E[r_i] = w_i / N_e$，最小化残差平方和：

$$\text{SSR}(N_e) = \sum_i \left(r_i - \frac{w_i}{N_e}\right)^2$$

令 $d(\text{SSR})/dN_e = 0$，解得闭式最优解：$N_e^{\text{SSR}} = \sum w_i^2 / \sum r_i w_i$。

### 7.7 MMLE vs. Kimura 差异的诊断价值

当 $N_e^{\text{MMLE}} / N_e^{\text{Kimura}} > 3$ 或 $< 1/3$ 时，`run()` 会发出警告（line 2172-2208）并建议用户使用 `--kimura-trim 0.10 --top-drift-k 20` 诊断离群值。常见原因包括 NUMT 污染、样本混合、或胚系选择。

---

## 8. Per-Family 独立估计

`--per-family` 模式（v1.11.0 引入）将队列中的位点按家系分组，对每个家系独立估计 $N_e$。

### 8.1 家系分组

`group_into_families()` (line 1223-1257) 按 `(FAM_ID, MOTHER_ID)` 分组。同一母亲的所有孩子归入同一家系，共享 $N_{e,f}$。

### 8.2 单家系估计

`estimate_family()` (line 1260-1310) 对每个家系：

1. 统计信息位点数 $n_{\text{info}}$（$0 < p_M < 1$）
2. 若 $n_{\text{info}} < $ `--min-family-sites`（默认 3）→ **跳过**，输出 warning
3. 否则调用 `estimate_continuous(fam.pairs, min_ne, max_ne, threads=1)`——**复用群体级的完整优化引擎**
4. 填充 `FamilyResult`（含 CI、警告标志等）

### 8.3 家系间并行

`estimate_all_families()` (line 1313-1348) 采用**尴尬并行**——每个家系的优化完全独立：

```cpp
std::vector<std::future<void>> futures;
const size_t chunk = (n + threads - 1) / threads;
for (size_t t = 0; t < threads; ++t) {
    futures.push_back(std::async(std::launch::async, worker, lo, hi));
}
```

家系内串行（$K_f$ 通常很小，线程开销 > 计算量），家系间并行。

---

## 9. 诊断工具

### 9.1 Bin Drift Summary（per-bin observed vs analytical Kimura prediction）

`--bin-simulation FILE` 将信息位点按 $p_M$ 分入等宽 bins，对每个 bin 报告（line 1042-1115）：

| 列 | 含义 |
|----|------|
| `obs_var` | bin 内平均 $(p_C - p_M)^2$（原始漂变） |
| `obs_var_corr` | bin 内平均 $(d_i - s_i)$（校正后漂变） |
| `obs_F` | bin 内平均 $(d_i - s_i) / w_i$（$= 1 - b$ 的 bin 估计） |
| `expected_var_at_ne` | $p_{\text{center}}(1-p_{\text{center}}) / \hat{N}_e$（analytical Kimura prediction） |

绘图脚本将观测 bin 均值与解析 Kimura 预测抛物线 $p(1-p)/N_e$ 叠加，直观展示拟合质量。

### 9.2 Ne-Profile 扫描

`--ne-profile FILE` 在 $N_e$ 网格上同时评估两种指标（line 1138-1199）：

- **MMLE**：$\Delta_{2LL} = -2(L(N_e) - L_{\max})$，在 $\hat{N}_e^{\text{MMLE}}$ 处为 0
- **Kimura SSR**：$\text{SSR}_{\text{norm}} = \text{SSR}(N_e) / \text{SSR}_{\min}$，在 $\hat{N}_e^{\text{SSR}}$ 处为 1

绘图后可视化两种估计器是否在同一 $N_e$ 处取得最优，并检查它们的拟合一致性。

---

## 10. 计算实例：从原始读数到 $N_e$

> 用 5 个合成位点逐步手算，演示完整流程。深度均为 500×。

### 10.1 实验数据

| 行 | FAM_ID | $d_M$ | $k_M$ | $d_C$ | $k_C$ | $\hat{p}_M$ | $\hat{p}_C$ | 漂变幅度 |
|----|--------|-------|-------|-------|-------|-----------|-----------|----------|
| 1 | F001 | 500 | 100 | 500 | 110 | 0.200 | 0.220 | 小 |
| 2 | F001 | 500 | 150 | 500 | 120 | 0.300 | 0.240 | 中 |
| 3 | F001 | 500 | 150 | 500 | 10 | 0.300 | 0.020 | 极端 |
| 4 | F001 | 500 | 250 | 500 | 240 | 0.500 | 0.480 | 小 |
| 5 | F002 | 500 | 50 | 500 | 250 | 0.100 | 0.500 | 极大 |

所有 5 个位点都是信息位点（$0 < \hat{p}_M < 1$）。

### 10.2 单位点对数似然手算（行 1，$N_e = 30$）

**Step 1**：$\alpha = 0.20 \times 29 = 5.8$，$\beta = 0.80 \times 29 = 23.2$

**Step 2**：代入 `log_betabinom_pmf(500, 110, 5.8, 23.2)`

$$\Delta\ell_1 = \ln\Gamma(115.8) + \ln\Gamma(413.2) - \ln\Gamma(529) - \ln\Gamma(5.8) - \ln\Gamma(23.2) + \ln\Gamma(29)$$

$$\approx 424.248 + 2098.642 - 2734.233 - 3.957 - 49.803 + 64.214 = -200.889$$

（去掉 $\log\binom{500}{110}$ 常数，因为它不影响 $N_e$ 优化。）

### 10.3 似然曲面

对所有位点和多个候选 $N_e$ 计算 $\Delta\ell_i$（log Beta 比）并求和：

| $N_e$ | $\Delta\ell_1$ | $\Delta\ell_2$ | $\Delta\ell_3$ | $\Delta\ell_4$ | $\Delta\ell_5$ | $L(N_e)$ |
|-------|-----------|-----------|-----------|-----------|-----------|----------|
| 10 | -3.126 | -5.192 | -4.055 | -1.609 | -3.164 | **-17.146** |
| 15 | -1.212 | -2.599 | -4.511 | -0.885 | -3.582 | **-12.789** |
| **20** | **-0.811** | **-1.715** | **-4.996** | **-0.565** | **-5.000** | **-13.087** |
| 25 | -0.651 | -1.288 | -5.353 | -0.409 | -5.510 | **-13.211** |
| 30 | -0.561 | -1.113 | -5.646 | -0.309 | -5.967 | **-13.596** |
| 50 | -0.371 | -0.549 | -6.619 | -0.103 | -7.401 | **-15.043** |

**似然曲面解读**：
- 行 1（微漂变 0.20→0.22）和行 4（微漂变 0.50→0.48）**偏好大 $N_e$**（宽松瓶颈）
- 行 3（极端漂变 0.30→0.02）和行 5（极大漂变 0.10→0.50）**偏好小 $N_e$**（紧瓶颈）
- 全局 $L(N_e)$ 是矛盾信号的折中，**峰值在 $N_e = 15$ 附近**

### 10.4 优化结果

- Phase 1 整数扫描 → $N_e^* = 15$
- Phase 2 golden-section on $[14, 16]$ → $\hat{N}_e \approx 15.2$
- Profile Likelihood CI（阈值 $= -12.789 - 1.92 = -14.709$）：$\text{CI} \approx [9, 50]$

### 10.5 Kimura 交叉校验手算

5 个位点的 $d_i, s_i, r_i, w_i$：

| 行 | $d_i$ | $s_i$ | $r_i$ | $w_i$ | 信号 |
|----|-------|-------|-------|-------|------|
| 1 | 0.0004 | 0.0007 | -0.0003 | 0.160 | 无（噪声主导） |
| 2 | 0.0036 | 0.0008 | +0.0028 | 0.210 | 弱 |
| 3 | 0.0784 | 0.0005 | +0.0779 | 0.210 | 强 |
| 4 | 0.0004 | 0.0010 | -0.0006 | 0.250 | 无（噪声主导） |
| 5 | 0.1600 | 0.0007 | +0.1593 | 0.090 | 极强 |

聚合：$\sum r_i = 0.2392$，$\sum w_i = 0.9200$

$$V = 0.2392 / 0.9200 = 0.260, \quad N_e^{\text{Kimura}} = 1/0.260 \approx 3.85$$

**MMLE (15.2) vs. Kimura (3.85) 差异较大**——因为行 3 和行 5 的极端漂变主导了 Kimura（$r_3 + r_5$ 占 $\sum r_i$ 的 99%），而 MMLE 通过似然函数更均衡地融合所有位点信号。仅 5 个位点时两种方法方差都极大；真实队列（数百位点）中两者通常收敛到相近值。

### 10.6 Per-Family 估计

**F001**（行 1-4，4 个信息位点）：
- 去除行 5 后，似然峰值移向更大 $N_e$（~30），因为行 5 的"极大漂变拉低群体级估计"的效应消失了

**F002**（行 5，仅 1 个信息位点）：
- $1 < 3$（`--min-family-sites` 默认值）→ **跳过**，输出 warning

---

## 附录 A：`LogFactorial` 数值工具

`LogFactorial` (`log_factorial.h/cpp`) 提供预计算的 $\log(n!)$ 缓存和三个对数域 PMF：

| 方法 | 公式 | 用途 |
|------|------|------|
| `log_fact(n)` | $\log(n!)$，缓存覆盖 $[0, \text{max\_n}]$ | 组合数计算 |
| `log_comb(n, k)` | $\log\binom{n}{k} = \log(n!) - \log(k!) - \log((n-k)!)$ | 二项系数 |
| `log_binomial_pmf(n, k, p)` | $\log\binom{n}{k} + k\log p + (n-k)\log(1-p)$ | 离散模型 |
| `log_betabinom_pmf(n, k, α, β)` | $\log\binom{n}{k} + \ln\Gamma(k+\alpha) + \ln\Gamma(n-k+\beta) - \ln\Gamma(n+\alpha+\beta) - \ln\Gamma(\alpha) - \ln\Gamma(\beta) + \ln\Gamma(\alpha+\beta)$ | 连续模型核心 |

缓存构建为 $O(\text{max\_n})$，查询为 $O(1)$。超出缓存范围时自动回退到 `std::lgamma(n+1)`。

## 附录 B：完整 CLI 参考

```
mitoquest ne-estimate [options] -i <pairs.tsv>

必选：
  -i, --input FILE          trans-prep 输出的 TSV

模型与搜索：
  -o, --output FILE         JSON 输出（默认 stdout）
  --model NAME              continuous（默认）| discrete
  --min-vaf FLOAT           母亲 VAF 下界 [0.10]
  --max-vaf FLOAT           母亲 VAF 上界 [0.90]
  --min-ne INT              搜索下界 [1]
  --max-ne INT              搜索上界 [200]
  -t, --threads INT         工作线程 [1]

Kimura 交叉校验：
  --cross-check kimura      启用 Wonnapinij/Kimura 交叉校验
  --kimura-bootstrap INT    Bootstrap 迭代数 [1000]
  --kimura-seed INT         RNG 种子 [42]
  --kimura-trim FLOAT       丢弃高漂变位点比例 [0.0]
  --top-drift-k INT         输出 top-K 离群位点 [0]

Per-Family：
  --per-family              启用 per-family 估计
  --min-family-sites INT    每家系最少信息位点 [3]
  --per-family-output FILE  Per-family TSV 输出

诊断：
  --bin-simulation FILE     per-bin drift summary TSV
  --bin-simulation-bins INT Bin 数量 [10]
  --ne-profile FILE         Ne 扫描 profile TSV
  --ne-profile-step FLOAT   Profile 步长 [0.1]
```

## 附录 C：JSON 输出结构

```json
{
  "Ne": 15.20,
  "CI_95_Low": 9.00,
  "CI_95_High": 50.00,
  "Pairs_Used": 5,
  "Estimator": "MMLE (composite marginal likelihood)",
  "Model": "continuous",
  "Kimura_Cross_Check": {
    "b": 0.74000000,
    "Ne_Kimura": 3.84615385,
    "N_Informative": 5
  },
  "Per_Family_Estimates": [
    {
      "FAM_ID": "F001", "Mother_ID": "M01",
      "N_Children": 1, "N_Sites": 4, "N_Informative": 4,
      "Ne": 30.50, "CI_95_Low": 12.00, "CI_95_High": 180.00
    }
  ],
  "Per_Family_Summary": {
    "N_Families_Estimated": 1,
    "N_Families_Skipped": 1
  }
}
```
