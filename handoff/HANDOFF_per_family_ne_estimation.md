# Handoff — Per-Family mtDNA Bottleneck $N_e$ Estimation

> 生成时间：2026-06-11  
> 主题：在 `ne-estimate` 中增加**以家系为单位**的线粒体遗传瓶颈 $N_e$ 估计  
> 前置阶段：[HANDOFF_v1.10.0_20250528.md](HANDOFF_v1.10.0_20250528.md)  
> 方法学定位：本文档为**方法学设计文档**，包含完整数学推导、算法设计、文献审查和实现规划

---

## 1. 问题背景与动机

### 1.1 现状

当前 `mitoquest ne-estimate` 将输入 TSV 中的所有母子传递对视为来自**同一群体**，拟合一个群体级别的共享 $N_e$（composite marginal likelihood）。这一估计反映的是**队列平均**瓶颈大小。

然而，多项研究（Li et al., 2016 Genome Res.; Rebolledo-Jaramillo et al., 2019 PNAS; Arnadóttir et al. / deCODE, Cell 2024）表明：

1. **瓶颈大小在家系间存在显著变异**：不同母-子对的 $N_e$ 可能相差数个数量级（~10 至 ~500）
2. **$N_e$ 与母系年龄相关**：高龄母亲传递更少的 mtDNA 有效拷贝（Li et al., 2016）
3. **群体均值掩盖个体差异**：队列级 $N_e$ 无法反映特定家系的瓶颈特征
4. **真实瓶颈可能远小于此前估计**：deCODE (Cell, 2024) 使用 116,663 对母子对估计 $N_e \approx 3$，发现胚系选择证据（详见 §2.3.1）

### 1.2 需求

在 `ne-estimate` 中新增功能：**以家系为单位**，利用每个家系内部所有异质性变异位点的传递信息，独立估计该家系自身的 $N_e$。

### 1.3 可行性论证

单个母-子对虽然仅经历**一次瓶颈事件**，但 mtDNA 基因组包含 ~16,569 bp，在高深度测序下（典型 $\geq 100\times$），母亲在多个位点携带可检测的异质性（heteroplasmy）。**每个异质性位点是同一瓶颈过程的一次独立实现**——因为不同 mtDNA 变异在传递时独立经历漂变。因此：

> **单一家系的 $K$ 个异质性位点 = 从同一 Beta-diffusion 过程中获得 $K$ 个独立样本**

这使得对单个家系进行 $N_e$ 估计在统计上是可行的（Li et al., 2016; Arnadóttir et al. / deCODE, Cell 2024）。

> **重要说明**：deCODE (Cell, 2024) 估计 $N_e \approx 3$，远小于 Li et al. (2016) 的 ~28-35。我们估计的 $N_e$ 应解读为「表观有效瓶颈大小」（apparent effective bottleneck size），即中性漂变假设下与实际传递方差匹配的等效值（详见 §2.3.1）。

---

## 2. 已有文档审查与勘误

在开始方法推导之前，我对已有 handoff 文档进行了逐份审查。以下为发现的问题：

### 勘误 1：v1.10.0 HANDOFF — D3 措辞歧义

**原文**（HANDOFF_v1.10.0, 决策 D3）：

> "混合队列估计单一共享 $N_e$（二代对 + 三代 trio 共用同一参数）"

**问题**：措辞暗示二代对和三代 trio 仅使用一次 bottleneck 事件。实际上，三代 trio 模型使用**两次** bottleneck 事件（G→M 和 M→C），每次均使用同一 $N_e$ 参数。正确的表述应为：

> "混合队列估计单一共享 $N_e$ 参数；二代对使用该参数建模一次 M→C 传递，三代 trio 使用该参数建模两次独立传递（G→M 和 M→C），均假设同一有效瓶颈大小。"

**严重程度**：措辞歧义，不影响代码正确性。

### 勘误 2：v1.8.0 HANDOFF — v1.8.0 离散模型公式不完整

**原文**（HANDOFF_v1.8.0, §2.4）：

> "离散 Beta-Binomial MLE：对每个 mother-child pair，在离散网格 $k \in \{0, 1, \ldots, N_e\}$ 上计算条件似然"

**问题**：未明确指出 $p_m$ 使用的是 Beta 后验（Uniform 先验 + 二项似然 → Beta(m_alt+1, m_ref+1)），而非简单的点估计。这一区别在理解离散模型与连续模型的差异时至关重要——离散模型对 $p_m$ 进行了完整的 Bayesian 后验积分，而连续模型（v1.8.2+）使用 plug-in $\hat{p}_m = m_{\text{alt}} / m_{\text{dp}}$。

**严重程度**：低，handoff 是概览文档。

### 勘误 3：v1.8.2 HANDOFF — "Beta-diffusion" 术语精度

**原文**（多处）：

> "连续 Beta-diffusion 模型"

**问题**：技术上讲，Kimura (1955) 扩散过程的 **瞬时** 转移密度并非 Beta 分布；Beta 分布是 Wright-Fisher 扩散的**平稳分布**（stationary distribution）。模型实际使用的是平稳 Beta 分布来近似经 bottleneck + 营养期分裂后的异质性分布。"Beta-diffusion" 作为简称是可接受的，但精确表述应为 "Kimura 平稳 Beta 近似"（stationary Beta approximation to the Kimura diffusion）。

**严重程度**：低，术语精度问题，不影响数学正确性。

### 勘误 4：v1.10.0 HANDOFF — 闭式 Trio 公式推导中的隐含假设

**原文**（HANDOFF_v1.10.0, §3）：闭式 G-M-C Trio 边际似然公式。

**问题**：公式推导中有一个**隐含的代数巧合**需要指出。给定 $p_M$，母亲读数的二项似然 $\text{Bin}(k_M | d_M, p_M)$ 与孩子的 BetaBinomial 似然 $\text{BetaBin}(k_C | d_C, p_M(N_e-1), (1-p_M)(N_e-1))$ 在 $p_M$ 积分时恰好能合并为 Beta 函数的闭式比。这依赖于 BetaBinomial 的核函数 $p_M^{k_C}(1-p_M)^{d_C-k_C} \cdot p_M^{\alpha_G-1}(1-p_M)^{\beta_G-1}$ 能被吸收为单一 Beta 积分——这一性质**仅对 Beta-Binomial 共轭族成立**。如果未来考虑非共轭模型（如含选择的传递模型），此闭式解将不再适用，需要回退到数值积分（Gauss-Legendre quadrature）。

**严重程度**：无，这是设计说明，非错误。

### 勘误 5：v1.9.2 HANDOFF — Bug 修复位置行号过时

**原文**（HANDOFF_v1.9.2, §4.2）：

> "BUG (line ~1322)"

**问题**：随着后续版本（v1.10.0 新增 trio 代码约 200 行）的添加，实际 bug 位置行号已偏移。当前代码中相关逻辑位于 `ne_estimate.cpp` 的 `run()` 方法（约 line ~1830 附近）。行号引用应标注为 "v1.9.2 时代 line ~1322"。

**严重程度**：低，仅影响可追溯性。

### 2.1 `gem_method.md` 方法审查与勘误

在融入新方法之前，我对 `gem_method.md` 中提出的两个方案进行了逐条数学审查。

#### 方案 1 审查结论：矩估计基本正确，但**缺少采样噪声校正**

**方案 1 核心公式**：$\hat{N}_e = 1/\bar{b}$，其中 $b_i = (p_{Oi}-p_{Mi})^2/[p_{Mi}(1-p_{Mi})]$。

**问题 1 — 缺失 Wonnapinij 采样噪声校正**：

$b_i$ 的期望值并非 $1/N_e$，而是：

$$E[b_i] = \frac{1}{N_e} + \underbrace{\frac{1}{d_{M,i}} + \frac{p_{Oi}(1-p_{Oi})}{p_{Mi}(1-p_{Mi}) \cdot d_{C,i}}}_{\text{测序二项抽样噪声}}$$

在有限深度下（如 $d = 100\times$），采样噪声项的量级为 $\sim 0.01$，而真实漂变信号在 $N_e = 30$ 时仅为 $1/30 \approx 0.033$。**未校正的矩估计会系统性高估 $1/N_e$，从而低估 $N_e$**。这正是 MitoQuest v1.8.1 修复的 Kimura 采样噪声校正 bug（见 [HANDOFF_v1.8.0-v1.8.6](HANDOFF_v1.8.0-v1.8.6_ne_estimation_pipeline.md) §3）。

> **修正**：必须使用 Wonnapinij (2010) 校正公式：
> $$r_i = d_i - s_i, \quad d_i = (p_{Oi}-p_{Mi})^2, \quad s_i = \frac{p_{Mi}(1-p_{Mi})}{d_{Mi}} + \frac{p_{Oi}(1-p_{Oi})}{d_{Ci}}$$
> $$\hat{N}_e^{\text{Kimura}} = \frac{\sum w_i}{\sum r_i}, \quad w_i = p_{Mi}(1-p_{Mi})$$

**问题 2 — 计算实例偏差**：

方案 1 计算实例中数值运算本身正确（$b_1 = 0.015625$, $\bar{b} = 0.0362$, $\hat{N}_e = 27.62$），但由于未做采样校正，结果偏高。在 500× 深度下，校正后 $\hat{N}_e \approx 31$（详见 §3.8 修正后计算实例）。

**问题 3 — MLE 描述不精确**：

方案 1 称 MLE 使用 "Kimura 分布" $f(p_O|p_M,N_e)$，并建议使用 L-BFGS 优化。事实上：
- MitoQuest 使用的是 **Beta-Binomial 边际似然**（已对 $p_C$ 积分），而非连续 $p_C$ 的 Kimura PDF
- 优化使用 **golden-section search**（非 L-BFGS），因为离散模型的似然曲面可能非单峰

#### 方案 2 审查结论：思路可取，但不是直接的 $N_e$ 估计

**方案 2 的核心价值**：两步残差法（vGWAS）是一个**下游应用框架**——将漂变偏离群体基准的残差定义为连续表型，用于核基因组 GWAS 关联分析。这与 Rath et al. (Science, 2019) 的母-子传递研究和 Laricchia et al. (Nature, 2023) 的 GWAS 路线一致。

**但方案 2 存在以下问题**：

| 问题 | 详情 |
|---|---|
| 不产生绝对 $N_e$ | 残差 $\epsilon$ 是相对度量，无法回答 "该家系的瓶颈是多少" |
| 缺失采样校正 | $D_j = \text{mean}[(p_O-p_M)^2/w]$ 未减去 $s_i$，同方案 1 |
| 分母硬平滑 | `denom = np.where(denom == 0, 0.001*0.999, denom)` 是 ad hoc 处理，不如 VAF 窗口过滤严谨 |
| 声称离散模型 "死路" 过于绝对 | 稀疏数据问题是 per-family 估计的通用挑战（连续模型同样受困于 $K < 3$），非离散/连续之别 |

**结论**：方案 2 作为 GWAS 下游分析框架**可取但需修正**；作为 per-family $N_e$ 估计方法**不充分**。

#### 综合判断

| 方案 | 适合 per-family $N_e$ 估计？ | 适合 GWAS 表型？ | 融入决策 |
|---|---|---|---|
| 方案 1（矩估计） | 需加入采样校正后可用 | 间接 | **融入**为校正后的 Kimura 交叉校验 |
| 方案 2（两步残差） | 否 | 是 | **融入**为下游 GWAS 应用模块 |
| Beta-Binomial MMLE | 最佳（已有推导） | 可转换 | 保持为主估计器 |

### 2.2 三篇参考文献审查与勘误

`gem_method.md` 方案 2 引用了三篇文献作为方法学背书。我逐篇阅读原文后发现**引用信息全部有误**，且第三篇与 mtDNA 瓶颈估计**方法学关联不成立**。

#### 勘误 6：deCODE 引用信息不完整

**gem_method.md 引用**：deCODE genetics (Cell, 2024), *The rate and nature of mitochondrial DNA mutations in human pedigrees.*

**实际论文**：Arnadóttir et al. (Cell, 2024), DOI: S0092-8674(24)00531-2。完整标题相同，但 gem_method.md 未给出作者和 DOI。

**关键发现与我们的方法的关系**（详见 §2.3）：
- 估计人类 mtDNA 瓶颈 $N_e \approx 3$（远小于此前 Li et al. 2016 的 ~28-35 估计）
- 发现强烈的**胚系选择（germline selection）**证据：有害 mtDNA 突变在卵子发生中被主动清除
- 2,548 个母系家系、116,663 对母子传递对

#### 勘误 7：Chinnery 引用信息严重错误——实际指向两篇不同论文

**gem_method.md 引用**：Patrick Chinnery (Nature Genetics, 2021), *Nuclear genetic control of mitochondrial heteroplasmy and mtDNA copy number in human lineages.*

**实际匹配**：此引用标题混合了两篇 Chinnery 团队相关论文，期刊、年份、作者均不匹配：

| 方面 | gem_method.md 引用 | 实际论文 A | 实际论文 B |
|---|---|---|---|
| 标题 | "Nuclear genetic control of..." | **Germline selection shapes human mitochondrial DNA diversity** | Nuclear genetic control of mtDNA copy number and heteroplasmy in humans |
| 作者 | Patrick Chinnery | **Rath, Gupta, Todres, ... Chinnery** | Laricchia, Gupta, Karczewski et al. |
| 期刊/年份 | Nature Genetics, 2021 | **Science, 2019** | Nature, 2023 |
| DOI | — | **10.1126/science.aau6520** | 10.1038/s41586-023-06426-5 |

**论文 A（Rath et al., Science, 2019）** 与我们的主题**最直接相关**：
- 分析了 12,975 个全基因组序列，其中 1,526 对母子传递对
- 发现 45.1% 的个体携带 mtDNA 异质性
- **首次大规模证明胚系选择（germline selection）塑造 mtDNA 多样性**：
  - 母-子传递中，已知变异比新发变异更易被传递
  - 特定基因组区域存在双向选择（对某些变异正向选择，对某些负向选择）
  - 新发异质性更倾向匹配核基因组成（而非 mtDNA 单倍群），证明**核基因控制胚系选择**
- 在 40,325 个个体中验证了发现

> **这篇论文是 gem_method.md 最可能想引用的文献**：Patrick Chinnery 是作者，主题正是母-子 mtDNA 传递中的选择与瓶颈——与我们 per-family $N_e$ 估计的核心问题完全一致。

**论文 B（Laricchia et al., Nature, 2023）** 与 GWAS 方法学相关：
- 274,832 个个体（UKB + AoU），使用异质性水平做 GWAS
- 发现 42 个核基因关联（SSBP1, TFAM, DGUOK, POLG2, PNP, LONP1, POLRMT 等）
- 使用 Inverse Rank Normalization (INT) 正态化
- 协变量校正使用残差法

#### 勘误 8：Campbell 引用信息错误且方法学关联不成立

**gem_method.md 引用**：Peter Campbell (Nature, 2021), *Lineage tracing and somatic mutation accumulation in human hematopoietic stem cells.*

**实际论文**：Mitchell, Spencer Chapman, Williams et al. (Nature, 2022), *Clonal dynamics of haematopoiesis across the human lifespan*, DOI: 10.1038/s41586-022-04786-y。标题、年份均与引用不符。

**更严重的问题**：此论文研究的是**造血干细胞（HSC）的体细胞突变谱系追踪**——从单细胞衍生的造血克隆中构建系统发育树，推断 HSC 群体大小和克隆动力学。这与 mtDNA 遗传瓶颈估计**没有直接方法学联系**。

gem_method.md 方案 2 将此论文解读为"证明了在个体或细胞克隆演化水平上，将方差偏离度转化为连续定量残差表型去敲核基因组，属于国际最高学术机构的方法学集体共识"，这一说法**过度外推**。Campbell 论文使用的是基于体细胞突变的系统发育推断（MPBoot + ABC 建模），而非异质性漂变方差的残差法。两者的统计框架完全不同：

| 方面 | Campbell (Nature, 2022) | MitoQuest Per-Family Ne |
|---|---|---|
| 数据 | 单细胞 WGS 体细胞突变 | 批量 WGS 母子 mtDNA 异质性 |
| 方法 | 系统发育树 + phylodyn + ABC | Beta-Binomial 复合边际似然 |
| 目标参数 | HSC 群体大小 $N\tau$ | mtDNA 瓶颈 $N_e$ |
| 统计框架 | Coalescent theory | Kimura diffusion + Beta conjugacy |

**结论**：Campbell 论文**不应作为 per-family $N_e$ 估计的方法学背书**，但作为"个体特异性表型提取"的广义概念参考仍可保留。

### 2.3 三篇论文对方法的关键启示

#### 2.3.1 deCODE $N_e \approx 3$：对我们方法的根本性挑战

deCODE (Cell, 2024) 使用 2,548 个冰岛母系家系（116,663 对母子传递对）估计 mtDNA 瓶颈大小，得到 $N_e \approx 3$。这与 Li et al. (2016) 的 ~28-35 估计以及我们现有方法在测试数据上得到的 ~30 范围**相差约一个数量级**。

**可能原因分析**：

1. **胚系选择偏差**：deCODE 发现有害突变在传递中被主动清除（germline selection）。我们的 Kimura 模型假设**中性漂变**，如果存在针对特定突变类别的选择压力，估计的 $N_e$ 会被系统性膨胀——我们观测到的变异传递变异包含了中性漂变 + 选择残余，模型将两者混淆为"更宽松的瓶颈"。

2. **位点组成差异**：deCODE 的估计基于大量母子传递对中的新发突变（de novo mutations），而我们依赖的是**母亲已有的异质性位点**的传递频率变化。两类数据对瓶颈的敏感度不同。

3. **数据量差异**：Li et al. (2016) 使用 ~200 对母子对（深度测序），deCODE 使用 116,663 对（大规模但可能深度较低）。估计量在大样本下可能更准确。

**对我们方法的影响**：

- 如果真实 $N_e \approx 3$，我们当前的 MMLE 和 Kimura 估计量可能存在**系统性向上偏差**
- 我们的 per-family 估计在 $N_e = 3$ 时，$\text{Var}(p_C | p_M) = p_M(1-p_M)/3 \approx 0.083$（$p_M = 0.5$ 时），方差极大
- 这意味着即使单个位点的一次传递就包含了大量信息（方差大 → 对 $N_e$ 敏感），但也意味着位点间方差极大 → CI 很宽
- **关键启示**：我们需要在结果中明确标注 $N_e$ 估计是基于**中性漂变假设**的表观瓶颈大小，而非真实的物理瓶颈单位数

> **D12（新增决策）**：$N_e$ 应被解读为"表观有效瓶颈大小（apparent effective bottleneck size）"——在中性漂变假设下，与实际传递方差最匹配的等效 $N_e$。如果存在胚系选择，估计值将反映"选择 + 漂变"的联合效应。

#### 2.3.2 Chinnery 团队研究路线：胚系选择 + 核基因 GWAS

gem_method.md 引用的 Chinnery 团队工作实际对应两篇论文（见勘误 7），分别提供不同的方法学启示：

**A. Rath et al. (Science, 2019) — 胚系选择证据**

此论文直接研究了母-子 mtDNA 传递过程，与我们的数据和方法最相关：

| 方面 | Rath et al. (Science, 2019) | MitoQuest Per-Family Ne |
|---|---|---|
| 数据 | 1,526 对母子传递对（WGS） | 任意母子传递对（WGS） |
| 核心发现 | 胚系选择塑造 mtDNA 传递 | 中性漂变模型估计 $N_e$ |
| 选择证据 | 已知变异 > 新发变异传递概率；核基因控制选择 | 未建模（D12） |
| 对我们的意义 | 证明 Kimura 中性假设被违反 → $N_e$ 为表观值 | |

> **关键启示**：Rath et al. 证明 mtDNA 传递不是纯中性漂变，而是受到核基因控制的胚系选择。这意味着我们的 $N_e$ 估计实际上是“漂变 + 选择”联合效应的表观值（详见 D12）。deCODE (Cell, 2024) 估计 $N_e \approx 3$ 与我们 ~30 的差异，可能正是因为 deCODE 的估计方法更直接捕捉了选择效应。

**B. Laricchia et al. (Nature, 2023) — GWAS 方法学**

| 方面 | Laricchia et al. (Nature, 2023) | MitoQuest 两步残差法 (§3.9) |
|---|---|---|
| GWAS 表型 | 个体**异质性水平** $h_i$（每个位点单独做 GWAS） | 家系**漂变残差** $\epsilon_f$（跨位点汇总后做 GWAS） |
| 单位 | 个体（无母子对要求） | 母子对（需要家系结构） |
| 核基因发现 | 42 个关联（SSBP1, TFAM, DGUOK, POLG2 等） | 预期（尚未实施） |
| 正态化 | Inverse Rank Normalization (INT) | 对数转换 $\log(D_f + \delta)$ |
| 协变量处理 | 残差法（先回归取残差再做 GWAS） | GLM 回归（在 GWAS 之前回归） |

**启示与建议**：

1. **两种表型互补**：Laricchia et al. 的异质性水平反映“核基因对特定 mtDNA 变异频率的控制”；我们的漂变残差反映“核基因对传递瓶颈强度的控制”。两者捕捉不同的生物学信号。

2. **正态化应使用 INT**：Laricchia et al. 使用 INT 而非对数转换。INT 是 GWAS 表型的标准做法（确保严格正态分布），我们应在 §3.9 中将 $\log(D_f + \delta)$ **替换为 INT** 作为默认选项，同时保留对数转换作为备选。

3. **协变量处理应分两步**：Laricchia et al. 的做法更严谨——先回归协变量取残差，再将残差做 INT 后输入 GWAS。我们的 §3.9 Step 3 应**拆分**为两步：(1) GLM 回归取残差 → (2) INT 转换 → (3) 输出。

> **D13（新增决策）**：两步残差法的正态化默认使用 Inverse Rank Normalization（与 Laricchia et al. 2023 一致），对数转换作为备选。

#### 2.3.3 Campbell 论文的方法学教训

虽然 Campbell (Nature, 2022) 与 mtDNA 瓶颈无直接关联，但其方法论仍有两点可借鉴：

1. **系统发育树作为信息放大器**：Campbell 利用多个单细胞克隆的系统发育关系（而非独立的 pairwise 比较）来推断群体参数。这暗示我们未来可考虑利用**多位点的联合似然**（而非独立位点的复合似然）来提高 $N_e$ 估计精度——前提是建模 mtDNA 位点间的连锁结构。

2. **ABC 建模的价值**：Campbell 使用 Approximate Bayesian Computation 来拟合复杂的群体动力学模型。当 Kimura 扩散的闭式解不再适用时（如引入选择压力），ABC 可能是一个可行的替代推断框架。

---

## 3. 方法学：Per-Family $N_e$ 估计的完整数学推导

### 3.1 统计模型

#### 3.1.1 核心假设

对于一个家系 $f$（包含一个母亲和 $J$ 个孩子），我们做以下假设：

**假设 A1**（共同瓶颈）：家系 $f$ 中所有孩子经历的 mtDNA 传递瓶颈大小相同，记为 $N_{e,f}$。这是该家系的固有属性。

**假设 A2**（位点独立性）：mtDNA 基因组上不同变异位点在传递过程中独立经历 Kimura 漂变。即对位点 $i \neq j$，$(p_{C,i} - p_{M,i})$ 与 $(p_{C,j} - p_{M,j})$ 条件独立（给定各自的 $p_{M,i}$ 和 $p_{M,j}$）。

**假设 A3**（独立测序噪声）：不同位点的测序抽样相互独立。

**假设 A4**（Kimura 平稳近似）：经 bottleneck + 营养期分裂后，孩子异质性的分布由 Kimura 扩散的平稳 Beta 分布描述：

$$
p_{C,j,i} \mid p_{M,i} \sim \text{Beta}\bigl(p_{M,i}(N_{e,f}-1),\; (1-p_{M,i})(N_{e,f}-1)\bigr)
$$

#### 3.1.2 两层生成模型

对于家系 $f$，其生成过程为：

**第一层：母亲异质性**
$$
p_{M,i} \text{ 为母亲在位点 } i \text{ 的真实等位基因频率（非随机，视为固定参数或 plug-in 估计）}
$$

**第二层：传递（对每个孩子 $j$ 的每个位点 $i$）**
$$
p_{C,j,i} \mid p_{M,i}, N_{e,f} \sim \text{Beta}\bigl(p_{M,i}(N_{e,f}-1),\; (1-p_{M,i})(N_{e,f}-1)\bigr)
$$

**第三层：测序**
$$
k_{C,j,i} \mid p_{C,j,i} \sim \text{Binomial}(d_{C,j,i},\; p_{C,j,i})
$$

### 3.2 Per-Family 边际似然推导

#### 3.2.1 单位点、单孩子边际似然

对位点 $i$、孩子 $j$，边际化潜变量 $p_{C,j,i}$：

$$
P(k_{C,j,i} \mid p_{M,i}, N_{e,f}) = \int_0^1 \binom{d_{C,j,i}}{k_{C,j,i}} p^k (1-p)^{d-k} \cdot \frac{p^{\alpha_i-1}(1-p)^{\beta_i-1}}{B(\alpha_i, \beta_i)} \, dp
$$

其中 $\alpha_i = p_{M,i}(N_{e,f}-1)$，$\beta_i = (1-p_{M,i})(N_{e,f}-1)$，$d = d_{C,j,i}$，$k = k_{C,j,i}$。

合并指数后：

$$
= \binom{d_{C,j,i}}{k_{C,j,i}} \frac{1}{B(\alpha_i, \beta_i)} \int_0^1 p^{\alpha_i + k - 1}(1-p)^{\beta_i + d - k - 1} \, dp
$$

积分恰为 Beta 函数 $B(\alpha_i + k, \beta_i + d - k)$，因此：

$$
\boxed{P(k_{C,j,i} \mid p_{M,i}, N_{e,f}) = \binom{d_{C,j,i}}{k_{C,j,i}} \frac{B(\alpha_i + k_{C,j,i},\; \beta_i + d_{C,j,i} - k_{C,j,i})}{B(\alpha_i,\; \beta_i)}}
$$

这就是 **Beta-Binomial 边际分布**，也是 MitoQuest 连续模型的核心似然函数。

#### 3.2.2 对数域实现

$$
\ell_{j,i}(N_{e,f}) = \log\binom{d_{C,j,i}}{k_{C,j,i}} + \log B(\alpha_i + k_{C,j,i},\; \beta_i + d_{C,j,i} - k_{C,j,i}) - \log B(\alpha_i, \beta_i)
$$

使用 `std::lgamma` 稳定计算：

```
log_B_ratio = lgamma(α_i + k) + lgamma(β_i + d - k) - lgamma(α_i + β_i + d)
            - lgamma(α_i) - lgamma(β_i) + lgamma(α_i + β_i)
```

#### 3.2.3 家系复合对数似然

对家系 $f$（母亲 + $J$ 个孩子，$K$ 个变异位点），在位点独立性和孩子独立性假设下：

$$
\boxed{\mathcal{L}_f(N_{e,f}) = \sum_{j=1}^{J} \sum_{i=1}^{K} \ell_{j,i}(N_{e,f})}
$$

其中仅**母亲异质性位点**（$0 < p_{M,i} < 1$，即 $0 < k_{M,i} < d_{M,i}$）对 $N_{e,f}$ 有信息量。母亲纯质位点（$p_{M,i} = 0$ 或 $1$）的边际似然退化为与 $N_{e,f}$ 无关的常数，对优化无贡献。

#### 3.2.4 多孩子家系的特殊处理

当一个母亲有多个孩子时，所有孩子共享同一个 $N_{e,f}$。第 $j$ 个孩子的贡献为：

$$
\mathcal{L}_{f,j}(N_{e,f}) = \sum_{i=1}^{K} \ell_{j,i}(N_{e,f})
$$

家系总似然为各孩子似然之和。这等价于将同一母亲-不同孩子的传递对视为来自同一 bottleneck 过程的**独立重复实验**。

### 3.3 三代家系（G-M-C Trio）扩展

当祖母数据可用时，三代边际似然已在 v1.10.0 中推导为闭式解。对家系 $f$ 中的 trio 行（位点 $i$）：

$$
I_{\text{trio},i}(N_{e,f}) = \binom{d_{M,i}}{k_{M,i}} \binom{d_{C,i}}{k_{C,i}} \cdot \frac{B(\alpha_{G,i} + k_{M,i} + k_{C,i},\; \beta_{G,i} + (d_{M,i}-k_{M,i}) + (d_{C,i}-k_{C,i}))}{B(\alpha_{G,i},\; \beta_{G,i})}
$$

其中 $\alpha_{G,i} = \hat{p}_{G,i}(N_{e,f}-1)$，$\beta_{G,i} = (1-\hat{p}_{G,i})(N_{e,f}-1)$，$\hat{p}_{G,i} = g_{\text{ad\_alt},i} / g_{\text{dp},i}$。

**推导过程**（完整）：

祖母 → 母亲的传递使用 Kimura 平稳 Beta：
$$
p_{M,i} \mid \hat{p}_{G,i} \sim \text{Beta}(\alpha_{G,i},\; \beta_{G,i})
$$

给定 $p_{M,i}$，母亲读数和孩子读数的条件联合为：
$$
P(k_{M,i}, k_{C,i} \mid p_{M,i}, N_{e,f}) = \text{Bin}(k_{M,i} \mid d_{M,i}, p_{M,i}) \times \text{BetaBin}(k_{C,i} \mid d_{C,i}, p_{M,i}(N_{e,f}-1), (1-p_{M,i})(N_{e,f}-1))
$$

对 $p_{M,i}$ 边际化：
$$
I_{\text{trio},i}(N_{e,f}) = \int_0^1 \text{Beta}(p_{M,i} \mid \alpha_{G,i}, \beta_{G,i}) \cdot \text{Bin}(k_{M,i} \mid d_{M,i}, p_{M,i}) \cdot \text{BetaBin}(k_{C,i} \mid d_{C,i}, \alpha_{M,i}, \beta_{M,i}) \, dp_{M,i}
$$

展开各因子：
- Beta 先验核：$p_{M,i}^{\alpha_{G,i}-1}(1-p_{M,i})^{\beta_{G,i}-1}$
- 二项似然核：$p_{M,i}^{k_{M,i}}(1-p_{M,i})^{d_{M,i}-k_{M,i}}$
- BetaBinomial 边际核（已在 §3.2.1 中推导）：吸收后等价于 $p_{M,i}^{k_{C,i}}(1-p_{M,i})^{d_{C,i}-k_{C,i}}$

合并 $p_{M,i}$ 的总指数：
$$
p_{M,i}^{\alpha_{G,i}+k_{M,i}+k_{C,i}-1} \cdot (1-p_{M,i})^{\beta_{G,i}+(d_{M,i}-k_{M,i})+(d_{C,i}-k_{C,i})-1}
$$

这恰是 $\text{Beta}(\alpha_{G,i}+k_{M,i}+k_{C,i},\; \beta_{G,i}+(d_{M,i}-k_{M,i})+(d_{C,i}-k_{C,i}))$ 的核函数，积分即得 Beta 函数比。$\blacksquare$

家系 $f$ 的混合似然（trio 行 + pair 行）：

$$
\mathcal{L}_f(N_{e,f}) = \sum_{i \in \text{trio}} \log I_{\text{trio},i}(N_{e,f}) + \sum_{i \in \text{pair}} \ell_i(N_{e,f})
$$

其中 trio 行自动回退到 pair 似然（当 `has_g == 0`）。

### 3.4 Per-Family MMLE 估计量

家系 $f$ 的 Maximum Marginal Likelihood Estimator（MMLE）定义为：

$$
\boxed{\hat{N}_{e,f} = \arg\max_{N_e \in [1+\epsilon,\; N_{e,\max}]} \mathcal{L}_f(N_e)}
$$

#### 3.4.1 优化策略

采用与现有群体级估计相同的**两阶段优化**：

1. **Phase 1 — 粗扫**：整数网格扫描 $N_e \in \{\text{min\_ne}, \ldots, \text{max\_ne}\}$，定位最佳整数 $N_e^*$
2. **Phase 2 — 精化**：在 $[N_e^*-1, N_e^*+1]$ 区间内 golden-section search

由于家系内位点数量 $K$ 远小于群体总位点数（$K \ll K \times F$），per-family 优化的计算量远小于群体级优化。

#### 3.4.2 搜索边界约束

- 下界 $N_{e,\min} = 1 + 10^{-3}$：避免 $\lgamma(0)$ 奇点（$N_e=1$ 时 Beta 退化）
- $N_e = 1$ 精确值：回退到离散模型（完全漂变 → 固定/丢失）
- 上界 $N_{e,\max}$：默认 200（用户可调）

### 3.5 Per-Family 95% 置信区间

使用 **profile likelihood** 方法（Wilks 定理）：

$$
\text{CI}_{95\%}(f) = \left\{ N_e : \mathcal{L}_f(N_e) \geq \mathcal{L}_f(\hat{N}_{e,f}) - \frac{\chi^2_{1,0.95}}{2} \right\}
$$

其中 $\chi^2_{1,0.95}/2 = 3.841/2 \approx 1.92$。

从 $\hat{N}_{e,f}$ 向两侧步进，找到对数似然下降 1.92 的边界点。

**注意**：当 $K$ 较小时（如 $K < 10$），Wilks 定理的渐近近似可能不精确。此时应考虑：
- 报告更宽的保守区间
- 标注小样本警告

### 3.6 Per-Family Kimura 交叉校验

对家系 $f$，Kimura 矩估计可直接在**家系内**计算：

#### 3.6.1 位点级统计量

对家系 $f$ 中每个信息位点 $i$（$0 < p_{M,i} < 1$）和孩子 $j$：

$$
\begin{aligned}
d_{j,i} &= (p_{C,j,i} - p_{M,i})^2 & &\text{（原始漂移平方）} \\
s_{j,i} &= \frac{p_{M,i}(1-p_{M,i})}{d_{M,i}} + \frac{p_{C,j,i}(1-p_{C,j,i})}{d_{C,j,i}} & &\text{（Wonnapinij 采样噪声校正）} \\
r_{j,i} &= d_{j,i} - s_{j,i} & &\text{（校正后漂移）} \\
w_{i}   &= p_{M,i}(1-p_{M,i}) & &\text{（方差归一化因子）}
\end{aligned}
$$

#### 3.6.2 家系 Kimura 估计量

$$
V_f = \frac{\displaystyle\sum_{j,i} r_{j,i}}{\displaystyle\sum_{j,i} w_i}
$$

$$
b_f = 1 - V_f
$$

$$
\boxed{N_{e,f}^{\text{Kimura}} = \frac{1}{1 - b_f} = \frac{1}{V_f}}
$$

$b_f$ 裁剪到 $(\epsilon, 1-\epsilon)$（$\epsilon = 10^{-9}$）以防护有限样本离群值。

#### 3.6.3 位点级 Bootstrap 置信区间

对家系 $f$：
1. 从家系内的信息位点中**有放回抽样** $B$ 次（每次抽 $K_f$ 个位点）
2. 对每次 bootstrap 样本重新计算 $b_f^{(b)}$
3. 取 2.5% 和 97.5% 分位数作为 95% CI

**注意**：当 $K_f$ 很小时（如 $< 5$），bootstrap 可能产生退化结果（所有 resample 相同），应标注警告。

#### 3.6.4 Trimmed Kimura（可选）

与群体级相同：按 $F_i = r_i / w_i$ 降序排列，丢弃 top `trim_frac` 的高漂移位点，在剩余位点上重新计算 $b_f$。用于诊断 NUMT / 测序错误污染。

### 3.7 信息量分析

#### 3.7.1 最少位点要求

理论上，至少需要 **2 个信息位点**才能估计 $N_e$（1 个自由度参数）。实践中：

| $K_f$（信息位点数） | 估计可靠性 | 建议 |
|---|---|---|
| $K_f < 3$ | 不可靠 | 报告警告；仅输出点估计，CI 可能退化 |
| $3 \leq K_f < 10$ | 可用但宽泛 | 正常输出；标注 "small sample" |
| $K_f \geq 10$ | 可靠 | 正常输出 |
| $K_f \geq 50$ | 高精度 | 与群体级估计可比 |

#### 3.7.2 Fisher 信息量

对单个位点 $i$（简化为单孩子），Fisher 信息量为：

$$
\mathcal{I}_i(N_e) = -E\left[\frac{\partial^2}{\partial N_e^2} \ell_i(N_e)\right]
$$

家系的总信息量为各位点信息量之和：

$$
\mathcal{I}_f(N_e) = \sum_{i=1}^{K_f} \mathcal{I}_i(N_e)
$$

$N_e$ 的 Cramér-Rao 下界为 $\text{Var}(\hat{N}_{e,f}) \geq 1/\mathcal{I}_f(N_e)$。信息量随 $K_f$ 线性增长，因此更多异质性位点 → 更精确的估计。

#### 3.7.3 信息量最大的位点

- $p_{M,i}$ 接近 0.5 的位点携带最多信息（$w_i = p(1-p)$ 最大）
- 高深度位点（$d_{M,i}$, $d_{C,j,i}$ 大）的采样噪声更小
- 这解释了 `--min-vaf 0.10 --max-vaf 0.90` 的默认 VAF 窗口：过滤低信息量位点

### 3.8 Per-Family Kimura 矩估计（含采样噪声校正）

> 本节修正并扩展 `gem_method.md` 方案 1 的矩估计方法，补充 Wonnapinij (2010) 采样噪声校正。

#### 3.8.1 未校正矩估计的偏差分析

`gem_method.md` 方案 1 定义的 $b_i = (p_{Oi}-p_{Mi})^2 / [p_{Mi}(1-p_{Mi})]$ 的期望值为：

$$E[b_i] = E\left[\frac{(p_{Oi}-p_{Mi})^2}{p_{Mi}(1-p_{Mi})}\right] = \underbrace{\frac{1}{N_e}}_{\text{真实漂变}} + \underbrace{\frac{1}{d_{Mi}} + \frac{1}{d_{Ci}} \cdot \frac{E[p_{Ci}(1-p_{Ci})]}{p_{Mi}(1-p_{Mi})}}_{\text{测序二项抽样噪声}}$$

**推导**：$p_{Oi} = p_{Ci}^{\text{true}} + \epsilon_C$（$\epsilon_C$ 为二项抽样噪声，$\text{Var}(\epsilon_C) = p_C^{\text{true}}(1-p_C^{\text{true}})/d_C$），类似 $p_{Mi} = p_{Mi}^{\text{true}} + \epsilon_M$。展开 $(p_{Oi}-p_{Mi})^2$ 的期望即得上式。

因此，**未校正估计 $\hat{N}_e^{\text{naive}} = 1/\bar{b}$ 系统性低估 $N_e$**（因为分母被采样噪声膨胀）。

#### 3.8.2 Wonnapinij 校正后的 Per-Family 矩估计

对家系 $f$ 中每个信息位点 $i$、孩子 $j$，计算：

$$
\begin{aligned}
d_{j,i} &= (p_{C,j,i} - p_{M,i})^2 \\
s_{j,i} &= \frac{p_{M,i}(1-p_{M,i})}{d_{M,i}} + \frac{p_{C,j,i}(1-p_{C,j,i})}{d_{C,j,i}} \\
r_{j,i} &= d_{j,i} - s_{j,i} \\
w_i &= p_{M,i}(1-p_{M,i})
\end{aligned}
$$

家系 Kimura 估计量为加权最小二乘解：

$$\boxed{\hat{N}_{e,f}^{\text{Kimura}} = \frac{\displaystyle\sum_{j,i} w_i^2}{\displaystyle\sum_{j,i} r_{j,i} \cdot w_i} = \frac{\displaystyle\sum_{j,i} w_i}{\displaystyle\sum_{j,i} r_{j,i}}}$$

（第二个等式在 $w_i$ 不随 $j$ 变化时成立——同一母亲的所有孩子共享 $w_i$。）

> **注**：此公式与 §3.6.2 等价。此处重新展开是为了与 `gem_method.md` 的未校正版本做显式对比。

#### 3.8.3 修正后的计算实例

复用 `gem_method.md` 的数据矩阵，假设 500× 均匀深度（$d_M = d_C = 500$）：

| 位点 | $p_M$ | $p_{C,j}$ | $d_i$ | $s_i$（500×） | $r_i$ | $w_i$ |
|---|---|---|---|---|---|---|
| 1-C1 | 0.20 | 0.25 | 0.002500 | 0.000640 | 0.001860 | 0.1600 |
| 1-C2 | 0.20 | 0.12 | 0.006400 | 0.000614 | 0.005786 | 0.1600 |
| 2-C1 | 0.10 | 0.05 | 0.002500 | 0.000990 | 0.001510 | 0.0900 |
| 3-C1 | 0.50 | 0.62 | 0.014400 | 0.000953 | 0.013447 | 0.2500 |
| 3-C2 | 0.50 | 0.40 | 0.010000 | 0.001000 | 0.009000 | 0.2500 |

采样校正后：
- $\sum r_i = 0.031603$，$\sum w_i = 0.9100$
- $\hat{N}_e^{\text{Kimura}} = 0.9100 / 0.031603 \approx 28.8$

对比未校正：$\hat{N}_e^{\text{naive}} = 1/0.0362 \approx 27.6$

> 差异约 4%（500× 深度下采样噪声较小）。在 100× 深度下，差异将扩大到 ~15-20%。

### 3.9 两步残差法：面向 GWAS 的下游表型提取

> 本节修正并扩展 `gem_method.md` 方案 2，将其定位为 **GWAS 下游分析模块**（非独立的 $N_e$ 估计器）。

#### 3.9.1 动机

当单个家系的信息位点数 $K_f$ 很少（如 $K_f < 3$，这在单子女家系中很常见）时，直接 per-family MMLE 的置信区间极宽，估计不可靠。此时，一种替代策略是：

1. 利用**全队列**的统计效能建立群体基准
2. 计算每个家系偏离群体基准的**标准化残差**
3. 将残差作为连续表型用于核基因组 GWAS

这与 Laricchia et al. (Nature, 2023) 研究核基因对 mtDNA 异质性控制的路线一致，但我们使用**传递漂变残差**而非**异质性水平**作为 GWAS 表型——两者捕捉不同的生物学信号（详见 §2.3.2）。

#### 3.9.2 修正后的两步残差法

**Step 1：群体基准估计**

使用 MitoQuest `ne-estimate`（群体级模式）拟合队列共享 $N_e^{\text{pop}}$。

**Step 2：Per-Family 采样校正漂变指数**

对家系 $f$，定义**采样校正后的标准化漂变指数**：

$$D_f = \frac{1}{K_f} \sum_{i=1}^{K_f} F_{f,i}, \quad F_{f,i} = \frac{r_{f,i}}{w_{f,i}} = \frac{d_{f,i} - s_{f,i}}{p_{M,f,i}(1-p_{M,f,i})}$$

> **关键修正 vs. `gem_method.md`**：使用 $r_i = d_i - s_i$（采样校正后漂移），而非原始 $(p_O-p_M)^2/w$。

$D_f$ 的理论期望为 $1/N_e^{\text{pop}}$。正偏离（$D_f > 1/N_e^{\text{pop}}$）表示该家系漂变强于群体平均（瓶颈更紧），负偏离表示漂变弱于平均（瓶颈更宽松）。

**Step 3：协变量回归与残差提取**

对 $D_f$ 进行正态化转换。默认使用 **Inverse Rank Normalization (INT)**（与 Laricchia et al. / Nature, 2023 一致）：

$$Y_f = \Phi^{-1}\left(\frac{\text{rank}(D_f)}{N_f + 1}\right)$$

其中 $\Phi^{-1}$ 为标准正态分布的逆函数，$N_f$ 为家系总数。备选方案为对数转换 $Y_f = \log(D_f + \delta)$（$\delta = 10^{-6}$）。

**Step 3a：协变量回归与残差提取**

拟合多元线性模型：

$$Y_f = \beta_0 + \beta_1 \cdot \overline{d}_{M,f} + \beta_2 \cdot \overline{d}_{C,f} + \beta_3 \cdot \text{Age}_f + \cdots + \epsilon_f$$

其中 $\overline{d}_{M,f}$、$\overline{d}_{C,f}$ 为家系平均测序深度，$\text{Age}_f$ 为母亲生育年龄。**残差 $\epsilon_f$** 即为去除了技术混杂和已知协变量后的**净线粒体遗传控制力表型（Mito-Control Phenotype）**。

> **注**：与 Laricchia et al. (Nature, 2023) 一致，先回归取残差，再做 INT。如果 $D_f$ 已经过 INT，则直接回归取残差即可。

**Step 3b：最终正态化（可选）**

如果残差 $\epsilon_f$ 偏离正态性，可对残差再做一次 INT，确保 GWAS 输入严格符合正态分布。

**Step 4：GWAS 关联分析**

将 $\epsilon_f$ 导出为 PLINK2 标准表型文件，执行核基因组 GWAS：

```
FID  IID  Mito_Control_Phenotype
F001 F001 0.234
F002 F002 -0.156
...
```

#### 3.9.3 与直接 per-family $N_e$ 的关系

| 方法 | 输出 | 统计功效 | 适用场景 |
|---|---|---|---|
| Per-Family MMLE (§3.4) | $\hat{N}_{e,f}$（绝对值 + CI） | 需 $K_f \geq 3$ | 生物学解读："该家系瓶颈是多少" |
| 两步残差法 (§3.9) | $\epsilon_f$（相对残差） | 利用全队列效能 | GWAS 表型："该家系偏离群体多少" |

**推荐用法**：两者**同时报告**。MMLE 提供可解读的绝对估计；残差法为 GWAS 提供正态化表型。

#### 3.9.4 VAF 窗口选择的数学依据

两步残差法中 VAF 窗口的选择（推荐 0.01-0.99 或 0.10-0.90）有三个数学原因：

1. **分母防零**：$w_i = p_M(1-p_M)$ 在 $p_M \to 0$ 或 $1$ 时趋零，导致 $F_i$ 发散
2. **测序错误混淆**：低于测序错误率（~0.1-1%）的信号无法区分真实异质性与技术噪声
3. **纯质位点无信息**：$p_M = 0$ 或 $1$ 时，$E[(p_C-p_M)^2] = 0$，对漂变估计无贡献

### 3.10 三种估计量的关系与选择指南

Per-Family $N_e$ 估计有三种互补的估计量，分别适用于不同场景：

```
                      Per-Family Ne Estimation
                              |
           +------------------+------------------+
           |                  |                  |
     Beta-Binomial      Kimura 矩估计       两步残差法
       MMLE (§3.4)     (§3.8, 校正后)      (§3.9, GWAS)
           |                  |                  |
     主估计器            交叉校验           下游应用
     Cramér-Rao 最优    计算快速           利用群体效能
     需 K >= 3          需 K >= 3          K < 3 也可用
     输出绝对 Ne        输出绝对 Ne        输出相对残差
```

**选择指南**：

| 场景 | 推荐估计量 | 理由 |
|---|---|---|
| $K_f \geq 10$，需精确 Ne | MMLE + Kimura 交叉校验 | 双重验证，最高可靠性 |
| $3 \leq K_f < 10$，需绝对 Ne | MMLE（标注宽 CI） | 唯一可行选择；CI 宽但无偏 |
| $K_f < 3$，需 GWAS 表型 | 两步残差法 | 直接 MMLE 不可靠；借用群体效能 |
| 大样本队列（>1000 家系） | 三者并行 | MMLE 做主表，残差法做 GWAS |

### 3.11 数学推导查漏补缺

#### 3.11.1 `PairData` 缺少家系标识字段

**现状**：`PairData` 结构仅包含 `m_dp, m_ad_alt, c_dp, c_ad_alt, g_dp, g_ad_alt, has_g`，不含 FAM_ID、MOTHER_ID、CHILD_ID。`load_pairs()` 也不读取这三列。

**修复**：

```cpp
struct PairData {
    // ... existing fields ...
    std::string fam_id;       // NEW: family identifier from TSV FAM_ID column
    std::string mother_id;    // NEW: mother sample ID
    std::string child_id;     // NEW: child sample ID
};
```

`load_pairs()` 需要新增对 `FAM_ID`、`MOTHER_ID`、`CHILD_ID` 列的读取（这三列已由 `trans-prep` 输出，见 `trans_prep.cpp` line 72）。

> **注意**：此修改对现有群体级估计完全向后兼容——`fam_id`/`mother_id`/`child_id` 字段仅被 per-family 逻辑使用，群体级代码路径忽略它们。

#### 3.11.2 Bootstrap 单位歧义：多孩子家系的位点级 bootstrap

**问题**：§3.6.3 称 bootstrap 以"位点"为单位有放回抽样。但当母亲有多个孩子时，每个位点 $i$ 对应 $J$ 个孩子的观测 $(k_{C,1,i}, \ldots, k_{C,J,i})$。如果仅按位点抽样，每次 resample 会同时包含该位点所有孩子的数据——这是正确的，但需要明确说明。

**明确定义**：

> **位点级 bootstrap**：从 $K_f$ 个信息位点中有放回抽样 $K_f$ 次。每次抽样以**位点**为原子单位——抽中位点 $i$ 意味着该位点所有孩子的 $(k_{C,j,i}, d_{C,j,i})$ 数据作为一个整体被选中。

这确保了 bootstrap 的 resample 单元与似然函数的独立单元一致（假设 A2：位点间独立）。

#### 3.11.3 Empirical Bayes 收缩估计（补充）

`gem_method.md` 方案 1 第 4 节提到"经验贝叶斯收缩（Empirical Bayes shrinkage）"，但 §3 方法推导中遗漏了这一重要技术。对于小 $K_f$ 家系，直接使用 MMLE 的方差极大，可通过群体先验进行收缩：

**模型**：假设群体级 $N_e$ 估计的先验为 $\log N_{e,f} \sim \mathcal{N}(\mu_{\text{pop}}, \sigma^2_{\text{pop}})$，其中 $\mu_{\text{pop}}, \sigma^2_{\text{pop}}$ 从群体级估计中获取。

**收缩公式**：

$$\hat{N}_{e,f}^{\text{EB}} = \exp\left( w_f \cdot \log \hat{N}_{e,f}^{\text{MMLE}} + (1 - w_f) \cdot \mu_{\text{pop}} \right)$$

其中收缩权重 $w_f = \frac{\mathcal{I}_f}{\mathcal{I}_f + 1/\sigma^2_{\text{pop}}}$，$\mathcal{I}_f$ 为家系 $f$ 的 Fisher 信息量（§3.7.2）。

| 方面 | 直接 MMLE | EB 收缩 |
|---|---|---|
| $K_f \geq 10$ | 可靠，不需收缩 | $w_f \approx 1$，等价于 MMLE |
| $3 \leq K_f < 10$ | CI 宽泛 | $w_f < 1$，向群体均值收缩 |
| $K_f < 3$ | 不可靠 | $w_f \approx 0$，几乎全收缩到群体均值 |

> **D14（新增决策）**：Empirical Bayes 收缩作为**可选后处理**步骤（`--eb-shrink`），不作为默认行为。当启用时，JSON 输出同时报告原始 MMLE 和收缩后的 $\hat{N}_{e,f}^{\text{EB}}$。

#### 3.11.4 INT 的 C++ 实现：逆正态 CDF

§3.9 的 INT 转换需要标准正态分布的逆函数 $\Phi^{-1}$。C++ 标准库不提供此函数，但可通过以下关系实现：

$$\Phi^{-1}(p) = \sqrt{2} \cdot \text{erfinv}(2p - 1)$$

其中 $\text{erfinv}$ 为逆误差函数。实现策略：

1. **使用 `std::erfinv`**（C++26 / 某些编译器扩展）
2. **Rational approximation**：使用 Wicherstein (1993) 或 Beasley-Springer-Moro 算法的高精度有理逼近
3. **Boost 库**：`boost::math::quantile(boost::math::normal_distribution<>(), p)`

**推荐方案**：自行实现 Wicherstein 算法（无外部依赖，精度 $< 10^{-9}$），代码量约 30 行。这在 C++17 标准下即可工作。

```cpp
// Inverse normal CDF via rational approximation (Wicherstein, 1993).
// Accuracy: |error| < 1.15e-9 for p in (0, 1).
static double inv_normal_cdf(double p);
```

#### 3.11.5 GLM 回归的 C++ 实现策略

§3.9 两步残差法的 Step 3a 需要拟合多元线性模型 $Y_f = X\beta + \epsilon$。C++ 实现选项：

| 方案 | 复杂度 | 适用场景 |
|---|---|---|
| **A. 内置 OLS（正规方程）** | 低 | 协变量数 $p \leq 10$；使用 Eigen 或手写 Cholesky 分解 |
| **B. 仅输出 $D_f$，GLM 交给 Python** | 最低 | C++ 端只计算 $D_f$ 和 INT，残差提取交给 `tools/mitoquest_pheno_extractor.py` |
| **C. 调用 R/Python 子进程** | 中 | 不推荐——引入运行时依赖 |

**推荐方案 B**：C++ 端负责统计量计算（$D_f$、INT），Python 脚本负责 GLM 回归和残差提取。理由：

1. OLS 在 C++ 中需要矩阵运算库，增加代码复杂度
2. 协变量可能因项目而异（年龄、BMI、 ethnicity 等），Python 端更灵活
3. `statsmodels` 提供完整的回归诊断（$R^2$、F 统计量、残差正态性检验）
4. 与 `gem_method.md` 已有的 Python 脚本 `mitoquest_pheno_extractor.py` 兼容

**C++ 端输出**：PLINK2 格式的中间表型文件，包含 $D_f$（或 INT($D_f$)）和协变量列：

```
FID  IID  Drift_Index  INT_Drift  Mean_Mother_DP  Mean_Child_DP
F001 F001 0.0312       0.452      487.2           512.1
F002 F002 0.0185       -0.231     501.3           498.7
```

**Python 端**：读取上述文件 + 用户提供的元数据 → GLM → 输出残差表型。

#### 3.11.6 跨家系并行策略

**现状**：现有 `compute_global_ll_parallel()` 在似然评估层并行（多个 pair 的 log-likelihood 求和用线程池）。

**Per-family 并行**：per-family 估计是**尴尬并行**——每个家系的优化完全独立。策略：

1. **家系间并行**（推荐）：将 $F$ 个家系分配到 `threads` 个 worker，每个 worker 串行处理分配到的家系
2. **家系内串行**：单个家系的似然评估不启用线程池（因为 $K_f$ 通常很小，线程开销 > 计算量）

```cpp
// Pseudocode for per-family parallel estimation
std::vector<FamilyResult> results(n_families);
auto task = [&](size_t start, size_t end) {
    for (size_t f = start; f < end; ++f) {
        results[f] = estimate_family(families[f]);
    }
};
parallel_for(task, n_families, config.threads);
```

#### 3.11.7 Kimura 公式中 $w_i$ 在多孩子情况下的求和澄清

§3.6.2 的家系 Kimura 估计量写为：

$$N_{e,f}^{\text{Kimura}} = \frac{\sum_{j,i} w_i}{\sum_{j,i} r_{j,i}}$$

**澄清**：$w_i = p_{M,i}(1-p_{M,i})$ 仅依赖于母亲数据，**不随孩子 $j$ 变化**。因此：

$$\sum_{j=1}^{J} \sum_{i=1}^{K_f} w_i = J \cdot \sum_{i=1}^{K_f} w_i$$

而 $r_{j,i}$ 随孩子 $j$ 变化（因为 $p_{C,j,i}$ 不同）。完整展开为：

$$\hat{N}_{e,f}^{\text{Kimura}} = \frac{J \cdot \sum_{i=1}^{K_f} p_{M,i}(1-p_{M,i})}{\sum_{j=1}^{J} \sum_{i=1}^{K_f} (d_{j,i} - s_{j,i})}$$

这一形式清晰地表明：多孩子家系中，Kimura 分母的权重被孩子数 $J$ 倍放大，而分子的漂变信号来自所有 $J \times K_f$ 个传递事件。

#### 3.11.8 $N_e$ 接近边界时的数值稳定性

**问题**：当 $N_e \to 1^+$ 时，$\alpha_i = p_M(N_e-1) \to 0$，导致 $\lgamma(\alpha_i) \to +\infty$。虽然 §3.4.2 设置了下界 $1 + 10^{-3}$，但在该边界附近：

- $\lgamma(10^{-3}) \approx 6.9$，量级可接受
- 但当 $p_M$ 也接近 0 或 1 时，$\alpha_i = p_M \cdot 10^{-3}$ 可能极小

**防护策略**：

1. 在 `compute_ll_single_continuous` 中，当 $\alpha_i < \epsilon_{\lgamma}$ 或 $\beta_i < \epsilon_{\lgamma}$ 时（$\epsilon_{\lgamma} = 10^{-12}$），返回 0.0（视为非信息位点）
2. CI 计算中，如果左边界触及 $1 + 10^{-3}$，标记 `ci_low_clipped = true`（已有）
3. 添加 stderr 警告：`"WARNING: Ne near lower boundary; Beta parameters may be numerically unstable"`

#### 3.11.9 母亲 plug-in 估计 vs 后验积分的一致性说明

**现状差异**：

| 模型 | $p_M$ 处理方式 |
|---|---|
| 连续模型（默认） | plug-in $\hat{p}_M = m_{\text{alt}} / m_{\text{dp}}$ |
| 离散模型 | Beta 后验积分 $p_M \sim \text{Beta}(m_{\text{alt}}+1, m_{\text{ref}}+1)$ |

这一差异在**per-family 场景中更加显著**——因为家系内位点数 $K_f$ 少，每个位点的 $p_M$ 估计不确定性对 $N_e$ 的影响更大。

**影响分析**：在典型 mtDNA 深度（$\geq 100\times$）下，plug-in 与后验积分的差异量级为 $O(1/d_M)$，可忽略。但在低深度（$< 50\times$）且 $K_f$ 很少时，差异可能不可忽略。

**决策**：保持与现有群体级连续模型一致（plug-in），不引入后验积分。在 JSON 输出中添加 `Mean_Mother_DP` 字段，供用户评估深度是否足够。

#### 3.11.10 Per-family bin-simulation 和 Ne-profile 规格

路线图 Phase 3 提到 per-family bin-simulation，但未给出详细规格。此处补充：

**Per-family bin-simulation**：
- 仅在 $K_f \geq 20$ 时有意义（否则每个 bin 内观测太少）
- 输出格式与群体级相同（`BinSimulationRow`），但仅包含该家系的数据
- 通过 `--per-family-bin-simulation` 标志启用，输出到指定目录
- 每个家系一个 TSV 文件：`{output_dir}/{fam_id}_bin_simulation.tsv`

**Per-family Ne-profile**：
- 对所有家系均可计算（计算量与 $K_f$ 成正比，通常很小）
- 通过 `--per-family-ne-profile` 标志启用
- 每个家系一个 TSV 文件：`{output_dir}/{fam_id}_ne_profile.tsv`
- 可合并为单个文件，增加 `FAM_ID` 列

> **优先级**：Per-family Ne-profile 的优先级高于 bin-simulation（计算量小、诊断价值高）。两者均为 Phase 3 可选功能。

---

## 4. 算法设计与数据结构

### 4.1 家系分组

从 `trans-prep` 输出的 TSV 中读取 `FAM_ID`、`MOTHER_ID`、`CHILD_ID` 列，按 **(FAM_ID, MOTHER_ID, CHILD_ID)** 分组。每个唯一组合定义一个**传递单元**（transmission unit）。

```
家系 F001:
  传递单元 1: FAM_F001, Mother_A, Child_B  → 15 个变异位点
  传递单元 2: FAM_F001, Mother_A, Child_C  → 15 个变异位点（同一母亲的不同孩子）

家系 F002:
  传递单元 1: FAM_F002, Mother_D, Child_E  → 8 个变异位点
```

**关键设计决策**：
- **同一母亲的多孩子**：归入同一家系，共享 $N_{e,f}$
- **祖母行（has_g=1）**：使用 trio 边际似然，与 pair 行混合

### 4.2 新增数据结构

```cpp
// One family's worth of transmission data.
struct FamilyData {
    std::string fam_id;
    std::string mother_id;
    std::vector<std::string> child_ids;  // may have >1 children
    std::vector<PairData>    pairs;      // all variant sites for this family
    
    // Per-family results
    struct FamilyResult {
        std::string fam_id;
        std::string mother_id;
        size_t      n_pairs     = 0;     // total variant sites used
        size_t      n_informative = 0;   // sites with 0 < p_M < 1
        double      ne          = 0.0;   // per-family MMLE Ne
        double      ci_low      = 0.0;
        double      ci_high     = 0.0;
        double      max_log_lik = 0.0;
        bool        ci_low_clipped  = false;
        bool        ci_high_clipped = false;
        KimuraCheck kimura;             // per-family Kimura cross-check
    };
};
```

### 4.3 核心算法流程

```
1. 读取 TSV（含 FAM_ID, MOTHER_ID, CHILD_ID 列）
2. 按 (FAM_ID, MOTHER_ID) 分组
3. 对每个家系 f:
   a. 过滤 QC != "PASS" 的位点
   b. 应用 VAF 窗口 [min_vaf, max_vaf]
   c. 统计信息位点数 n_informative
   d. 如果 n_informative < min_family_sites（默认 3）:
      → 跳过该家系，输出 WARNING
   e. 两阶段优化: Phase 1 整数扫 + Phase 2 golden-section
   f. 计算 95% profile-likelihood CI
   g. [可选] 计算家系内 Kimura 交叉校验 + bootstrap CI
4. 输出 JSON（群体级 + per-family 数组）
```

### 4.4 CLI 设计

```bash
# 现有用法（不变）：群体级 Ne
mitoquest ne-estimate -i mc_pairs.tsv --cross-check kimura -o ne_pop.json

# 新增用法：per-family Ne
mitoquest ne-estimate -i mc_pairs.tsv --per-family --cross-check kimura -o ne_family.json

# 新增用法：per-family Ne + 两步残差法 GWAS 表型
mitoquest ne-estimate -i mc_pairs.tsv --per-family --cross-check kimura \
    --residual-phenotype plink_pheno.txt --metadata maternal_metadata.tsv \
    -o ne_family.json

# 新增参数：
#   --per-family           启用 per-family 估计（默认仍为群体级）
#   --min-family-sites INT 每个家系最少信息位点数 [3]
#   --residual-phenotype FILE  启用两步残差法，输出 PLINK2 表型文件
#   --metadata       FILE  协变量元数据 TSV（含 FAM_ID, Maternal_Age 等）
#   --family-output  FILE  per-family 结果单独输出到指定文件（可选）
```

### 4.5 JSON 输出格式

```json
{
  "Ne":              32.50,
  "CI_95_Low":       28.10,
  "CI_95_High":      38.20,
  "Pairs_Used":      1250,
  "Estimator":       "MMLE (composite marginal likelihood)",
  "Model":           "continuous",
  "Per_Family_Estimates": [
    {
      "FAM_ID":       "F001",
      "Mother_ID":    "Mother_A",
      "N_Children":   2,
      "N_Sites":      30,
      "N_Informative": 18,
      "Ne":           25.30,
      "CI_95_Low":    15.80,
      "CI_95_High":   42.10,
      "Max_Marginal_LogLik": -142.35,
      "Kimura_Cross_Check": {
        "b":           0.0396,
        "Ne_Kimura":   25.25,
        "Sampling_Corrected": true,
        "N_Informative": 18
      }
    },
    {
      "FAM_ID":       "F002",
      "Mother_ID":    "Mother_D",
      "N_Children":   1,
      "N_Sites":      12,
      "N_Informative": 8,
      "Ne":           48.70,
      "CI_95_Low":    22.50,
      "CI_95_High":   100.00,
      "CI_High_Clipped": true,
      "Max_Marginal_LogLik": -58.12,
      "Kimura_Cross_Check": {
        "b":           0.0205,
        "Ne_Kimura":   48.78,
        "Sampling_Corrected": true,
        "N_Informative": 8
      }
    }
  ],
  "Residual_Phenotype_Summary": {
    "N_Families":      4200,
    "Mean_D_f":        0.0312,
    "GLM_Covariates":  ["Mean_Mother_DP", "Mean_Child_DP", "Maternal_Age"],
    "Output_File":     "plink_pheno.txt"
  }
}
```

---

## 5. 技术决策汇总

| # | 决策 | 理由 |
|---|---|---|
| D1 | 家系内位点视为独立观测 | mtDNA 不同位点在传递中独立漂变（假设 A2）；与群体级复合似然框架一致 |
| D2 | 同一母亲的多孩子共享 $N_{e,f}$ | 瓶颈是母系生殖细胞的固有属性，对所有孩子相同 |
| D3 | 使用现有连续 Beta-diffusion 模型 | 与 v1.8.2+ 的群体级模型保持一致；避免离散模型的网格偏差 |
| D4 | Profile likelihood CI（与群体级相同） | Wilks 定理在 $K \geq 10$ 时可靠；小样本时标注警告 |
| D5 | 家系内 Kimura 交叉校验使用位点级 bootstrap（非 pair-level） | 家系内无"pair"概念；bootstrap 单位是位点 |
| D6 | 最少 3 个信息位点才输出估计 | 少于 3 个位点时估计极不稳定，不如不报告 |
| D7 | Per-family 与群体级 $N_e$ 同时输出 | 提供完整信息；用户可自行比较 |
| D8 | `--per-family` 作为显式开关 | 向后兼容；默认行为不变 |
| D9 | 读取 FAM_ID/MOTHER_ID/CHILD_ID 列（已存在于 TSV 中） | 无需修改 `trans-prep` 输出格式 |
| D10 | Kimura 矩估计使用 Wonnapinij 采样校正 | 未校正估计系统性低估 $N_e$（§3.8） |
| D11 | 两步残差法作为下游 GWAS 模块（非主估计器） | 残差法不产生绝对 $N_e$，但与 GWAS 兼容（§3.9） |
| D12 | $N_e$ 解读为「表观有效瓶颈大小」 | deCODE (Cell, 2024) 发现 $N_e \approx 3$ 与我们的 ~30 差异巨大，可能因胚系选择偏差（§2.3.1） |
| D13 | 两步残差法默认使用 INT 正态化 | 与 Laricchia et al. (Nature, 2023) 的 GWAS 表型处理一致（§2.3.2） |
| D14 | EB 收缩作为可选后处理（`--eb-shrink`） | 小 $K_f$ 家系可通过群体先验稳定估计，但不作为默认行为（§3.11.3） |
| D15 | GLM 回归交给 Python 脚本，C++ 端仅输出 $D_f$ + INT | C++ 标准库无 OLS/GLM，且协变量因项目而异，Python 更灵活（§3.11.5） |
| D16 | Per-family 并行采用家系间尴尬并行 | $K_f$ 通常很小，家系内线程开销 > 计算量（§3.11.6） |

---

## 6. 已知限制与风险

| 问题 | 缓解策略 |
|---|---|
| 小 $K_f$ 时估计不可靠 | 设置 `--min-family-sites` 阈值；小样本标注警告；$K_f < 3$ 时自动切换到两步残差法 |
| Profile likelihood CI 在小 $K_f$ 时不准确 | 标注 "asymptotic CI may be unreliable" 警告 |
| 位点间可能存在连锁（LD） | mtDNA 重组率极低，LD 主要来自共同祖先而非重组；复合似然仍一致 |
| 母亲异质性位点数量受测序深度限制 | 高深度 WGS（$\geq 100\times$）通常可检测 10-100 个异质性位点 |
| 母亲纯质（homoplasmic）位点无信息 | 自动跳过；已在信息位点计数中体现 |
| 多位点复合似然假设位点独立 | 若有强烈选择压力，独立性假设可能被违反；实践中 mtDNA 选择信号弱 |
| Kimura 模型假设中性漂变 | deCODE (Cell, 2024) 发现胚系选择证据；$N_e$ 估计为表观值，含选择+漂变联合效应（§2.3.1） |
| 真实 $N_e$ 可能远小于估计值 | deCODE 估计 $N_e \approx 3$，我们得到 ~30；差异可能来自胚系选择、位点组成和数据量（§2.3.1） |
| 两步残差法的正态化选择 | INT 作为默认（与 Laricchia 2023 一致）；对数转换作为备选；大 $D_f$ 时差异可忽略 |
| $N_e$ 边界附近的数值稳定性 | $N_e \to 1^+$ 时 Beta 参数极小，$\lgamma$ 可能不稳定；设置 $\epsilon_{\lgamma}$ 防护（§3.11.8） |
| 母亲 plug-in 估计在低深度下偏差 | 低深度（$< 50\times$）且 $K_f$ 很少时，plug-in 与后验积分差异不可忽略；保持 plug-in 并在输出中报告深度（§3.11.9） |
| EB 收缩依赖群体先验质量 | 如果群体级 $N_e$ 估计本身偏差，收缩会传播偏差；仅作为可选功能（§3.11.3） |

---

## 7. 代码升级方案：详细实现计划

> 基于现有 `ne_estimate.h` (475 行) / `ne_estimate.cpp` (1885 行) 的代码结构，以下为详细的代码修改方案。

### 7.1 Phase 1：数据结构与加载（核心基础）

#### Task 1.1：扩展 `PairData` 结构

**文件**：`src/ne_estimate.h` (line 100-108)

```cpp
struct PairData {
    int m_dp;
    int m_ad_alt;
    int c_dp;
    int c_ad_alt;
    int g_dp     = 0;
    int g_ad_alt = 0;
    int has_g    = 0;
    // NEW: family identifiers (read from TSV, used by per-family mode)
    std::string fam_id;
    std::string mother_id;
    std::string child_id;
};
```

#### Task 1.2：新增 `FamilyData` 和 `FamilyResult` 结构

**文件**：`src/ne_estimate.h` (在 `KimuraCheck` 结构后新增)

```cpp
struct FamilyData {
    std::string fam_id;
    std::string mother_id;
    std::vector<std::string> child_ids;
    std::vector<PairData>    pairs;    // all variant sites for this family
};

struct FamilyResult {
    std::string fam_id;
    std::string mother_id;
    size_t      n_children    = 0;
    size_t      n_pairs       = 0;     // total variant sites used
    size_t      n_informative = 0;     // sites with 0 < p_M < 1
    double      ne            = 0.0;   // per-family MMLE Ne
    double      ci_low        = 0.0;
    double      ci_high       = 0.0;
    double      max_log_lik   = 0.0;
    bool        ci_low_clipped  = false;
    bool        ci_high_clipped = false;
    double      mean_mother_dp  = 0.0; // mean depth for quality assessment
    double      mean_child_dp   = 0.0;
    KimuraCheck kimura;                // per-family Kimura cross-check
    // Optional EB shrinkage
    double      ne_eb_shrunk  = 0.0;   // 0.0 when EB not applied
    bool        eb_applied    = false;
    // Warnings
    std::string warning;               // e.g. "small sample" or "Ne near boundary"
};
```

#### Task 1.3：修改 `load_pairs()` 读取家系列

**文件**：`src/ne_estimate.cpp` (line 643-730)

在 `load_pairs()` 中新增对 `FAM_ID`, `MOTHER_ID`, `CHILD_ID` 列的读取。使用与现有 trio 列相同的 `opt_col()` 策略——列为可选，不影响旧格式 TSV 的加载。

```cpp
// After the existing optional trio column detection:
const int idx_fam_id    = opt_col("FAM_ID");
const int idx_mother_id = opt_col("MOTHER_ID");
const int idx_child_id  = opt_col("CHILD_ID");

// Inside the row parsing loop:
if (idx_fam_id >= 0)    pd.fam_id    = tk[idx_fam_id];
if (idx_mother_id >= 0) pd.mother_id = tk[idx_mother_id];
if (idx_child_id >= 0)  pd.child_id  = tk[idx_child_id];
```

#### Task 1.4：新增 `group_into_families()` 静态方法

**文件**：`src/ne_estimate.h` (声明) + `src/ne_estimate.cpp` (实现)

```cpp
static std::vector<FamilyData>
group_into_families(const std::vector<PairData>& data);
```

实现逻辑：按 `(fam_id, mother_id)` 分组。对于 `fam_id` 为空的行（旧 TSV 格式），退化为单一大组（即群体级模式）。

### 7.2 Phase 2：Per-Family MMLE 估计

#### Task 2.1：新增 `estimate_family()` 静态方法

**文件**：`src/ne_estimate.h` + `src/ne_estimate.cpp`

```cpp
static FamilyResult estimate_family(const FamilyData& fam,
                                     int min_ne, int max_ne,
                                     int min_family_sites);
```

实现核心：
1. 过滤家系内 `QC != PASS` 和非信息位点
2. 检查 `n_informative >= min_family_sites`，否则跳过并设置 warning
3. 调用 `estimate_continuous(fam.pairs, min_ne, max_ne, 1)` （线程=1，家系内串行）
4. 填充 `FamilyResult` 字段
5. 检查 CI clipped、Ne 边界等警告条件

**关键复用**：`estimate_continuous()` 已实现完整的两阶段优化 + profile CI，直接对家系的 `pairs` 子集调用即可。无需重复实现优化逻辑。

#### Task 2.2：新增 `estimate_all_families()` 静态方法

```cpp
static std::vector<FamilyResult>
estimate_all_families(const std::vector<FamilyData>& families,
                      int min_ne, int max_ne,
                      int min_family_sites,
                      int threads);
```

实现：按 §3.11.6 的策略，家系间尴尬并行。使用现有 `ThreadPool` 将家系分配到 worker。

### 7.3 Phase 3：Per-Family Kimura 交叉校验

#### Task 3.1：新增 `compute_family_kimura_check()` 静态方法

```cpp
static KimuraCheck compute_family_kimura_check(
    const FamilyData& fam,
    int n_bootstrap, uint64_t seed,
    double trim_frac);
```

实现与现有 `compute_kimura_check()` 基本相同，但：
1. Bootstrap 单位是**位点**而非 pair（§3.11.2）
2. 多孩子家系的 $w_i$ 求和需按 §3.11.7 的澄清公式计算

### 7.4 Phase 4：CLI 与输出

#### Task 4.1：扩展 `Config` 结构

**文件**：`src/ne_estimate.h` (line 175-209)

```cpp
struct Config {
    // ... existing fields ...
    // NEW: per-family options
    bool        per_family         = false;   // --per-family
    int         min_family_sites   = 3;       // --min-family-sites
    bool        eb_shrink          = false;   // --eb-shrink
    std::string residual_phenotype_file;      // --residual-phenotype
    std::string metadata_file;                // --metadata
    std::string per_family_ne_profile_dir;    // --per-family-ne-profile
};
```

#### Task 4.2：新增 CLI 选项

**文件**：`src/ne_estimate.cpp` (`_parse_args` 和 `usage`)

新增 `getopt_long` 条目：
- `--per-family` (flag 15)
- `--min-family-sites INT` (flag 16)
- `--eb-shrink` (flag 17)
- `--residual-phenotype FILE` (flag 18)
- `--metadata FILE` (flag 19)
- `--per-family-ne-profile DIR` (flag 20)

#### Task 4.3：修改 `run()` 方法

**文件**：`src/ne_estimate.cpp` (line 1544-1884)

在现有群体级估计之后，添加 per-family 分支：

```cpp
Result NeEstimator::run() {
    // ... existing: load_pairs, estimate, kimura_check ...

    // NEW: per-family estimation
    std::vector<FamilyResult> family_results;
    if (_config.per_family) {
        auto families = group_into_families(data);
        family_results = estimate_all_families(
            families, _config.min_ne, _config.max_ne,
            _config.min_family_sites, _config.threads);

        if (_config.kimura_check) {
            for (auto& fr : family_results) {
                // Find corresponding FamilyData and compute Kimura
                // ...
            }
        }

        if (_config.eb_shrink) {
            apply_eb_shrinkage(family_results, r.ne, r.ci_low, r.ci_high);
        }
    }

    // ... existing: bin-simulation, ne-profile ...
    // ... modified: _write_json now includes family_results ...

    // NEW: residual phenotype output
    if (!_config.residual_phenotype_file.empty()) {
        write_residual_phenotype(family_results, data);
    }
}
```

#### Task 4.4：修改 `_write_json()` 输出

**文件**：`src/ne_estimate.cpp` (line 1439-1542)

在现有 JSON 输出的 `Search_Max_Ne` 之后，添加 `Per_Family_Estimates` 数组：

```json
"Per_Family_Estimates": [
  {
    "FAM_ID": "F001",
    "Mother_ID": "Mother_A",
    "N_Children": 2,
    "N_Sites": 30,
    "N_Informative": 18,
    "Ne": 25.30,
    "CI_95_Low": 15.80,
    "CI_95_High": 42.10,
    "Mean_Mother_DP": 487.2,
    "Mean_Child_DP": 512.1,
    "Max_Marginal_LogLik": -142.35,
    "Warning": "",
    "Kimura_Cross_Check": { ... }
  },
  ...
],
"Per_Family_Summary": {
  "N_Families_Estimated": 3850,
  "N_Families_Skipped": 350,
  "Mean_Ne": 31.2,
  "Median_Ne": 28.5,
  "Min_Ne": 2.1,
  "Max_Ne": 198.0
}
```

### 7.5 Phase 5：两步残差法（GWAS 表型提取）

#### Task 5.1：新增 `write_residual_phenotype()` 方法

```cpp
void write_residual_phenotype(
    const std::vector<FamilyResult>& family_results,
    const std::vector<PairData>& data,
    const std::string& output_file);
```

实现：
1. 计算每个家系的 $D_f$（采样校正后的平均漂变指数，§3.9 Step 2）
2. 计算 INT($D_f$)（使用 `inv_normal_cdf`，§3.11.4）
3. 计算家系平均深度等协变量
4. 输出 PLINK2 格式的中间表型文件

#### Task 5.2：更新 `tools/mitoquest_pheno_extractor.py`

更新现有 Python 脚本，使其：
1. 读取 C++ 端输出的中间表型文件
2. 合并用户提供的元数据协变量
3. 执行 GLM 回归（statsmodels OLS）
4. 提取残差 $\epsilon_f$，再 INT 正态化
5. 输出最终 PLINK2 表型文件

### 7.6 Phase 6：测试与验证

#### Task 6.1：新增合成测试夹具

**文件**：`tests/data/ne_pipeline/`

创建 `smoke_per_family.tsv`：包含 3 个家系（F001: 2 孩子 15 位点, F002: 1 孩子 8 位点, F003: 1 孩子 2 位点——测试小样本警告），true Ne=20。

#### Task 6.2：新增单元测试

**文件**：`tests/test_ne_estimate.cpp`

新增 `NeEstFamily` 测试套件：

| 测试用例 | 验证内容 |
|---|---|
| `GroupIntoFamilies.Basic` | 家系分组逻辑正确 |
| `GroupIntoFamilies.EmptyFamId` | 空 FAM_ID 退化为单组 |
| `EstimateFamily.MultiChild` | 多孩子家系 Ne 估计与手算一致 |
| `EstimateFamily.SmallSample` | $K_f < 3$ 时跳过并输出 warning |
| `EstimateFamily.KimuraCrossCheck` | Per-family Kimura 与 MMLE 一致性 |
| `EstimateFamily.TrioSupport` | 三代家系的 trio 似然正确调用 |
| `PerFamilyJson.Output` | JSON 输出格式完整 |
| `ResidualPhenotype.DriftIndex` | $D_f$ 计算与手算一致 |
| `InvNormalCdf.Accuracy` | $\Phi^{-1}$ 精度 $< 10^{-9}$ |
| `EBShrinkage.Convergence` | 大 $K_f$ 时 EB 收敛到 MMLE |

#### Task 6.3：真实数据端到端验证

- 使用 BIGCS_II 队列数据运行完整 per-family 估计
- 验证群体级 Ne 与 per-family Ne 中位数的一致性
- 检查 MMLE 与 Kimura 交叉校验的一致性
- 评估运行时间（预期：4200 家系在 8 线程下 < 5 分钟）

### 7.7 实施优先级与版本规划

| Phase | 内容 | 预估工作量 | 目标版本 |
|---|---|---|---|
| Phase 1 | 数据结构 + 加载 | 1-2 天 | v1.11.0 |
| Phase 2 | Per-Family MMLE | 1-2 天 | v1.11.0 |
| Phase 3 | Per-Family Kimura | 1 天 | v1.11.0 |
| Phase 4 | CLI + 输出 | 1-2 天 | v1.11.0 |
| Phase 5 | 两步残差法 + Python | 2-3 天 | v1.12.0 |
| Phase 6 | 测试 + 验证 | 2-3 天 | v1.11.0 (Phase 1-4) |
| **Total** | | **~2 周** | |

> **建议**：Phase 1-4 + Phase 6 合并为 v1.11.0 发布，Phase 5 (两步残差法) 推迟到 v1.12.0。这样可先交付核心的 per-family Ne 估计功能，GWAS 表型提取作为后续扩展。

---

## 8. 与 Li et al. (2016)、Rath et al. (2019)、deCODE (2024)、Laricchia et al. (2023) 及 `gem_method.md` 方法的对比

| 方面 | Li et al. (2016) | Rath et al. (Science, 2019) | deCODE (Cell, 2024) | Laricchia et al. (Nature, 2023) | `gem_method.md` 方案 1 | `gem_method.md` 方案 2 | MitoQuest Per-Family MMLE |
|---|---|---|---|---|---|---|---|
| 模型 | 可变 bottleneck | 胚系选择 + 传递偏差 | Wright-Fisher + 选择 | 核基因-异质性关联 | Kimura 方差 | 残差法 (vGWAS) | Kimura 平稳 Beta-diffusion |
| 估计方法 | 复合似然 | 母子频率比较 | 母子频率分布 | GWAS（异质性水平） | 矩估计 | GLM 残差 | 复合边际似然 (MMLE) |
| 估计 $N_e$ | 是（~28-35） | 否（研究选择） | 是（$\approx 3$） | 否（估计异质性水平） | 是（未校正） | 否（残差） | 是（表观 $N_e$） |
| 胚系选择 | 未考虑 | **发现强证据** | **发现强证据** | N/A | 未考虑 | 未考虑 | 未建模（D12） |
| 核基因控制 | 未考虑 | **首次证明** | 未明确 | 42 个核基因关联 | 未考虑 | 未考虑 | 可导出 GWAS 表型 |
| 母子对 | ~200 对 | **1,526 对** | 116,663 对 | 无（横截面） | 概念描述 | 概念描述 | 任意 |
| 采样校正 | 未明确 | 未明确 | 未明确 | N/A | **缺失** | **缺失** | Wonnapinij (2010) |
| 三代 | 不支持 | 不支持 | 支持（多代家系） | 不支持 | 不支持 | 不支持 | 支持（G-M-C trio） |
| GWAS 表型 | 无 | 无 | 无 | 异质性水平 + INT | 无 | 残差 $\epsilon$ | 漂变残差 + INT |
| 软件 | 独立脚本 | 未公开 | 未公开 | 未公开 | 概念描述 | Python 脚本 | `ne-estimate` |

---

## 9. 关键文件索引

| 文件 | 当前状态 | 需修改内容 |
|---|---|---|
| `src/ne_estimate.h` | 475 行 | ① `PairData` 新增 3 个 ID 字段 ② 新增 `FamilyData`/`FamilyResult` 结构 ③ 新增 `Config` 字段（per_family 等） ④ 新增方法声明 |
| `src/ne_estimate.cpp` | 1885 行 | ① `load_pairs()` 新增 FAM_ID 读取 (~20 行) ② `group_into_families()` (~40 行) ③ `estimate_family()` (~60 行) ④ `estimate_all_families()` (~40 行) ⑤ `compute_family_kimura_check()` (~80 行) ⑥ `run()` per-family 分支 (~50 行) ⑦ `_write_json()` 输出 (~80 行) ⑧ `_parse_args()` CLI (~40 行) ⑨ `inv_normal_cdf()` (~30 行) |
| `src/trans_prep.cpp` | 无需修改 | FAM_ID/MOTHER_ID/CHILD_ID 已由 `tsv_header()` 输出 |
| `src/trans_prep.h` | 无需修改 | `PairRecord` 已含 fam_id/mother_id/child_id 字段 |
| `tests/test_ne_estimate.cpp` | 新增测试套件 | `NeEstFamily` 测试套件（~10 个测试用例） |
| `tests/data/ne_pipeline/` | 新增测试夹具 | `smoke_per_family.tsv`（3 个家系，true Ne=20） |
| `tools/mitoquest_pheno_extractor.py` | 需更新 | 集成采样校正的 $D_f$ 计算 + INT + GLM 残差提取 |
| `CMakeLists.txt` | 无需修改 | 无新源文件（所有改动在现有文件中） |

---

## 10. 下一阶段入口

→ 实现 per-family 功能后，生成 `HANDOFF_v1.11.0_per_family_ne.md`
