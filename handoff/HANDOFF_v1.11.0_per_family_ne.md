# Handoff — v1.11.0：Per-Family $N_e$ 估计

> 生成时间：2026-06-11  
> 版本标签：`v1.11.0`（已 tag 并推送至 GitHub）  
> 关键 commit：`76b0a07` → `32f3b78` → `9a214f5`  
> 前置阶段：[HANDOFF_per_family_ne_estimation.md](HANDOFF_per_family_ne_estimation.md)（设计文档）、[HANDOFF_v1.10.0_20250528.md](HANDOFF_v1.10.0_20250528.md)

---

## 1. 阶段总览

| 版本 | 日期 | 主要变更 |
|---|---|---|
| v1.11.0 | 2026-06-11 | **新增** `ne-estimate --per-family` 模式：按 (FAM_ID, MOTHER_ID) 分组，对每个家系独立估计 $N_e$ |

---

## 2. 实现概述

### 2.1 功能描述

在群体级 $N_e$ 估计的基础上，新增 **per-family 模式**：

1. 按 `(FAM_ID, MOTHER_ID)` 将传输对分组为家系
2. 对每个家系独立运行连续 Beta-diffusion MMLE 估计
3. 家系间尴尬并行（`std::async`）
4. 信息位点不足的家系自动跳过并输出警告
5. 可选 per-family Kimura 交叉校验

### 2.2 CLI 新增选项

```
      --per-family             Enable per-family Ne estimation.
      --min-family-sites INT   Minimum informative sites per family [3].
      --per-family-output FILE Write per-family Ne results as a TSV file.
```

### 2.3 输出格式

#### Per-Family TSV（18 列）

| 列名 | 说明 |
|---|---|
| `FAM_ID` | 家系标识 |
| `MOTHER_ID` | 母亲样本 ID |
| `N_Children` | 孩子数 |
| `N_Sites` | 总变异位点数 |
| `N_Informative` | 信息位点数（$0 < p_M < 1$） |
| `Ne` | 家系 MMLE $N_e$ 估计 |
| `CI_95_Low` | 95% CI 下界 |
| `CI_95_High` | 95% CI 上界 |
| `Max_LogLik` | 最大对数边际似然 |
| `Mean_Mother_DP` | 母亲平均深度 |
| `Mean_Child_DP` | 孩子平均深度 |
| `CI_Low_Clipped` | CI 下界是否被裁剪 |
| `CI_High_Clipped` | CI 上界是否被裁剪 |
| `Skipped` | 是否因位点不足被跳过 |
| `Warning` | 警告信息 |
| `Kimura_b` | Kimura 矩估计 $b$ 值（如启用） |
| `Kimura_Ne` | Kimura $N_e$ 估计（如启用） |
| `Kimura_N_Informative` | Kimura 使用的信息位点数 |

#### JSON 输出

- 新增 `Per_Family_Estimates` 数组（每家系一条记录）
- 新增 `Per_Family_Summary` 块（估计/跳过的家系统计）

---

## 3. 新增函数

在 `src/ne_estimate.h` / `src/ne_estimate.cpp` 中新增 6 个静态方法：

| 函数 | 位置 | 说明 |
|---|---|---|
| `group_into_families()` | ne_estimate.cpp:1222-1257 | 按 (FAM_ID, MOTHER_ID) 分组；空 FAM_ID 时退化为单组 "ALL" |
| `estimate_family()` | ne_estimate.cpp:1259-1310 | 单家系 MMLE：过滤 + 信息位点计数 + 调用 `estimate_continuous()` |
| `estimate_all_families()` | ne_estimate.cpp:1312-1348 | 批量估计 + `std::async` 并行 |
| `compute_family_kimura_check()` | ne_estimate.cpp:1350-1359 | 家系内 Kimura 矩估计（Wonnapinij 采样校正） |
| `_write_family_tsv()` | ne_estimate.cpp:1361-1415 | 输出 per-family TSV 文件 |
| `load_pairs()` 扩展 | ne_estimate.cpp:643-730 | 新增 FAM_ID/MOTHER_ID/CHILD_ID 列读取 |

---

## 4. 技术决策

| # | 决策 | 理由 |
|---|---|---|
| D1 | 复用现有连续 MMLE 估计器 | 无需重复实现优化逻辑；家系的 `pairs` 子集直接输入 |
| D2 | 家系间 `std::async` 尴尬并行 | 家系优化完全独立；家系内位点少，不需线程池 |
| D3 | `min_family_sites` 默认 3 | 少于 3 个信息位点时 MMLE 极不稳定 |
| D4 | 跳过家系输出 `ne=0, skipped=true` | 避免报告无意义的估计值 |
| D5 | 小样本（3-9 位点）输出 "small sample" 警告 | Profile likelihood CI 在小样本时渐近近似不精确 |
| D6 | `load_pairs()` 中家系列为可选 | 完全向后兼容旧格式 TSV |
| D7 | `Config cfg{}` 值初始化 | 避免 POD 成员未定义值导致的 UB |

---

## 5. 数据结构变更

### PairData 新增字段

```cpp
struct PairData {
    // ... existing fields ...
    std::string fam_id;       // 从 TSV FAM_ID 列读取
    std::string mother_id;    // 从 TSV MOTHER_ID 列读取
    std::string child_id;     // 从 TSV CHILD_ID 列读取
};
```

### 新增结构

```cpp
struct FamilyData {
    std::string fam_id;
    std::string mother_id;
    std::vector<std::string> child_ids;
    std::vector<PairData>    pairs;
};

struct FamilyResult {
    std::string fam_id, mother_id;
    size_t n_children, n_pairs, n_informative;
    double ne, ci_low, ci_high, max_log_lik;
    double mean_mother_dp, mean_child_dp;
    bool ci_low_clipped, ci_high_clipped, skipped;
    std::string warning;
    KimuraCheck kimura;
};
```

---

## 6. 测试

### 6.1 新增测试套件 `NeEstFamily`（23 项）

| 测试用例 | 验证内容 |
|---|---|
| `GroupIntoFamiliesCorrectGrouping` | 按 (FAM_ID, MOTHER_ID) 正确分组 |
| `GroupIntoFamiliesLegacyFallback` | 无家系列列时退化为 "ALL" |
| `GroupIntoFamiliesUniqueChildIds` | child_ids 去重 |
| `GroupIntoFamiliesEmptyInput` | 空输入返回空 |
| `GroupIntoFamiliesSameFamDiffMother` | 同 FAM 不同母亲分为独立家系 |
| `EstimateFamilyRecoversTrueNe` | 合成数据（true Ne=15, 200 位点）恢复 Ne ±40% |
| `EstimateFamilySkipsSmallFamilies` | 信息位点 < min_family_sites 时跳过 |
| `EstimateFamilySkipsHomoplasmicOnly` | 全纯合母亲 → 0 信息位点 → 跳过 |
| `EstimateFamilyMeanDepth` | 平均深度计算正确 |
| `EstimateFamilyWarningSmallSample` | 3-9 个信息位点时输出 "small sample" 警告 |
| `EstimateAllFamiliesSerialParallelAgree` | 串行/并行结果完全一致（4 家系） |
| `EstimateAllFamiliesMultipleFamilies` | 3 个不同 true Ne 的家系正确排序 |
| `EstimateAllFamiliesMixedSkipped` | 混合有足够/不足位点的家系 |
| `FamilyKimuraConsistency` | 与直接调用 `compute_kimura_check` 结果一致 |
| `FamilyKimuraOnSkippedFamily` | 小样本仍可计算 Kimura |
| `WriteFamilyTsvFormat` | 表头、数据行格式正确，无 Kimura 时不含 Kimura 列 |
| `WriteFamilyTsvWithKimura` | 有 Kimura 时输出额外列 |
| `WriteFamilyTsvSkippedKimuraNA` | 跳过家系输出 NA for Kimura |
| `LoadPairsReadsFamilyColumns` | FAM_ID/MOTHER_ID/CHILD_ID 正确读取 |
| `LoadPairsLegacyNoFamilyColumns` | 无家系列列时字段为空 |
| `EndToEndValidationTsv` | 完整管线：加载 → 分组 → 估计 → F003 被跳过 |
| `EndToEndWithKimura` | 完整管线 + Kimura 交叉校验 |
| `ConfigPerFamilyDefaults` | per-family 配置默认值正确 |

### 6.2 合成测试数据

| 文件 | 说明 |
|---|---|
| `tests/data/ne_pipeline/per_family_validation.tsv` | 3 家系合成数据（F001: 2 孩子 15 位点 true Ne≈10, F002: 1 孩子 20 位点 true Ne≈30, F003: 1 孩子 2 位点——应被跳过） |
| `tests/data/ne_pipeline/smoke_per_family.tsv` | 快速冒烟测试数据 |

### 6.3 CI 测试结果

**全部 130 项测试通过**（107 已有 + 23 新增）。

---

## 7. CI 修复记录

### 问题 1：`NeEstFamily.EndToEndValidationTsv` 和 `EndToEndWithKimura` 失败

**原因**：测试使用路径 `tests/data/ne_pipeline/per_family_validation.tsv`，仅在从项目根目录运行时有效。CI 中 `ctest` 从 `build/tests/` 运行（CMake `file(COPY)` 将测试数据复制到 `build/tests/data/`）。

**初次修复（错误方向）**：添加 `WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}` 到 `add_test()`，使 ctest 从项目根运行。

**后果**：破坏了旧测试（BamTest, MtCopyNumRunE2E 等），它们使用 `data/...` 路径，仅在 `build/tests/` 下有效。

**最终修复**：
1. 还原 `WORKING_DIRECTORY` 修改
2. 将 NeEstFamily 测试路径改为 `data/ne_pipeline/per_family_validation.tsv`，与旧测试一致

**教训**：项目中所有测试统一使用 `data/...` 路径（相对于 `build/tests/`），这是 CMake `file(COPY tests/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})` 的约定。

### 问题 2：Config POD 成员未初始化

**原因**：`NeEstimator::Config cfg;` 默认初始化时 POD 成员（如 `min_vaf`, `kimura_bootstrap`）值未定义。

**修复**：使用 `Config cfg{}` 值初始化（所有 POD 成员零初始化）。

---

## 8. 文件变更清单

| 文件 | 变更 | 行数 |
|---|---|---|
| `src/ne_estimate.h` | 新增 `FamilyData`/`FamilyResult` 结构、`PairData` 扩展、6 个方法声明、`Config` 新字段 | +73 |
| `src/ne_estimate.cpp` | 新增 6 个函数实现、`load_pairs()` 扩展、CLI 选项、JSON 输出、`_write_family_tsv()` | +234 |
| `tests/test_ne_estimate.cpp` | 新增 `NeEstFamily` 测试套件（23 项）、2 个辅助函数 | +633 |
| `tests/data/ne_pipeline/per_family_validation.tsv` | 3 家系合成测试数据 | +44 |
| `tests/data/ne_pipeline/smoke_per_family.tsv` | 快速冒烟测试数据 | 已有 |
| `README.md` | per-family 参数说明、用法章节、TSV 格式表、验证指南 | +156 |
| `CMakeLists.txt` | 版本 1.10.0 → 1.11.0 | ±1 |
| `tests/CMakeLists.txt` | CI 修复（最终还原，无净变更） | 0 |

**总计**：~1140 行新增代码 + 测试 + 文档。

---

## 9. 与设计文档的对照

v1.11.0 实现了 [HANDOFF_per_family_ne_estimation.md](HANDOFF_per_family_ne_estimation.md) 中 Phase 1-4 + Phase 6 的核心功能：

| 设计文档 Phase | 状态 | 说明 |
|---|---|---|
| Phase 1: 数据结构 + 加载 | ✅ 完成 | PairData 扩展、FamilyData/FamilyResult、load_pairs() 扩展 |
| Phase 2: Per-Family MMLE | ✅ 完成 | group_into_families()、estimate_family()、estimate_all_families() |
| Phase 3: Per-Family Kimura | ✅ 完成 | compute_family_kimura_check() |
| Phase 4: CLI + 输出 | ✅ 部分完成 | --per-family、--min-family-sites、--per-family-output、JSON 输出；未实现 --eb-shrink、--residual-phenotype、--metadata |
| Phase 5: 两步残差法 | ⏳ 推迟至 v1.12.0 | INT、GLM 残差提取、PLINK2 表型输出 |
| Phase 6: 测试 + 验证 | ✅ 完成 | 23 项单元测试 + 合成数据验证 |

### 未实现功能（v1.12.0 计划）

- [ ] Empirical Bayes 收缩估计（`--eb-shrink`）
- [ ] 两步残差法 GWAS 表型提取（`--residual-phenotype`、`--metadata`）
- [ ] Per-family Ne-profile（`--per-family-ne-profile`）
- [ ] Per-family bin-simulation（`--per-family-bin-simulation`）
- [ ] INT 正态化的 `inv_normal_cdf()` 实现
- [ ] `tools/mitoquest_pheno_extractor.py` GLM 集成

---

## 10. 已知限制

| 问题 | 缓解策略 |
|---|---|
| 小 $K_f$（< 10）时 CI 不精确 | 标注 "small sample" 警告 |
| Kimura 中性假设可能违反 | $N_e$ 为表观有效瓶颈大小（deCODE, Cell, 2024） |
| 位点间可能存在 LD | mtDNA 重组率极低，实践中影响可忽略 |
| 母亲 plug-in 估计在低深度偏差 | 保持与群体级一致；报告 Mean_Mother_DP 供用户评估 |

---

## 11. 下一步计划

### 短期

- [ ] **发布 GitHub Release v1.11.0**：tag `v1.11.0` 已推送；通过 GitHub Web UI 创建 Release
- [ ] **两步残差法**（Phase 5）：INT 正态化 + PLINK2 表型输出

### 中期（v1.12.0）

- [ ] Empirical Bayes 收缩估计
- [ ] 真实家系数据端到端验证（BIGCS_II 队列）
- [ ] Per-family Ne-profile 和 bin-simulation

### 长期

- [ ] 核基因 GWAS 关联分析（漂变残差表型）
- [ ] 多代谱系扩展（4-gen+）
- [ ] 家系结构自动推断

---

## 12. 关键文件索引

| 文件 | 用途 |
|---|---|
| `src/ne_estimate.h` | NeEstimator 类定义，含 PairData、FamilyData、FamilyResult、Config |
| `src/ne_estimate.cpp` | ne-estimate 完整实现（~2100 行），含 per-family 6 个新函数 |
| `tests/test_ne_estimate.cpp` | NeEstPair + NeEstTrio + NeEstFamily 测试套件（56 项） |
| `tests/data/ne_pipeline/per_family_validation.tsv` | 3 家系合成验证数据 |
| `README.md` | 主文档：per-family 用法、TSV 格式、验证指南 |
| `HANDOFF_per_family_ne_estimation.md` | 完整方法学设计文档（1429 行） |
