## The functions and usage of scripts

#### mito_vqsr.py (原: mito_classifier.py)
##### 1. 功能：通过每个样本中每个位点的AD，HF和HQ值构建高斯混合模型（GMM）对mtDNA的变异VCF文件进行位点上的分类和筛选
##### 2. 构建流程介绍：
- 1. 通过公开数据库中的已报导的变异位点、群体变异频率AF>0.01以及自身数据集AF>0.05的位点作为高质量训练集，构建Good GMM模型。
  - 训练集的相关处理：
    - 1. 对全部数据都进行了Z-score标准化处理
    - 2. 为了使未变异的位点的训练集和变异位点的训练集尽量平衡，对未变异位点的训练集进行了随机采样，使得变异位点的数量和未变异位点的数量一致
    - 3. 剔除训练数据集中在平均值的±6个标准差外的数据
- 2. 使用Good GMM模型对自身数据集进行评分（log10为底），筛选出一个阈值，使得98%高质量训练集中位点的评分大于该阈值，低于这个阈值的位点被认为是低质量的。
- 3. 使用低质量训练集构建Bad GMM模型，并对自身数据集进行评分（log10为底）。
  - 训练集的相关处理：
    - 1. 对全部数据都进行了Z-score标准化处理
    - 2. 剔除训练数据集中在平均值的±6个标准差外的数据
- 4. 使用Good GMM模型和Bad GMM模型的评分结果差值，随后计算每个阈值差对应的高质量集（由Good GMM模型评分大于 ii步骤 中阈值的位点组成）的百分比和低质量集（由Good GMM模型评分小于 ii步骤 中阈值的位点组成）的百分比。
- 5. 确定了3个阈值，分别使得iv步骤中的高质量集的百分比为95%，98%，99%。
- 6. 最终指定阈值为99%，对自身数据集进行分类，输出分类结果。
- 7. 对于位点的分类：如果该位点中50%的评分结果差值都大于 vi步骤 中的阈值则认为该位点是高质量（PASS）, 反之则认为该位点是低质量（LowQual）。

##### 3. 使用方法：
```python
$ cd /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/20241225_BIGCS_mtDNA/output_bq20/03.annotation
$ dir=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/software/mitoquest/tools
$ python $dir/mito_classifier.py -i 328.samples.filtered.annotated.vcf.gz -o 328.samples.filtered.annotated.GMM98_diy.vcf.gz -c 328.samples.filtered.annotated.GMM98_diy.csv -f 328.samples.filtered.annotated.GMM98_diy.pdf -rgm 0.98 -gnc 5 -bnc 7 &

$ python mito_classifier.py -h
usage: mito_classifier.py [-h] -i INPUT_VCF_FILE -o OUTPUT_VCF_FILE -c OUTPUT_CSV_FILE -f FIGURE_FILE [-mnc MAX_N_COMPONENTS]
                          [-gnc GOOD_MODULE_N_COMPONENTS] [-bnc BAD_MODULE_N_COMPONENTS] [-rgm GOODRATIO_IN_GOODMODULE] [-pr PASS_RATIO]
                          [-fmr FINAL_RES_RATIO]

Build GMM model using AD, HF, HQ information of mtDNA data to filter out bad sites values.

options:
  -h, --help            show this help message and exit
  -i INPUT_VCF_FILE, --input_vcf_file INPUT_VCF_FILE # 输入VCF文件路径
                        Input VCF file path
  -o OUTPUT_VCF_FILE, --output_vcf_file OUTPUT_VCF_FILE # 输出VCF文件路径
                        Output VCF file path
  -c OUTPUT_CSV_FILE, --output_csv_file OUTPUT_CSV_FILE # 输出CSV文件路径
                        Output CSV file path
  -f FIGURE_FILE, --figure_file FIGURE_FILE # 输出图文件路径
                        Output figure file path
  -mnc MAX_N_COMPONENTS, --max_n_components MAX_N_COMPONENTS # 构建GMM模型的最大组件数量
                        Maximum components for GMM model training (default: 10)
  -gnc GOOD_MODULE_N_COMPONENTS, --good_module_n_components GOOD_MODULE_N_COMPONENTS # 指定构建Good GMM模块的组件数量
                        Number of components for Good GMM module if not specified, default is auto-selected
  -bnc BAD_MODULE_N_COMPONENTS, --bad_module_n_components BAD_MODULE_N_COMPONENTS # 指定构建Bad GMM模块的组件数量
                        Number of components for Bad GMM moduleif not specified, default is auto-selected
  -rgm GOODRATIO_IN_GOODMODULE, --goodratio_in_goodmodule GOODRATIO_IN_GOODMODULE # Good模块中Good比例阈值
                        Ratio of good values to total good values after Good module to select subdataset (Bad values) to train Bad GMM model (default: 0.98)
  -pr PASS_RATIO, --pass_ratio PASS_RATIO # 输出PASS比例阈值
                        Ratio of Good Allels to total Allels per POS to label a site as PASS or LowQual (default: 0.5)
  -fmr FINAL_RES_RATIO, --final_res_ratio FINAL_RES_RATIO # 输出最终结果的Good比例阈值
                        The Good ratio of the final result is used to determine whether the value is good or bad (default: 0.99)
```
##### 4. 输出文件说明：
- 1. 输出VCF文件：输出的VCF文件中，低质量的位点被标记为LowQual，高质量的位点被标记为PASS。
- 2. 输出CSV文件：包括每个样本每个位点的相关信息。
- 3. 输出Figure文件：包括GMM模型的训练过程和分类结果。
- 4. 输出Log文件：包括每个步骤中数据的相关数量信息。

#### plot_bottleneck_simulation.py
##### 1. Purpose
Reproduces the deCODE 2024 Cell **"Figure 5. Observed and simulated
means for bottleneck parameter, b"** panel from the per-bin TSV emitted
by `mitoquest ne-estimate --bin-simulation`. The figure shows two panels
side by side:

* **Left** — observed mean drift `(p_c - p_m)^2` per maternal-VAF bin
  with error bars, overlaid on the simulated parabolas
  `p_m(1 - p_m) / Ne` at the fitted MLE Ne (red, with 95% CI ribbon) and
  the Wonnapinij/Kimura cross-check Ne (blue dashed).
* **Right** — per-bin estimate of `1 - b = F_i = (d_i - s_i) / [p_m(1 - p_m)]`
  with error bars, overlaid on the horizontal lines `1/Ne_MLE` (red
  with CI band) and `1/Ne_Kimura` (blue dashed).

Marker size scales with the number of pairs in each bin so visual
weight is proportional to information content, exactly as in the
deCODE figure.

##### 2. Upstream input — generate the TSV with `mitoquest ne-estimate`
```bash
# Build mitoquest first (cmake --build build -j8), then:
mitoquest ne-estimate \
    -i  cohort.transmission_pairs.tsv \
    -o  cohort.ne.json \
    --cross-check kimura \
    --kimura-bootstrap 1000 \
    --bin-simulation       cohort.bin_sim.tsv \
    --bin-simulation-bins  10
```
The TSV stores per-bin `obs_var`, `obs_var_corr`, `obs_F`, `obs_F_se`
plus the theoretical curves at the fitted Ne and its 95% CI; a
commented header records the fitted Ne, CI bounds, and (when Kimura
is enabled) the Wonnapinij b and Ne_Kimura.

##### 3. Plot it
```bash
python tools/plot_bottleneck_simulation.py \
    -i  cohort.bin_sim.tsv \
    -o  cohort.bottleneck.png
```
Key options: `--dpi` (default 300), `--figsize W,H` (default `13,5.2`),
`--title` (overrides the auto-generated suptitle).

##### 4. How to read it
* If the observed bin means **track the red MLE parabola**, the data
  are consistent with the fitted Ne and the single-generation
  Wright-Fisher model is a good fit.
* If the observed means **fall between the red and blue curves**, the
  MLE and Wonnapinij b are picking up different aspects of the same
  variance (typically the Kimura b is pulled by a few high-drift
  outliers; rerun `mitoquest ne-estimate` with `--kimura-trim 0.10
  --top-drift-k 20` to inspect).
* If the observed means **deviate systematically from the parabolic
  shape** (e.g. a linear or U-shaped pattern), the single-generation
  model is misspecified — likely candidates are deep pedigrees (`g > 1`),
  NUMT contamination, or strongly heterogeneous Ne across sites.


#### plot_ne_profile.py
##### 1. Purpose
Reproduces the deCODE 2024 *Cell* paper's "best-fit Ne" exercise: the
paper scanned candidate Ne values, simulated the Kimura distribution at
each one, and chose the Ne that best fits the observed allele-frequency
change distribution (Ne ≈ 3 across 137 variants × 53,041 mother-child
pairs).  This script renders the equivalent diagnostic for the cohort
fed to `mitoquest ne-estimate`.

For every candidate Ne in `[--min-ne, --max-ne]` (step `--ne-profile-step`)
the upstream C++ command scores two **independent** goodness-of-fit
metrics on the *same* informative pair set:

* **MLE log-likelihood** under the configured model (continuous
  Beta-diffusion or discrete Beta-Binomial).  Maximised at the fitted
  `Ne_MLE`.
* **Kimura per-pair SSR** = Σᵢ ((dᵢ − sᵢ) − p_mᵢ (1 − p_mᵢ) / Ne)².
  Minimised at the analytic
      `Ne_Kimura_SSR = Σ w² / Σ rw`,
  which is the closed-form least-squares fit of the one-generation
  Wright-Fisher prediction.

The figure has two panels (MLE on the left, Kimura on the right) so the
user can directly see whether the two estimators agree on the location
of the best Ne, or whether the data are pulling them in different
directions (a strong indicator of high-drift outliers).

##### 2. Upstream input — generate the TSV with `mitoquest ne-estimate`
```bash
mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --cross-check kimura --kimura-bootstrap 200 \
    --ne-profile cohort.ne_profile.tsv --ne-profile-step 0.1 \
    --max-ne 30 \
    -o cohort.ne.json
```
The TSV starts with `#key=value` provenance lines (fitted Ne, model,
VAF window, Wonnapinij b/Ne, Kimura bootstrap CI, …) followed by a
standard 5-column body:
```
ne_candidate  mle_log_lik  mle_delta_2ll  kimura_ssr  kimura_norm_ssr
```

##### 3. Plot it
```bash
python tools/plot_ne_profile.py \
    -i  cohort.ne_profile.tsv \
    -o  cohort.ne_profile.png
```
Key options: `--dpi` (default 150), `--figsize W,H` (default `13,5.2`),
`--title` (optional super-title).

##### 4. How to read it
* The **left (red) panel** shows the standard `−2(logL − logL_max)`
  profile-likelihood curve.  The fitted `Ne_MLE` sits at the global
  minimum (= 0); the dashed horizontal line at 3.841 is the χ² (1 df,
  0.95) threshold whose intercepts define the 95% CI bracket.
* The **right (blue) panel** shows the Kimura per-pair SSR normalised
  by its minimum.  The solid vertical line marks `Ne_Kimura_SSR`
  (least-squares best fit), the dashed line marks the
  Wonnapinij/method-of-moments `Ne_Kimura = 1 / (1 − b)`, the shaded
  band is the bootstrap 95% CI for the Wonnapinij Ne, and the green
  dotted line at Ne = 3 is the deCODE 2024 reference.
* **Both panels agree** ⇒ the data are well-described by the
  single-generation Wright-Fisher model and either estimator is fine.
* **MLE bowl is far to the left of the Kimura bowl** (e.g. Ne_MLE ≈
  2.6 vs Ne_Kimura ≈ 4.8 in the demo cohort) ⇒ the cohort contains
  high-drift outlier pairs that pull the variance-of-moments Kimura
  upward but do not hurt the MLE.  Re-run `mitoquest ne-estimate` with
  `--kimura-trim 0.10 --top-drift-k 20` to identify and inspect them.
* The **Ne = 1 grid point is dropped** automatically: under the
  continuous Beta-diffusion model Ne = 1 is a degenerate point
  (complete drift to fixation) and gives `−∞` log-likelihood for any
  pair whose child VAF is strictly inside (0, 1).

