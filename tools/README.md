## The functions and usage of scripts

#### mito_classifier.py
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