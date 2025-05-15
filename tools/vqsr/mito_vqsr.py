"""
========================================================================
= Variant quality score recalibrator (VQSR) for mitochondrial variants =
========================================================================
Date: 2025-05-15
"""
import argparse
import gzip
import re

import logging
import warnings
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_curve, auc
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('default')
import seaborn as sns

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def open_file(file_path):
    """
    Open a file, regardless of whether it is gzipped or not.
    """
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file


def vcf_info2dataframe(vcf_file):
    """
    Read vcf file and convert some information to pandas dataframe.
    """
    vcf_list = []
    colname = ["Sample_name","allGT", "snv_type", "CHROM", "POS", "REF", "ALT", "GT", "AD", "HF", "HQ", "LHF", "SB", "VT",'AF',
                'gnomad_af_hom','gnomad_af_het','in_phylotree', 'helix_af_hom','helix_af_het', 'mitomap_af','mitomap_status']
    vcf_list.append(colname)
    vcff = open_file(vcf_file)
    for line in vcff:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            _CHROM,_POS,_ID,_REF,_ALT,_QUAL,_FILTER,_INFO,_FORMAT,*SAMPLES = line.strip().split('\t')
        else:
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*samples = line.strip().split('\t')
            for sample_name, sample_data in zip(SAMPLES, samples):
                if sample_data.startswith('.'):
                    continue
                af_list = re.findall(r'AF=(.*?);', INFO)[0].split(',')
                gnomad_af_hom_list = re.findall(r'gnomad_af_hom=(.*?);', INFO)[0].split(',')
                gnomad_af_het_list = re.findall(r'gnomad_af_het=(.*?);', INFO)[0].split(',')
                in_phylotree_list = re.findall(r'in_phylotree=(.*?);', INFO)[0].split(',')
                helix_af_hom_list = re.findall(r'helix_af_hom=(.*?);', INFO)[0].split(',')
                helix_af_het_list = re.findall(r'helix_af_het=(.*?);', INFO)[0].split(',')
                mitomap_af_list = re.findall(r'mitomap_af=(.*?);', INFO)[0].split(',')
                mitomap_status_list = re.findall(r'mitomap_status=(.*?);mitomap_plasmy', INFO)[0].split(',')
                GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT = sample_data.split(':')
                GT_list = GT.split('/')
                if len(GT_list) == 1 and GT != '0':
                    snv_type = 'SNP'
                else:
                    snv_type = 'HF'
                for i in range(len(GT_list)):
                    if GT_list[i] == '0':
                        gt = '0'
                        alt = REF
                    else:
                        gt = GT_list[i]
                        alt = ALT.split(',')[int(gt)-1]
                    ad = int(AD.split(',')[i])
                    hf = float(HF.split(',')[i])
                    hq = HQ.split(',')[i]
                    lhf = float(LHF.split(',')[i])
                    sb = SB.split(';')[i]
                    vt = VT.split(',')[i]
                    af = float(af_list[int(gt)-1])
                    gnomad_af_hom = float(gnomad_af_hom_list[int(gt)-1])
                    gnomad_af_het = float(gnomad_af_het_list[int(gt)-1])
                    in_phylotree = in_phylotree_list[int(gt)-1]
                    helix_af_hom = float(helix_af_hom_list[int(gt)-1])
                    helix_af_het = float(helix_af_het_list[int(gt)-1])
                    mitomap_af = float(mitomap_af_list[int(gt)-1])
                    mitomap_status = mitomap_status_list[int(gt)-1]
                    vcf_list.append([sample_name,GT, snv_type,CHROM,POS,REF,alt,gt,ad,hf,hq,lhf,sb,vt,af,
                                        gnomad_af_hom,gnomad_af_het,in_phylotree,helix_af_hom,helix_af_het,mitomap_af,mitomap_status])
    vcfdf = pd.DataFrame(vcf_list[1:], columns=vcf_list[0])
    return vcfdf


def combine_df(df1, df2):
    """
    Combine two pandas dataframes, randomly sampling from one dataframe if necessary to ensure they have the same number of rows
    """
    if len(df1) > len(df2):
        df1_sampled = df1.sample(n=len(df2), random_state=42)
        final_df = pd.concat([df1_sampled, df2], ignore_index=True)
    else:
        df2_sampled = df2.sample(n=len(df1), random_state=42)
        final_df = pd.concat([df1, df2_sampled], ignore_index=True)
    return final_df


def clean_data(data):
    """
    Retain values within ± 6 standard deviations of the mean
    """
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    lower_bound = mean - 6 * std
    upper_bound = mean + 6 * std
    mask = np.all((data >= lower_bound) & (data <= upper_bound), axis=1)
    filtered_data = data[mask]
    return filtered_data


def GMM_fit(train_data, test_data, num_components, max_num_components):
    """
    Fit a GMM model to the data and return the predicted labels and probabilities
    """
    bics_A = []
    best_gmm_A = None
    best_bic_A = np.inf
    n_range = range(1, max_num_components+1)
    for n in n_range:
        logging.info(f"Fitting GMM with {n} components")
        gmm_A = GaussianMixture(n_components=n, random_state=0).fit(train_data)
        gmm_bic_A = gmm_A.bic(train_data)
        bics_A.append(gmm_bic_A)
        bic_A = gmm_bic_A
        if bic_A < best_bic_A:
            best_bic_A = bic_A
            best_gmm_A = gmm_A
    if num_components == '':
        scores_A = best_gmm_A.score_samples(test_data) / np.log(10)
    else:
        gmm_A = GaussianMixture(n_components=num_components, random_state=0).fit(train_data)
        best_gmm_A = gmm_A
        gmm_bic_A = gmm_A.bic(train_data)
        scores_A = gmm_A.score_samples(test_data) / np.log(10)
    return scores_A, n_range, bics_A, best_gmm_A


def find_threshold(df, colnameA, colnameB, target_ratio=0.98):
    """
    Find a threshold, if score is greater than this threshold, then the quantity of specified elements in columnA 
    will account for 0.95 of the overall elements of that type.
    
    df: pandas DataFrame
    colnameA: value type, like 1 or 0; 'Good' or 'Bad'
    colnameB: GMM score diff corresponding to each element in columnA
    target_ratio: float, ratio of ones in colnameA to total number of ones in the DataFrame
    """
    if not {colnameA, colnameB}.issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns {colnameA} and {colnameB}")
    total_ones = (df[colnameA] == 1).sum()
    if total_ones == 0:
        raise ValueError(f"No rows with {colnameA} == 1 found")
    
    target_ones = int(total_ones * target_ratio)
    sorted_df = df[[colnameA, colnameB]].sort_values(by=colnameB, ascending=False)
    
    current_ones = 0
    threshold = None
    for index, row in sorted_df.iterrows():
        if row[colnameA] == 1:
            current_ones += 1
        if current_ones >= target_ones:
            threshold = row[colnameB]
            break
    if threshold is None:
        raise ValueError("No threshold found to satisfy the condition")
    return threshold


def label_uncertain(row, colname, scores_A, threshold):
    if row[scores_A] < threshold:
        return 'Bad'
    elif row[colname] == 1:
        return 'Good'
    else:
        return 'Unlabelled'


def caculate_good_bad_ratio(df, label_col, LDD_col, good_ratio, bad_ratio ):
    """
    Calculate the ratio of good and bad elements for each score_diff value and add two new columns to the input dataframe.
    Parameters:
    -----------
    df: pandas.DataFrame
        The input dataframe containing the label_col, LDD_col columns.
    label_col: str
        The column name of labels gererated by good module, like 'Good' or 'Bad' or 'Unlabelled'.
    LDD_col: str
        The column name of LDD values, like 'LDD_score'.
    good_ratio: str
        The column name of good ratio values, like 'Good_ratio'.
    bad_ratio: str
        The column name of bad ratio values, like 'Bad_ratio'.

    Returns:
    --------
    df: pandas.DataFrame
        The dataframe containing the good_ratio, bad_ratio columns.
    """
    labels = df[label_col].values
    values = df[LDD_col].values
    total_good = (labels == "Good").sum()
    total_bad = (labels == "Bad").sum()
    sort_idx = np.argsort(-values) 
    sorted_labels = labels[sort_idx]
    sorted_values = values[sort_idx]
    cum_good = np.cumsum(sorted_labels == "Good")
    cum_bad = np.cumsum(sorted_labels == "Bad")
    good_ratios = np.zeros(len(df))
    bad_ratios = np.zeros(len(df))
    for i in range(len(df)):
        good_count = cum_good[i] if i > 0 else 0
        bad_count = cum_bad[i] if i > 0 else 0
        good_ratios[sort_idx[i]] = good_count / total_good if total_good > 0 else 0
        bad_ratios[sort_idx[i]] = bad_count / total_bad if total_bad > 0 else 0
    
    df[good_ratio] = good_ratios
    df[bad_ratio] = bad_ratios
    
    return df


def plt_gmm_results(n_range, bics_A, best_gmm_A , bics_B, best_gmm_B, cutoff_dict,
                    df_sorted, output_figure):
    fig, ax = plt.subplots(figsize=(16, 4), nrows=1, ncols=3, constrained_layout=True)
    sns.set_style('ticks')
    sns.scatterplot(x= n_range, y= bics_A, ax=ax[0])
    sns.lineplot(x=n_range, y=bics_A, ax=ax[0])
    ax[0].set(xlabel='Number of Components', ylabel='BIC', title='Good Models: BIC')
    ax[0].scatter(n_range[best_gmm_A.n_components-1],bics_A[best_gmm_A.n_components-1],color='red', 
                marker='x', s = 100, label='Optimal Number of Components')
    ax[0].legend(loc='right')

    sns.scatterplot(x= n_range, y= bics_B, ax=ax[1])
    sns.lineplot(x=n_range, y=bics_B, ax=ax[1])
    ax[1].set(xlabel='Number of Components', ylabel='BIC', title='Bad Model: BIC')
    ax[1].scatter(n_range[best_gmm_B.n_components-1],bics_B[best_gmm_B.n_components-1],color='red', 
                marker='x', s = 100, label='Optimal Number of Components')
    ax[1].legend(loc='right')
    

    sns.scatterplot(data=df_sorted, x="Bad_ratio", y="Good_ratio", label="different score_diff", s=5, alpha=0.4, ax=ax[2])
    colors = [ 'blue','green','red', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'lime', 'pink']
    for i, cutoff_idx in enumerate(cutoff_dict):
        ax[2].scatter(
            cutoff_dict[cutoff_idx][1],
            cutoff_dict[cutoff_idx][0],
            color=colors[i],
            s=100,
            marker="x",
            label=f"LODR={cutoff_dict[cutoff_idx][2]:.4f}\nGood ratio={cutoff_dict[cutoff_idx][0]:.2%}, Bad ratio={cutoff_dict[cutoff_idx][1]:.2%}")

    ax[2].set(xlabel='Bad value ratio', ylabel='Good value ratio', title='Frequency of Good/Bad for per LODR', ylim=(0.94, 1.01))
    ax[2].legend(loc='right')
    plt.savefig(output_figure)


def output_info(df, colname, train_A_data, train_B_data_filtered, scores_A_cutoff, cutoff_dict, optimal_threshold_idx, logfilename):
    """
    Output the log information.
    """
    # Initial classification
    raw_good_count = (df[colname] == 1).sum()
    raw_good_POS = df[df[colname] == 1]['POS'].unique().__len__()
    HQ_alt_count = df[(df['GT'] !='0') & (df[colname] == 1)].__len__()
    HQ_ref_count = df[(df['GT'] =='0') & (df[colname] == 1)].__len__()
    HQ_alt_POS = df[(df['GT'] !='0') & (df[colname] == 1)]['POS'].unique().__len__()
    HQ_ref_POS = df[(df['GT'] =='0') & (df[colname] == 1)]['POS'].unique().__len__()
    
    selected_good_count = df[df['good_train_value'] == 1].__len__()
    selected_good_POS = df[df['good_train_value'] == 1]['POS'].unique().__len__()
    # Retain values within ± 6 standard deviations of the mean for training Good module data
    filtered_good_count = train_A_data.shape[0]

    # Num of greater than threshold / num of good data = 0.95
    # Classify all dataset to Good, Uncertain, or Bad using scores-good-cutoff
    Good_count = (df['Good_module_res']=='Good').sum()
    Good_POS = df[df['Good_module_res'] == 'Good']['POS'].unique().size
    Uncertain_count = (df['Good_module_res']=='Unlabelled').sum()
    Uncertain_POS = df[df['Good_module_res'] == 'Unlabelled']['POS'].unique().size
    Bad_count = (df['Good_module_res']=='Bad').sum()
    Bad_POS = df[df['Good_module_res'] == 'Bad']['POS'].unique().size

    filtered_Bad_count = train_B_data_filtered.shape[0]
    with open(logfilename, 'w') as f:
        print("############### datasets information ###############", file=f)
        print (f"Total count:\t{df.shape[0]}\t unique POS: {df['POS'].unique().size}", file=f)
        print(f"Good count:\t{raw_good_count}\tunique POS: {raw_good_POS}", file=f)
        print(f"Good alt count:\t{HQ_alt_count}\tunique POS: {HQ_alt_POS}", file=f)
        print(f"Good ref count:\t{HQ_ref_count}\tunique POS: {HQ_ref_POS}", file=f)
        
        print("\n############### information of Good module ###############", file=f)
        print(f"Good training data:\t{selected_good_count}\tunique POS: {selected_good_POS}" , file=f)
        # print(f"Raw Uncertain count:\t {raw_uncertain_count}\tunique POS: {raw_uncertain_POS}", file=f)
        print(f"Filtered Good training data:\t{filtered_good_count}\t", file=f)
        
        print("\n############### information of Bad module ###############", file=f)
        print(f"Score cutoff to classify data to good, unlabelled, or bad:\t{scores_A_cutoff}", file=f)
        print(f'Good count:\t{Good_count}\tunique POS: {Good_POS}', file=f)
        print(f"Unlabelled count:\t{Uncertain_count}\tunique POS: {Uncertain_POS}", file=f)
        print(f"Bad training data:\t{Bad_count}\tunique POS: {Bad_POS}", file=f)
        print(f"Filtered Bad training data:\t{filtered_Bad_count}", file=f)
        
        print("\n############### Choose optimal threshold ###############", file=f)
        for idx, (good_ratio, bad_ratio, threshold) in cutoff_dict.items():
            print(f"Good ratio: {good_ratio}\tBad ratio: {bad_ratio}\tThreshold: {threshold}", file=f)
        print(f"Optimal threshold:\tGood ratio: {df.loc[optimal_threshold_idx, 'Good_ratio']}\tBad ratio: {df.loc[optimal_threshold_idx, 'Bad_ratio']}\tThreshold: {df.loc[optimal_threshold_idx, 'LODR']}", file=f)
        
        print("\n############### GMM module result (LowQual POS has not been removed) ###############", file=f)
        good_count = (df['GMM_res']== 'Good').sum()
        good_pos = df[df['GMM_res'] == 'Good']['POS'].unique().size
        bad_count = (df['GMM_res']== 'Bad').sum()
        bad_pos = df[df['GMM_res'] == 'Bad']['POS'].unique().size
        print (f"Final result :\nGood count:{good_count}\t unique POS: {good_pos}\nBad count:{bad_count}\t unique POS: {bad_pos}", file=f)
        
        print("\n############### Other infomation ###############", file=f)
        _value_count = df[(df[colname]==1) & (df['GMM_res'] == 'Bad')].__len__()
        _pos_count = df[(df[colname]==1) & (df['GMM_res'] == 'Bad')]['POS'].unique().__len__()
        print(f"Good training data in final bad result:\t{_value_count}\tunique POS\t{_pos_count}", file=f)


def df_info2dict(final_df):
    """
    Convert the information in the final dataframe to dictionaries
    param: final_df: the final dataframe containing the final predicted results
    return: dictionaries containing the information of the final dataframe
    """
    # make a dictionary of the good count for each POS
    tmp_df = final_df.groupby(['POS'])['GMM_res'].agg('value_counts').reset_index(name='count')
    POS_good_count_dict = defaultdict(int)
    for idx, row in tmp_df.iterrows():
        if row['GMM_res'] == 'Good':
            POS_good_count_dict[row['POS']] = int(row['count'])
    
    # make a dictionary for each POS, Sample_name, and GT, store the LODR and GMM_res
    VQSR_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for idx, row in final_df.iterrows():
        VQSR_dict[row['POS']][row['Sample_name']][row['GT']].extend([round(float(row['LODR']), 4), row['GMM_res']])
    return POS_good_count_dict, VQSR_dict


def main():
    parser = argparse.ArgumentParser(description='Build GMM model using AD, HF, HQ information of mtDNA data to filter out bad sites values.')
    parser.add_argument('-i','--input_vcf_file', type=str, help='Input VCF file path', required=True)
    parser.add_argument('-o', '--output_vcf_file', type=str, help='Output VCF file path', required=True)
    parser.add_argument('-c','--output_csv_file', type=str, help='Output CSV file path', required=True)
    parser.add_argument('-f','--figure_file', type=str, help='Output figure file path', required=True)
    parser.add_argument('-mnc','--max_n_components', type=int, action='store',default=10, help='Maximum components for GMM model training (default: 10)')
    parser.add_argument('-gnc', '--good_module_n_components', type=int, action='store', help='Number of components for Good GMM module if not specified, default is auto-selected')
    parser.add_argument('-bnc', '--bad_module_n_components', type=int, action='store',  help='Number of components for Bad GMM moduleif not specified, default is auto-selected')
    parser.add_argument('-rgm', '--goodratio_in_goodmodule', type=float, action='store', default=0.98, help='Ratio of good values to total good values after Good module to select subdataset (Bad values) to train Bad GMM model (default: 0.98)')
    parser.add_argument('-pr', '--pass_ratio', type=float, action='store', default=0.5, help='Ratio of Good Allels to total Allels per POS to label a site as PASS or LowQual (default: 0.5)')
    parser.add_argument('-fmr', '--final_res_ratio', type=float, action='store', default=0.99, help='The Good ratio of the final result is used to determine whether the value is good or bad (default: 0.99)')
    args = parser.parse_args()
    
    # load data
    logging.info('Loading data ...')
    dataset = vcf_info2dataframe(args.input_vcf_file)
    dataset['gnomad_af'] = dataset['gnomad_af_het'] + dataset['gnomad_af_hom']
    dataset['helix_af'] = dataset['helix_af_het'] + dataset['helix_af_hom']
    dataset[['AF', 'gnomad_af', 'helix_af', 'mitomap_af', 'AD', 'HF', 'HQ']] = dataset[['AF', 'gnomad_af', 'helix_af', 'mitomap_af', 'AD', 'HF', 'HQ']].apply(pd.to_numeric, errors='coerce')
    dataset['Sample_name_POS_GT'] = dataset['Sample_name'] + '_' + dataset['POS'] + '_' + dataset['GT']

    # select sites with good quality: 1 means good quality, 0 means bad quality
    logging.info('Select high good quality datasets, like reported or AF>0.01 in other database ...')
    HQcolname = 'Reported_AFgt005'
    dataset[HQcolname] = dataset.apply(
    lambda row: 1 if ((row['mitomap_status'] != '') or 
                        row['gnomad_af']>=0.01 or 
                        row['helix_af']>=0.01 or
                        row['mitomap_af']>=0.01 or 
                        row['AF']>=0.05) else 0,
                        axis=1)
    
    logging.info('Z-score normalize data and select good training data')
    dataset[['ADz', 'HFz', 'HQz']] = dataset[['AD', 'HF', 'HQ']].apply(zscore)
    all_data_scaled = np.array(dataset[['ADz', 'HFz', 'HQz']])
    
    HQ_dataset = dataset[dataset[HQcolname] == 1]
    # LQ_dataset = dataset[dataset[HQcolname] == 0]
    HQ_dataset_alt = HQ_dataset[HQ_dataset['GT'] !='0']
    HQ_dataset_ref = HQ_dataset[HQ_dataset['GT'] =='0']
    tmp_good_train_data = combine_df(HQ_dataset_alt, HQ_dataset_ref)
    good_train_value_list = tmp_good_train_data["Sample_name_POS_GT"].tolist()
    dataset['good_train_value'] = dataset['Sample_name_POS_GT'].isin(good_train_value_list).astype(int)
    good_train_data = np.array(tmp_good_train_data[['ADz', 'HFz', 'HQz']])

    logging.info('Retain values of good training data within ± 6 standard deviations of the mean')
    filtered_good_train_data = clean_data(good_train_data)

    logging.info('Build Good GMM model using filtered good training data and evaluate on test data')
    max_n_components = args.max_n_components
    if args.good_module_n_components:
        good_module_n_components = args.good_module_n_components
    else:
        good_module_n_components = ''
    scores_Good, n_range, bics_Good, best_gmm_Good = GMM_fit(filtered_good_train_data, all_data_scaled, good_module_n_components, max_n_components)
    dataset['scores_Good'] = scores_Good
    
    logging.info('Select the Bad training data set to build the Bad GMM model and evaluate on test data')
    scores_good_cutoff = find_threshold(dataset, HQcolname, 'scores_Good', target_ratio=args.goodratio_in_goodmodule)
    dataset['Good_module_res'] = dataset.apply(lambda row: label_uncertain(row, HQcolname, 'scores_Good', scores_good_cutoff), axis=1)
    Bad_train_data = np.array(dataset[(dataset['Good_module_res']=='Bad')][['ADz', 'HFz', 'HQz']])
    Bad_train_data_filtered = clean_data(Bad_train_data)
    if args.bad_module_n_components:
        bad_module_n_components = args.bad_module_n_components
    else:
        bad_module_n_components = ''
    scores_Bad, n_range, bics_Bad, best_gmm_Bad = GMM_fit(Bad_train_data_filtered, all_data_scaled, bad_module_n_components, max_n_components)
    score_diff = scores_Good - scores_Bad
    dataset['scores_Bad'] = scores_Bad
    dataset['LODR'] = score_diff

    logging.info('Calculate the ratio of Good and Bad values for each LODR')
    dataset = caculate_good_bad_ratio(dataset, 'Good_module_res', 'LODR', "Good_ratio", "Bad_ratio")
    
    logging.info('Specify the Good ratio (0.95, 0.98, 0.99) corresponding to the final threshold')
    closest_95_idx = (dataset['Good_ratio'] - 0.95).abs().idxmin()
    closest_95_Good_ratio = dataset.loc[closest_95_idx, 'Good_ratio']
    closest_95_Bad_ratio = dataset.loc[closest_95_idx, 'Bad_ratio']
    closest_95_LODR = dataset.loc[closest_95_idx, 'LODR']

    closest_98_idx = (dataset['Good_ratio'] - 0.98).abs().idxmin()
    closest_98_Good_ratio = dataset.loc[closest_98_idx, 'Good_ratio']
    closest_98_Bad_ratio = dataset.loc[closest_98_idx, 'Bad_ratio']
    closest_98_LODR = dataset.loc[closest_98_idx, 'LODR']

    closest_99_idx = (dataset['Good_ratio'] - 0.99).abs().idxmin()
    closest_99_Good_ratio = dataset.loc[closest_99_idx, 'Good_ratio']
    closest_99_Bad_ratio = dataset.loc[closest_99_idx, 'Bad_ratio']
    closest_99_LODR = dataset.loc[closest_99_idx, 'LODR']
    
    cutoff_dict = {closest_95_idx:[closest_95_Good_ratio, closest_95_Bad_ratio, closest_95_LODR],
                closest_98_idx:[closest_98_Good_ratio, closest_98_Bad_ratio, closest_98_LODR],
                closest_99_idx:[closest_99_Good_ratio, closest_99_Bad_ratio, closest_99_LODR]}

    optimal_threshold_idx = (dataset['Good_ratio'] - args.final_res_ratio).abs().idxmin()
    optimal_threshold = dataset.loc[optimal_threshold_idx, 'LODR']
        
    logging.info('Output the final result using GMM model')
    dataset['GMM_res'] = dataset.apply(lambda x: 'Good' if x['LODR'] > optimal_threshold else 'Bad', axis=1)
    
    logging.info('Plot the results....')
    plt_gmm_results(n_range, bics_Good, best_gmm_Good, bics_Bad,best_gmm_Bad,
                    cutoff_dict,dataset, args.figure_file)
    
    logging.info('Output the log information....')
    output_info(dataset, HQcolname, filtered_good_train_data, Bad_train_data_filtered, scores_good_cutoff,
                cutoff_dict, optimal_threshold_idx, args.output_csv_file.strip('.csv')+'.log')
        
    logging.info('Output the final result to VCF file....')
    POS_good_count_dict, VQSR_dict = df_info2dict(dataset)
    
    good_allele_ratio = args.pass_ratio
    output_vcf_file = args.output_vcf_file.strip('.gz')
    POS_count = 0
    POS_pass_count = 0
    POS_lowqual = []
    with open(output_vcf_file, 'w') as output_vcf:
        for line in open_file(args.input_vcf_file):
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    output_vcf.write(f'##FORMAT=<ID=LODR,Number=R,Type=Float,Description="An ordered, comma delimited Log-Likelihood Difference for non-reference allele alleles in the order listed">\n')
                    # output_vcf.write(f'##FORMAT=<ID=FILTER,Number=R,Type=String,Description="An ordered, comma delimited list of non-reference allele filtration situation ">\n')
                    output_vcf.write(f'##mito_classifier_command=python GMM_classifier.py -i {args.input_vcf_file} -o {args.output_vcf_file} -c {args.output_csv_file} -f {args.figure_file}\n')
                    _CHROM,_POS,_ID,_REF,_ALT,_QUAL,_FILTER,_INFO,_FORMAT,*_SAMPLES = line.strip().split('\t')
                output_vcf.write(line)
            else:
                POS_count += 1
                CHROM,POS,ID,REFs,ALTs,QUAL,FILTER,INFO,FORMAT, *SAMPLES = line.strip().split('\t')
                FORMAT = 'GT:GQ:DP:AD:HF:CI:HQ:LHF:SB:FS:SOR:VT:LODR'
                AN = re.findall(r'AN=(.*?);', INFO)[0].split(',')[0]
                if POS_good_count_dict[POS] / int(AN) >= good_allele_ratio:
                    FILTER = 'PASS'
                    POS_pass_count += 1
                else:
                    FILTER = 'LowQual'
                    POS_lowqual.append(POS)
                newSAMPLES = []
                for sample_name, sample_data in zip(_SAMPLES, SAMPLES):
                    if sample_data.startswith('.'):
                        continue
                    GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT = sample_data.split(':')
                    LODR_list = []
                    GT_list = []
                    for gt in GT.split('/'):
                        LODR_list.append(str(VQSR_dict[POS][sample_name][gt][0]))
                        if VQSR_dict[POS][sample_name][gt][1] == 'Bad':
                            GT_list.append('.')
                        else:
                            GT_list.append(gt)
                    GT = '/'.join(GT_list)
                    LODR = ','.join(LODR_list)
                    new_sample_data = f'{GT}:{GQ}:{DP}:{AD}:{HF}:{CI}:{HQ}:{LHF}:{SB}:{FS}:{SOR}:{VT}:{LODR}'
                    newSAMPLES.append(new_sample_data)
                print(CHROM,POS,ID,REFs,ALTs,QUAL,FILTER,INFO,FORMAT, *newSAMPLES, sep='\t', file=output_vcf)
                
    logging.info('Output the final result to CSV file....')
    dataset.drop(columns=['Reported_AFgt005','ADz', 'HFz', 'HQz', 'scores_Good', 'Good_module_res', 'scores_Bad',
                        'Good_ratio', 'Bad_ratio',], inplace=True)
    with open(args.output_csv_file.strip('.csv')+'.log', 'a') as f:
        print("\n############### POS statistics after GMM Filtering ###############", file=f)
        print("Total number of positions: ", POS_count, file=f)
        print("Number of positions PASS: ", POS_pass_count, file=f)
        print("Ratio of PASS positions to total positions: ", POS_pass_count/POS_count, file=f)
        
        print("\n############### GMM module result (LowQual POS has been removed) ###############", file=f)
        filtered_df = dataset[~dataset['POS'].isin(POS_lowqual)]
        good_count = (filtered_df['GMM_res']== 'Good').sum()
        good_pos = filtered_df[filtered_df['GMM_res'] == 'Good']['POS'].unique().size
        bad_count = (filtered_df['GMM_res']== 'Bad').sum()
        bad_pos = filtered_df[filtered_df['GMM_res'] == 'Bad']['POS'].unique().size
        print (f"Final result :\nGood count:{good_count}\t unique POS: {good_pos}\nBad count:{bad_count}\t unique POS: {bad_pos}", file=f)
    filtered_df.to_csv(args.output_csv_file, index=False)
    
    logging.info('Compress and index the VCF file....')
    subprocess.run(["bgzip", "-f", output_vcf_file])
    subprocess.run(["tabix", "-p", "vcf", output_vcf_file+".gz"])
    logging.info('Done!')


if __name__ == '__main__':
    main()
