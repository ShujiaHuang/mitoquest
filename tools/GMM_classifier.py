import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from scipy.stats import zscore
from scipy.signal import savgol_filter
import gzip
import warnings
import re
import subprocess
from collections import defaultdict
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

def open_file(file_path):
    """
    Open a file, regardless of whether it is gzipped or not.
    """
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file

def vcf_info2dataframe(vcf_file, colname = 'all'):
    """
    Read vcf file and convert some information to pandas dataframe.
    """
    vcf_list = []
    if colname == 'all':
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
                        continue
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

def GMM_fit(train_data, test_data, num_components='', max_num_components=15):
    """
    Fit a GMM model to the data and return the predicted labels and probabilities
    """
    bics_A = []
    best_gmm_A = None
    best_bic_A = np.inf
    n_range = range(1, max_num_components)
    for n in n_range:
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

def find_threshold(df, colnameA, colnameB, target_ratio=0.95):
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
        return 'Uncertain'

def caculate_good_bad_ratio(df, label_col, LDD_col, good_ratio, bad_ratio ):
    """
    Calculate the ratio of good and bad elements for each score_diff value and add two new columns to the input dataframe.
    Parameters:
    -----------
    df: pandas.DataFrame
        The input dataframe containing the label_col, LDD_col columns.
    label_col: str
        The column name of labels gererated by good module, like 'Good' or 'Bad'.
    LDD_col: str
        The column name of LDD values, like 'LDD_score'.
    good_ratio: str
        The column name of good ratio values, like 'Good_ratio'.
    bad_ratio: str
        The column name of bad ratio values, like 'Bad_ratio'.

    Returns:
    --------
    df: pandas.DataFrame
        The input dataframe containing the good_ratio, bad_ratio columns.
    """
    good_total = (df[label_col] == "Good").sum()
    bad_total = (df[label_col] == "Bad").sum()
    for i, row in df.iterrows():
        current_score = row[LDD_col]
        subset1 = df[df[LDD_col] > current_score]
        # calculate good_freq1 and bad_freq1 for each score_diff value
        if len(subset1) > 0:
            good_count1 = (subset1[label_col] == "Good").sum()
            bad_count1 = (subset1[label_col] == "Bad").sum()
            df.at[i, good_ratio] = good_count1 / good_total
            df.at[i, bad_ratio] = bad_count1 / bad_total
        else:
            df.at[i, good_ratio] = 0  
            df.at[i, bad_ratio] = 0
    return df

def calcu_inflection_points(df, x_col, y_col, value_col, raw_good_col, raw_bad_col):
    """
    Find the best inflection point for Good/Bad classification based on GMM classifier.
    Parameters:
    -----------
    df: pandas.DataFrame
        The input dataframe containing the x_col, y_col, value_col, raw_good_col, raw_bad_col columns.
    x_col: str
        The column name of x-axis values, like Bad ratio.
    y_col: str
        The column name of y-axis values, like Good ratio.
    value_col: str
        The column name of values to be classified, like LLD values.
    raw_good_col: str
        The column including the raw good values that used to build the good model, like HQcolname 
    raw_bad_col: str
        The column including the raw bad values that used to build the bad model, like Good_module_res

    Returns:
    --------
    best_inflection: float
        The best inflection point for Good/Bad classification.
    inflection_points: pandas.DataFrame
        The dataframe containing the inflection points and their corresponding values.
    log_info: list
        The list containing the log information of the inflection point selection process.
    """
    df_sorted = df.sort_values(x_col).reset_index(drop=True)
    df_sorted["smoothed_y"] = savgol_filter(df_sorted[y_col], window_length=11, polyorder=1)

    # # 2. caculate slope and slope_change
    df_sorted["slope"] = np.gradient(df_sorted["smoothed_y"], df_sorted[x_col])
    df_sorted["slope_change"] = np.abs(np.gradient(df_sorted["slope"]))

    # 3. extract inflection points
    inflection_points = df_sorted.nlargest(3, "slope_change")
    inflection_points_list = []
    for idx, row in inflection_points.iterrows():
        inflection_points_list.append(row[value_col])
    _tmp_list = []
    log_info = []
    for cutoff in sorted(inflection_points_list):
        df_sorted['predicted'] = df_sorted.apply(lambda x: 'Good' if x[value_col] > cutoff else 'Bad', axis=1)
        good2bad_count = df_sorted[(df_sorted["predicted"] == 'Bad') & (df_sorted[raw_good_col] == 1)].__len__()
        good_total = df_sorted[df_sorted[raw_good_col] == 1].shape[0]
        log_info.append(f"The ratio of Good to Bad under {cutoff:.4f} cutoff: {(good2bad_count/good_total):.2%}")
        bad2good_count = df_sorted[(df_sorted["predicted"] == 'Good') & (df_sorted[raw_bad_col] == 'Bad')].__len__()
        bad_total = (df_sorted[raw_bad_col] == 'Bad').sum()
        log_info.append(f"The ratio of Bad to Good under {cutoff:.4f} cutoff: {(bad2good_count/bad_total):.2%}")
        _tmp_list.append([cutoff, (good2bad_count/good_total), (bad2good_count/bad_total)])

    best_inflection = None
    best_ratio_diff = np.inf
    for cutoff, good2bad_ratio, bad2good_ratio in _tmp_list:
        ratio_diff = abs(good2bad_ratio - bad2good_ratio)
        if ratio_diff < best_ratio_diff:
            best_inflection = cutoff
    return best_inflection, inflection_points, log_info

def calcu_optimal_threshold(df, cutoff, value_col, Goodmodule_res_col):
    """
    To calculate the optimal threshold for the GMM classifier based on the ROC curve.
    :param df: the dataframe containing the predicted scores and true labels
    :param cutoff: the cutoff value for the GMM classifier, default: the best inflection point
    :param value_col: the column name of the predicted scores, default: 'LLD'
    :return: the optimal threshold, ROC AUC, and the ROC curve data
    """
    df['predicted'] = df.apply(lambda x: 'Good' if x[value_col] > cutoff else 'Bad', axis=1)
    def inflection_res(df):
        if df[Goodmodule_res_col] == "Good" and df['predicted'] == "Good":
            return 1
        elif df[Goodmodule_res_col] == "Good" and df['predicted'] == "Bad":
            return 1
        elif df[Goodmodule_res_col] == "Bad" and df['predicted'] == "Good":
            return 0
        elif df[Goodmodule_res_col] == "Bad" and df['predicted'] == "Bad":
            return 0
        elif df[Goodmodule_res_col] == "Uncertain" and df['predicted'] == "Good":
            return 1
        elif df[Goodmodule_res_col] == "Uncertain" and df['predicted'] == "Bad":
            return 0

    df['inflection_res'] = df.apply(inflection_res, axis=1)
    fpr, tpr, thresholds = roc_curve(df['inflection_res'], df['LLD'])
    roc_auc = auc(fpr, tpr)
    youden_index = tpr - fpr
    optimal_idx = np.argmax(youden_index)
    optimal_threshold = thresholds[optimal_idx]
    return optimal_threshold, roc_auc, fpr, tpr, df

def plt_gmm_results(n_range, bics_A,best_gmm_A , bics_B, best_gmm_B, 
                    df_sorted, inflection_points, fpr, tpr, roc_auc,optimal_threshold, output_figure):
    fig, ax = plt.subplots(figsize=(10, 8), nrows=2, ncols=2, constrained_layout=True)
    sns.set_style('ticks')
    sns.scatterplot(x= n_range, y= bics_A, ax=ax[0,0])
    sns.lineplot(x=n_range, y=bics_A, ax=ax[0,0])
    ax[0,0].set(xlabel='Number of Components', ylabel='BIC', title='Good Models: BIC')
    ax[0,0].scatter(n_range[best_gmm_A.n_components-1],bics_A[best_gmm_A.n_components-1],color='red', 
                marker='x', s = 100, label='Optimal Number of Components')
    ax[0,0].legend(loc='right')

    sns.scatterplot(x= n_range, y= bics_B, ax=ax[0,1])
    sns.lineplot(x=n_range, y=bics_B, ax=ax[0,1])
    ax[0,1].set(xlabel='Number of Components', ylabel='BIC', title='Bad Model: BIC')
    ax[0,1].scatter(n_range[best_gmm_B.n_components-1],bics_B[best_gmm_B.n_components-1],color='red', 
                marker='x', s = 100, label='Optimal Number of Components')
    ax[0,1].legend(loc='right')

    sns.scatterplot(data=df_sorted, x="Bad_ratio", y="Good_ratio", label="different score_diff", s=5, alpha=0.4, ax=ax[1,0])
    colors = ['green', 'blue', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'lime', 'pink']
    for i, (idx, row) in enumerate(inflection_points.iterrows()):
        ax[1,0].scatter(
            row["Bad_ratio"],
            row["Good_ratio"],
            color=colors[i],
            s=100,
            marker="x",
            label=f"Inflection {i+1} (LLD={row['LLD']:.4f})"
        )
    ax[1,0].scatter(df_sorted[(df_sorted['LLD'] == optimal_threshold)]['Bad_ratio'].values[0], 
                    df_sorted[(df_sorted['LLD'] == optimal_threshold)]['Good_ratio'].values[0], 
                    color='red', s=100, marker="x",label=f"Optimal threshold (LLD={optimal_threshold:.4f})")
    ax[1,0].set(xlabel='Bad value ratio', ylabel='Good value ratio', title='Frequency of Good/Bad Labels for larger than LLD', ylim=(0.96, 1.01))
    ax[1,0].legend(loc='right')

    sns.lineplot(x=fpr, y=tpr, color='darkorange', linewidth=2, label=f'ROC curve (AUC = {roc_auc:.4f})', ax=ax[1,1])
    sns.lineplot(x=[0, 1], y=[0, 1], color='navy', linewidth=2, linestyle='--',  ax=ax[1,1])
    ax[1,1].set(xlabel='False Positive Rate', ylabel='True Positive Rate', title='ROC Curve', ylim=(0.0, 1.01), xlim=(-0.01, 1.0))
    ax[1,1].legend(loc='right')
    plt.savefig(output_figure)
    
def output_info(df, colname, train_A_data, train_B_data_filtered, scores_A_cutoff, 
                best_inflection, inflection_res_info, optimal_threshold, logfilename):
    """
    Output the log information.
    """
    # Initial classification
    raw_good_count = (df[colname] == 1).sum()
    raw_uncertain_count = (df[colname] == 0).sum()

    raw_good_POS = df[df[colname] == 1]['POS'].unique().size
    raw_uncertain_POS = df[df[colname] == 0]['POS'].unique().size

    # Retain values within ± 6 standard deviations of the mean for training Good module data
    filtered_good_count = train_A_data.shape[0]

    # Num of greater than threshold / num of good data = 0.95
    # Classify all dataset to Good, Uncertain, or Bad using scores-good-cutoff
    Good_count = (df['Good_module_res']=='Good').sum()
    Good_POS = df[df['Good_module_res'] == 'Good']['POS'].unique().size
    Uncertain_count = (df['Good_module_res']=='Uncertain').sum()
    Uncertain_POS = df[df['Good_module_res'] == 'Uncertain']['POS'].unique().size
    Bad_count = (df['Good_module_res']=='Bad').sum()
    Bad_POS = df[df['Good_module_res'] == 'Bad']['POS'].unique().size

    filtered_Bad_count = train_B_data_filtered.shape[0]
    with open(logfilename, 'w') as f:
        print("############### datasets information ###############", file=f)
        print(f"Raw Good count:\t{raw_good_count}\tunique POS: {raw_good_POS}" , file=f)
        print(f"Raw Uncertain count:\t {raw_uncertain_count}\tunique POS: {raw_uncertain_POS}", file=f)
        print(f"Filtered Good count:\t{filtered_good_count}", file=f)
        print(f"Score cutoff to classify data to good, uncertain, or bad:\t{scores_A_cutoff}", file=f)
        print(f'Good count:\t{Good_count}\tunique POS: {Good_POS}', file=f)
        print(f"Uncertain count:\t{Uncertain_count}\tunique POS: {Uncertain_POS}", file=f)
        print(f"Bad count:\t{Bad_count}\tunique POS: {Bad_POS}", file=f)
        print(f"Filtered Bad count:\t{filtered_Bad_count}", file=f)
        
        print("############### optimal inflection point caculation ###############", file=f)
        good_ratio1 = df[df['LLD'] == best_inflection]['Good_ratio'].values[0]
        bad_ration2 = df[df['LLD'] == best_inflection]['Bad_ratio'].values[0]
        print(*inflection_res_info, sep='\n', file=f)
        print(f"optimal inflection point:\tscore_diff:{best_inflection}\tgood_ratio:{good_ratio1:.2%}\tbad_ratio:{bad_ration2:.2%}", file=f)
        
        print("############### optimal thresholdcaculation ###############", file=f)
        good_ratio1 = df[df['LLD'] == optimal_threshold]['Good_ratio'].values[0]
        bad_ration2 = df[df['LLD'] == optimal_threshold]['Bad_ratio'].values[0]
        print(f"optimal threshold:\t{optimal_threshold}\tgood_ratio:{good_ratio1:.2%}\tbad_ratio:{bad_ration2:.2%}", file=f)
        
        print("############### GMM module result ###############", file=f)
        good_count = (df['GMM_res']== 'Good').sum()
        good_pos = df[df['GMM_res'] == 'Good']['POS'].unique().size
        bad_count = (df['GMM_res']== 'Bad').sum()
        bad_pos = df[df['GMM_res'] == 'Bad']['POS'].unique().size
        print (f"Final result:\nGood count:{good_count}\t unique POS: {good_pos}\nBad count:{bad_count}\t unique POS: {bad_pos}", file=f)

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
    
    # make a dictionary for each POS, Sample_name, and GT, store the LLD and GMM_res
    VQSR_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for idx, row in final_df.iterrows():
        VQSR_dict[row['POS']][row['Sample_name']][row['GT']].extend([round(float(row['LLD']), 4), row['GMM_res']])
    return POS_good_count_dict, VQSR_dict

def main():
    parser = argparse.ArgumentParser(description='Build GMM model using AD, HF, HQ information of mtDNA data to filter out bad sites values.')
    parser.add_argument('-i','--input_file', type=str, help='Input VCF file path', required=True)
    parser.add_argument('-o', '--output_vcf_file', type=str, help='Output VCF file path', required=True)
    parser.add_argument('-c','--output_csv_file', type=str, help='Output CSV file path', required=True)
    parser.add_argument('-f','--figure_file', type=str, help='Output figure file path', required=True)
    parser.add_argument('-mnc','--max_n_components', type=int, action='store',default=15, help='Maximum components for GMM model training')
    parser.add_argument('-gnc', '--good_module_n_components', type=int, action='store', help='Number of components for Good GMM module if not specified, default is auto-selected')
    parser.add_argument('-bnc', '--bad_module_n_components', type=int, action='store',  help='Number of components for Bad GMM moduleif not specified, default is auto-selected')
    parser.add_argument('-r', '--good_ratio', type=float, action='store', default=0.95, help='Ratio of good values to total good values after Good module to select subdataset (Bad values) to train Bad GMM model')
    args = parser.parse_args()
    
    # load data
    dataset = vcf_info2dataframe(args.input_file)
    dataset['gnomad_af'] = dataset['gnomad_af_het'] + dataset['gnomad_af_hom']
    dataset['helix_af'] = dataset['helix_af_het'] + dataset['helix_af_hom']
    dataset[['AF', 'gnomad_af', 'helix_af', 'mitomap_af', 'AD', 'HF', 'HQ']] = dataset[['AF', 'gnomad_af', 'helix_af', 'mitomap_af', 'AD', 'HF', 'HQ']].apply(pd.to_numeric, errors='coerce')

    # select sites with good quality: 1 means good quality, 0 means bad quality
    HQcolname = 'Reported_AFgt005'
    dataset[HQcolname] = dataset.apply(
    lambda row: 1 if ((row['mitomap_status'] != '') or 
                        row['gnomad_af']>=0.05 or 
                        row['helix_af']>=0.05 or
                        row['mitomap_af']>=0.05 or 
                        row['AF']>=0.05) else 0,
                        axis=1)
    
    # z-score normalize data and split into Good train and test sets
    dataset[['ADz', 'HFz', 'HQz']] = dataset[['AD', 'HF', 'HQ']].apply(zscore)
    all_data_scaled = np.array(dataset[['ADz', 'HFz', 'HQz']])
    good_train_data = np.array(dataset[dataset[HQcolname] == 1][['ADz', 'HFz', 'HQz']])

    # retain values of good training data within ± 6 standard deviations of the mean
    filtered_good_train_data = clean_data(good_train_data)

    # Build Good GMM model using filtered_good_train_data and evaluate on test data
    if args.max_n_components:
        max_n_components = args.max_n_components
    else:
        max_n_components = 15
    if args.good_module_n_components:
        good_module_n_components = args.good_module_n_components
    else:
        good_module_n_components = ''
    if args.bad_module_n_components:
        bad_module_n_components = args.bad_module_n_components
    else:
        bad_module_n_components = ''
    scores_Good, n_range, bics_Good, best_gmm_Good = GMM_fit(filtered_good_train_data, all_data_scaled, good_module_n_components, max_n_components)
    dataset['scores_Good'] = scores_Good
    # select the Bad training data set to build the Bad GMM model and evaluate on test data
    if args.good_ratio:
        scores_good_cutoff = find_threshold(dataset, HQcolname, 'scores_Good', target_ratio=args.good_ratio)
    else:
        scores_good_cutoff = find_threshold(dataset, HQcolname, 'scores_Good', target_ratio=0.95)
    dataset['Good_module_res'] = dataset.apply(lambda row: label_uncertain(row, HQcolname, 'scores_Good', scores_good_cutoff), axis=1)
    Bad_train_data = np.array(dataset[(dataset['Good_module_res']=='Bad')][['ADz', 'HFz', 'HQz']])
    Bad_train_data_filtered = clean_data(Bad_train_data)
    if args.bad_module_n_components:
        scores_Bad, n_range, bics_Bad, best_gmm_Bad = GMM_fit(Bad_train_data_filtered, all_data_scaled, bad_module_n_components, max_n_components)
    else:
        scores_Bad, n_range, bics_Bad, best_gmm_Bad = GMM_fit(Bad_train_data_filtered, all_data_scaled)
    score_diff = scores_Good - scores_Bad
    dataset['scores_Bad'] = scores_Bad
    # Log-Likelihood Difference
    dataset['LLD'] = score_diff

    # calculate the ratio of Good and Bad values in each score_diff(Log-Likelihood Difference, LLD) 
    dataset = caculate_good_bad_ratio(dataset, 'Good_module_res', 'LLD', "Good_ratio", "Bad_ratio")

    # calculate the inflection points and optimal threshold
    best_inflection, inflection_points, inflection_res_info = calcu_inflection_points(dataset, "Bad_ratio", "Good_ratio", 'LLD', HQcolname, 'Good_module_res')
    optimal_threshold, roc_auc, fpr, tpr, df_sorted = calcu_optimal_threshold(dataset, best_inflection, 'LLD', 'Good_module_res')
    df_sorted['GMM_res'] = df_sorted.apply(lambda x: 'Good' if x['LLD'] > optimal_threshold else 'Bad', axis=1)

    # plot the results
    plt_gmm_results(n_range, bics_Good, best_gmm_Good, bics_Bad, best_gmm_Bad, df_sorted, 
                    inflection_points, fpr, tpr, roc_auc,optimal_threshold, args.figure_file)
    
    # output the log information
    output_info(df_sorted, HQcolname, good_train_data, Bad_train_data_filtered, scores_good_cutoff,
                best_inflection, inflection_res_info, optimal_threshold, args.output_csv_file.strip('.csv')+'.log')
    
    # output the final result to CSV file
    df_sorted.drop(columns=['Reported_AFgt005','ADz', 'HFz', 'HQz', 'scores_Good', 'Good_module_res', 'scores_Bad',
                        'Good_ratio', 'Bad_ratio', 'predicted', 'inflection_res'], inplace=True)
    dataset.to_csv(args.output_csv_file, index=False)
    
    # output the final result to VCF file
    POS_good_count_dict, VQSR_dict = df_info2dict(df_sorted)
    
    output_vcf_file = args.output_vcf_file.strip('.gz')
    with open(output_vcf_file, 'w') as output_vcf:
        for line in open_file(args.input_file):
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    output_vcf.write(f'##FORMAT=<ID=LLD,Number=R,Type=Float,Description="An ordered, comma delimited Log-Likelihood Difference for non-reference allele alleles in the order listed">\n')
                    output_vcf.write(f'##FORMAT=<ID=FILTER,Number=R,Type=String,Description="An ordered, comma delimited list of non-reference allele filtration situation ">\n')
                    output_vcf.write(f'##GMM_filtering_command=python GMM_classifier.py -i {args.input_file} -o {args.output_vcf_file} -c {args.output_csv_file} -f {args.figure_file}\n')
                    _CHROM,_POS,_ID,_REF,_ALT,_QUAL,_FILTER,_INFO,_FORMAT,*_SAMPLES = line.strip().split('\t')
                output_vcf.write(line)
            else:
                CHROM,POS,ID,REFs,ALTs,QUAL,FILTER,INFO,FORMAT, *SAMPLES = line.strip().split('\t')
                FORMAT = 'GT:GQ:DP:AD:HF:CI:HQ:LHF:SB:FS:SOR:VT:LLD:FILTER'
                AC = re.findall(r'AC=(.*?);', INFO)[0].split(',')[0]
                if POS_good_count_dict[POS] / int(AC) > 0.5:
                    FILTER = 'PASS'
                else:
                    FILTER = 'LowQual'
                newSAMPLES = []
                for sample_name, sample_data in zip(_SAMPLES, SAMPLES):
                    if sample_data.startswith('.'):
                        continue
                    GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT = sample_data.split(':')
                    LLD_list = []
                    LLD_res_list = []
                    for gt in GT.split('/'):
                        if gt == '0':
                            LLD_list.append("")
                            LLD_res_list.append('')
                        else:
                            try:
                                LLD_list.append(str(VQSR_dict[POS][sample_name][gt][0]))
                                LLD_res_list.append(VQSR_dict[POS][sample_name][gt][1])
                            except:
                                print(gt)
                    LLD = ','.join(LLD_list)
                    LLD_res = ','.join(LLD_res_list)
                    new_sample_data = f'{GT}:{GQ}:{DP}:{AD}:{HF}:{CI}:{HQ}:{LHF}:{SB}:{FS}:{SOR}:{VT}:{LLD}:{LLD_res}'
                    newSAMPLES.append(new_sample_data)
                print(CHROM,POS,ID,REFs,ALTs,QUAL,FILTER,INFO,FORMAT, *newSAMPLES, sep='\t', file=output_vcf)
    subprocess.run(["bgzip", "-f", output_vcf_file])
    subprocess.run(["tabix", "-p", "vcf", output_vcf_file+".gz"])

if __name__ == '__main__':
    main()