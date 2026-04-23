"""Perform mtDNA variant quality control analysis from VCF file.
This script parses mtDNA VCF file, estimates background noise from normal (non-mutated) samples, 
computes KL divergence for variant allele frequencies, and applies a Bayesian filter to determine 
the posterior probability of good calls.

Returns:
    It generates a filtered VCF file and a CSV report of the results.

Usage:
    python mtDNA_variant_QC.py --vcf input.vcf --output output.csv --output_vcf filtered_output.vcf
    [--bins 100] [--pi 5e-8] [--threshold 0.9] 
    
Author: Shujia Huang
Date: 2025-07-10
"""
import argparse
import sys
import os
import pysam

import numpy as np
import pandas as pd

from scipy.stats import beta, betabinom, binom
from scipy.optimize import minimize
from scipy.special import betaln

from typing import Dict, Generator, Optional
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend suitable for cluster/terminal
import matplotlib.pyplot as plt


GLOBAL_BLACKLIST_SITES_SET = set()
_GLOBAL_BLACKLISTED_REGIONS = [
    # 都是 chrM 所以简化为 (start, end) 即可
    (299, 317),
    (511, 525),
    (564, 571),
    (952, 955),
    (3105, 3108),
    (5895, 5899),
    (8268, 8279),
    (13645, 13650),
    (16180, 16187),
]
for start, end in _GLOBAL_BLACKLISTED_REGIONS:
    for pos in range(start, end + 1):
        GLOBAL_BLACKLIST_SITES_SET.add(pos)

def plot_beta_fit_convergence(alpha_hist, beta_hist, diff_hist, save_path='beta_fit_convergence.png'):
    plt.figure(figsize=(8, 5))
    ax1 = plt.gca()
    ax2 = ax1.twinx()

    # Plot alpha and beta on left y-axis
    ax1.plot(alpha_hist, label=r'${\alpha}$', color='tab:blue', marker='o')
    ax1.plot(beta_hist, label=r'${\beta}$', color='tab:orange', marker='s')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Beta parameter value')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Plot diff on right y-axis
    ax2.plot(diff_hist, label=r'$\delta$ (Mutation call change)', color='tab:red', marker='x', linestyle='--')
    ax2.set_ylabel('Mutation call change ratio', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    plt.title('Convergence of Beta Parameters and Mutation Call Change')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


def plot_beta_fit_and_vaf(alpha_hist, beta_hist, vaf_true_list, vaf_false_list, save_path='beta_fit_vaf.png'):
    plt.figure(figsize=(10, 6))
    x = np.linspace(0, 1, 200)
    ax1 = plt.gca()
    ax2 = ax1.twinx()

    # Plot Beta distributions for each iteration (faded lines) on ax1
    for i, (a, b) in enumerate(zip(alpha_hist, beta_hist)):
        ax1.plot(x, beta.pdf(x, a, b), color='tab:blue', alpha=0.2 + 0.6 * (i+1)/len(alpha_hist), lw=1, ls='--',
                 label='Beta dist (iter {})'.format(i+1) if i == len(alpha_hist)-1 else None)

    # Plot final Beta distribution (bold)
    ax1.plot(x, beta.pdf(x, alpha_hist[-1], beta_hist[-1]), color='tab:blue', lw=2, label=f'Final Beta({alpha_hist[-1]:.4f}, {beta_hist[-1]:.4f})')
    ax1.set_ylabel('Beta PDF', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Plot VAF histograms on ax2
    bins = 100
    ax2.hist(vaf_false_list, bins=bins, range=(0,1), color='tab:gray', alpha=0.4, label='VAF (is_mutation=False)')
    ax2.hist([v for v in vaf_true_list  if v < 0.99], bins=bins, range=(0,1), color='tab:orange', alpha=0.4, label='VAF (is_mutation=True)')
    
    # ax2.set_yscale('log')
    max_count = np.percentile(np.concatenate([np.histogram(vaf_true_list, bins=bins, range=(0,1))[0],
                                              np.histogram(vaf_false_list, bins=bins, range=(0,1))[0]]), 98)
    ax2.set_ylim(top=int(max_count))
    ax2.set_ylabel('VAF Count', color='tab:orange')

    # Plot VAF histograms as density
    # ax2.hist(vaf_false_list, bins=bins, range=(0,1), color='tab:gray', alpha=0.4, label='VAF (is_mutation=False)', density=True)
    # ax2.hist(vaf_true_list, bins=bins, range=(0,1), color='tab:orange', alpha=0.4, label='VAF (is_mutation=True)', density=True)
    # ax2.set_yscale('log')
    # ax2.set_ylabel('Density', color='tab:orange')
    ax2.tick_params(axis='y', labelcolor='tab:orange')

    # Legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    plt.xlabel('Variant Allele Frequency (VAF)')
    plt.title('Beta Distribution Fit and VAFs by Mutation Status')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

# Histogram of variant counts for all samples
def plot_variant_counts(sample_variant_count, save_path='variant_counts_per_sample.png'):
    all_counts = sample_variant_count['all']
    mut_counts = sample_variant_count['mutations']
    
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    ax.hist(all_counts, bins=50, color='tab:gray', alpha=0.4, label=f'Total Variants (Mean: {np.mean(all_counts):.0f})')
    ax.hist(mut_counts, bins=50, color='tab:orange', alpha=0.4, label=f'Good Variants (Mean: {np.mean(mut_counts):.0f})')
    ax.set_xlabel('Variant Counts per Sample')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Distribution of Variant Counts per Sample')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    

def check_index_file(variant_file_path):
    """
    Check if index file exists for the variant file
    
    Args:
        variant_file_path (str): Path to the VCF/BCF file
    
    Returns:
        tuple: (bool, str) - (index exists, index command if needed)
    """
    path = Path(variant_file_path)
    is_bcf = path.suffix.lower() == '.bcf'
    is_vcf_gz = path.suffix.lower() == '.gz' and path.stem.lower().endswith('.vcf')
    
    index_files = []
    if is_bcf:
        index_files = [f"{path}.csi", f"{path}.tbi"]
    elif is_vcf_gz:
        index_files = [f"{path}.tbi", f"{path}.csi"]
    else:
        raise ValueError("File must be BCF or compressed VCF (*.vcf.gz)")
    
    index_exists = any(os.path.exists(idx) for idx in index_files)
    
    # Prepare indexing command if needed
    if not index_exists:
        if is_bcf:
            cmd = f"bcftools index {variant_file_path}"
        else:
            cmd = f"tabix -p vcf {variant_file_path}"
        return False, cmd
    
    return True, ""


def variant_generator(
    file_path: str, 
    region: Optional[str]=None, 
    create_index: Optional[bool]=False
) -> Generator[Dict, None, None]:
    """
    Generator function to yield variants from a VCF/BCF file
    
    Args:
        file_path (str): Path to the VCF/BCF file
        region (str, optional): Genomic region to process (e.g., "chr1:1000-2000")
        create_index (bool): Whether to create index if missing
        
    Yields:
        dict: Dictionary containing variant information
    
    Raises:
        ValueError: If file format is invalid or index is missing
    """
    # Check file extension
    path = Path(file_path)
    if not (path.suffix.lower() == '.bcf' or 
            (path.suffix.lower() == '.gz' and path.stem.lower().endswith('.vcf'))):
        raise ValueError("File must be BCF or compressed VCF (*.vcf.gz)")
    
    # Check index
    index_exists, index_cmd = check_index_file(file_path)
    if not index_exists:
        if create_index:
            print(f"Creating index file using: {index_cmd}")
            os.system(index_cmd)
        else:
            raise ValueError(f"Index file missing. Create using command:\n{index_cmd}")
    
    vcf_file = None
    try:
        vcf_file = pysam.VariantFile(file_path, 'r')
        # Determine whether to fetch from specific region or entire file
        records = vcf_file.fetch(region=region) if region else vcf_file
        for record in records:
            variant_info = {
                'chrom': record.chrom,
                'pos': record.pos,
                'id': record.id,
                'ref': record.ref,
                'alts': record.alts,
                'qual': record.qual,
                'filter': list(record.filter.keys()),
                'info': dict(record.info),
                'format': record.format,
            }
            
            # Add samples information if present
            if record.samples:
                # Get all sample data: sample => record.samples[sample] is a dictionary of FORMAT fields for that sample
                variant_info['samples'] = {sample: record.samples[sample] for sample in record.samples}

            yield variant_info
            
    except Exception as e:
        sys.stderr.write(f"{str(e)} found when it's processing the file: {file_path}.\n")
        
    finally:
        if vcf_file is not None:
            vcf_file.close()
     

def write_vcf(input_vcf_path, output_vcf_path, results, pos_kl_div):
    """
    Write a new VCF file with PP and IS_MUT added to each sample's FORMAT.
    Args:
        input_vcf_path: str, original VCF file path
        output_vcf_path: str, output VCF file path
        results: list of dict, each dict contains sample, chrom, pos, pp, is_mutation
    """
    result_lookup = {}
    for r in results:
        key = (r['chrom'], r['pos'], r['sample'])
        # GT, PP, IS_MUT
        result_lookup[key] = (r['gt'], r['posterior'], r['is_mutation'])

    df_results = pd.DataFrame(results)
    qc_status = df_results.groupby(['chrom', 'pos'])['is_mutation'].max().reset_index().set_index(['chrom', 'pos']).to_dict()['is_mutation']
    with pysam.VariantFile(input_vcf_path) as vcf_in:
        # Add new FORMAT fields to header
        vcf_in.header.add_line('##FORMAT=<ID=PP,Number=1,Type=Float,Description="Posterior probability of true genotype">')
        vcf_in.header.add_line('##FORMAT=<ID=GOOD_CALL,Number=1,Type=String,Description="Good call (True=Good, False=bad)">')
        vcf_in.header.add_line('##FILTER=<ID=PASS,Description="Passed mtDNA QC filter">')
        vcf_in.header.add_line('##FILTER=<ID=BLACKLISTED_SITE,Description="Site is in the blacklisted regions">')
        vcf_in.header.add_line('##FILTER=<ID=LOW_QUALITY,Description="Low quality site based on KL divergence from background noise distribution">')
        vcf_in.header.add_line(f'##QC_command=python {" ".join(sys.argv)}')
        with pysam.VariantFile(output_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                for sample in record.samples:
                    key = (record.chrom, record.pos, sample)
                    if key in result_lookup:
                        gt, pp, is_mut = result_lookup[key]
                        record.samples[sample].update(
                            {'GT': gt, 
                             'PP': float(pp), 
                             'GOOD_CALL': 'True' if is_mut else 'False'}
                        )

                # Use record.filter.clear() and record.filter.add() to update FILTER.
                record.filter.clear()
                chrom_pos = (record.chrom, record.pos)
                if record.pos in GLOBAL_BLACKLIST_SITES_SET:
                    record.filter.add('BLACKLISTED_SITE')
                elif qc_status.get(chrom_pos, False):  
                    # If any sample at this position is called as mutation, 
                    # mark the site as PASS
                    record.filter.add('PASS')
                else:
                    record.filter.add('LOW_QUALITY')
                
                record.qual = round(pos_kl_div.get(chrom_pos, 0))
                vcf_out.write(record)
        
        # If output VCF is compressed, index it
        if output_vcf_path.endswith('.gz'):
            pysam.tabix_index(output_vcf_path, preset='vcf', force=True)


# A function to collect sample information of each variant by parsing FORMAT fields in VCF
def collect_samples_info(variant, samples):
    samples_info = []
    for sample in samples:
        dp   = variant['samples'][sample].get('DP')  # None if DP is missing for this sample or an int if DP is present
        ad   = variant['samples'][sample].get('AD')  # None if AD is missing for this sample or a tuple if AD is present
        hq   = variant['samples'][sample].get('AQ')  # None if AQ is missing for this sample or a tuple if AQ is present
        fss  = variant['samples'][sample].get('FS')  # None if FS is missing for this sample or a tuple if FS is present
        sors = variant['samples'][sample].get('SOR') # None if SOR is missing for this sample or a tuple if SOR is present
        gt   = list(variant['samples'][sample].get('GT'))
        hf   = list(variant['samples'][sample].get('AF'))
        
        if any(g is None for g in gt) if gt else True:
            sys.stderr.write(f"[WARNING] GT is missing for sample {sample} at {variant['chrom']}:{variant['pos']}. "
                             f"Setting GT to empty list.\n")
            gt = []
        if any(h is None for h in hf) if hf else True:
            sys.stderr.write(f"[WARNING] AF is missing for sample {sample} at {variant['chrom']}:{variant['pos']}. "
                             f"Setting AF to empty list.\n")
            hf = []

        # pysam doesn't support directly getting a tuple of tuples for the specific label of SB (format problem), 
        # so we need to parse it from the string format. Set back to be a string of format in 'a1_fwd,a1_reverse;a2_forward,a2_reverse;...' 
        # for each sample, where a1, a2, ... are ref and alt alleles in order.
        sb_s = ','.join(variant['samples'][sample].get('SB'))
        sb_split = sb_s.split(';') if sb_s != '.' else []
        
        # A list of tuples [(a1_forward, a1_reverse), (a2_forward, a2_reverse), ...]
        sbs = [tuple(map(int, sb.split(','))) for sb in sb_split] # could be empthy if SB is missing
        
        # Strand ratio factor can be calculated as np.min([fwd, rev])/np.max([fwd, rev]) for each allele.
        # The range of SRF is [0, 1], where values close to 0 indicate strong strand bias and values 
        # close to 1 indicate balanced strand representation.
        srf = [(np.min([fwd, rev]) + 1e-10)/(np.max([fwd, rev]) + 1e-10) for fwd, rev in sbs]
        samples_info.append({
            'sample': sample,
            'ploidy': len(gt) if gt else 0, # Determine ploidy based on GT length, if GT is missing, set ploidy to 0
            'GT': gt, # A list of allele indices for the sample, if GT is missing, set to empty list
            'AF': hf, # A list of allele frequencies for ref and alts, if AF is missing, set to empty list
            'HQ': hq, # A list of allele qualities for ref and alts, if AQ is missing, set to None
            'DP': dp,
            'AD': ad,
            'FS': fss,
            'SOR': sors,
            'SB': sbs,
            'SRF': srf  # SRF calculated from SB, can be used as a panelty factor in the Bayesian filter if needed
        })
        
    return samples_info

# A function to fetch background error rates from suitable samples for a given variant based on 
# their FORMAT fields and pre-filtering criteria
def fetch_background_error_rates(samples_info, args):
    background_error_rates = []
    for sample in samples_info:
        ploidy = sample['ploidy']
        hf = sample['AF']
        hq = sample['HQ']
        dp = sample['DP']
        fss = sample['FS']
        sors = sample['SOR']

        if ploidy != 1:
            continue  # Only use haploid samples for background error estimation, skip diploid or polyploid samples
        
        # Pre-filtering based on DP and HQ thresholds
        if dp is None or dp < args.DP_to_ploidy_threshold:
            continue
        
        if hq is None or any(q < args.HQ_threshold for q in hq):
            continue
        
        # 收集用于背景误差估计的合适样本时，我直接将 FS > 100 或 SOR > 5 的样本将被排除在外。
        # 因为对于 ploidy == 1 的样本，AF 肯定是高的，这时若 FS、SOR 也高，则它更可能是由于 
        # NUMTs 引起的系统性偏差，而非测序随机噪音（背景噪声）。
        if (fss is not None) and any(fs > 100 for fs in fss):
            continue
        if (sors is not None) and any(sor > 5 for sor in sors):
            continue

        # If ploidy is not 1 or AF is missing for this sample, skip it for background error estimation
        if (ploidy == 1) and (hf is not None) and all(h is not None for h in hf): 
            # Only use homozygous samples to estimate background error rate.
            background_error_rates.append(1.0 - sum(hf))
    
    return background_error_rates


def estimate_background_noise(background_error_rates, bins=100):
    """Estimate background noise distribution using normal samples (GT='0', '1', '2', ...)."""
    # Remove None and extreme values to avoid fitting errors
    vaf_arr = np.array([v for v in background_error_rates if v is not None and 0 < v < 1])
    if len(vaf_arr) < 2 or np.all(vaf_arr == 0):
        # If no valid background error rates, return a default Beta distribution that is heavily skewed 
        # towards 0 (indicating low error rate) and a uniform binning for KL divergence calculation.
        return 2.5, 1680, np.linspace(0, 1, bins + 1)  #两个默认值根据我的实际数据计算所得, np.linspace(0, 1, bins + 1)
    
    hist, bin_edges = np.histogram(vaf_arr, bins=bins, density=True, range=(0, 1))

    # Fit a Beta distribution to the background error rates using MLE. The `beta.fit` function 
    # returns the parameters of the fitted Beta distribution, which include the shape parameters 
    # alpha (a) and beta (b), as well as location and scale parameters. Since we are fitting a 
    # distribution for error rates that are between 0 and 1, we set the location to 0 (`floc=0`)
    # and scale to 1 (`fscale=1`) to ensure that the fitted distribution is valid for our data.
    a, b, loc, scale = beta.fit(vaf_arr, floc=0, fscale=1)

    return a, b, bin_edges  # return alpha and beta parameters for Beta distribution


# Fetch the high-quality VAFs for samples that are likely to have true mutations 
# (e.g., those with high VAFs and good quality metrics)
def ff(samples_info, p_error, args, ep17_alpha=0.01):
    vafs  = []  # The observed VAFs for samples that are likely to have true mutations
    k_obs = []  # The read counts supporting the ALT allele(s) for samples that are likely to have true mutations
    n_obs = []  # The total read depths for samples that are likely to have true mutations
    for sample in samples_info:
        gt = sample['GT']
        dp = sample['DP']
        hf = sample['AF']
        hq = sample['HQ']
        ad = sample['AD']
        
        if gt is None or hf is None or hq is None or dp is None:
            continue
        
        # ignore only ref samples
        if all(g == 0 for g in gt):
            continue

        if dp < args.DP_to_ploidy_threshold:
            continue
        
        if any(q < args.HQ_threshold for q in hq):
            continue
        
        # The limit of blanket read count for a given error rate p_error under the null hypothesis (no true mutation) 
        # can be calculated using the binomial distribution. We can use the percent point function (ppf) of the binomial 
        # distribution to find the read count threshold above which we would reject the null hypothesis at a significance
        # level of ep17
        lob_read_count = binom.ppf(1 - ep17_alpha, dp, p_error).astype(int)
        for g, a, vaf in zip(gt, ad, hf):
            if (a is not None) and (g > 0) and (a > lob_read_count):
                k_obs.append(a)
                n_obs.append(dp)
                vafs.append(vaf)
                
    return vafs, k_obs, n_obs


# Estimate the paramters of the true mutation distribution (Beta) using the background error rates 
# and method of MLE. The samples used for MLE should be those that are likely to have true mutations 
# (e.g., those with high VAFs and good quality metrics).
def estimate_true_mutation_params_betabinom(samples_info, p_error, args, ep17_alpha=0.01):
    """Estimate the parameters of the true mutation distribution (Betabinom) using the background error rates and method of MLE."""
    # Collect VAFs from samples that are likely to have true mutations 
    # (e.g., those with high VAFs and good quality metrics)
    vafs, k_obs, n_obs = ff(samples_info, p_error, args, ep17_alpha)
    return estimate_params_betabinom(vafs, k_obs, n_obs)       


def estimate_params_betabinom(vafs, k_obs, n_obs):
    """Estimate the parameters of the true mutation distribution (Betabinom) using the background error rates and method of MLE."""
    if len(vafs) < 1:
        sys.stderr.write(f"WARNING: Not enough variants ({len(vafs)}) for MLE. "
                         f"Using defaults (0.5, 0.2).\n")
        # Return default parameters by my own experiment if not enough data to fit
        return 0.5, 0.2 

    # Initial estimate of alpha and beta parameters using method of moments
    mu = np.mean(vafs)
    var = np.var(vafs)
    if var == 0:
        var = 1e-10  # Add a small value to avoid division by zero
    
    # 防止方差过大导致的非物理意义负数 (Beta分布要求方差 < mu*(1-mu))
    max_var = mu * (1.0 - mu)
    if var >= max_var:
        # Set var to be slightly less than the maximum possible variance 
        # for the given mean to ensure valid Beta parameters
        var = max_var * 0.99 
    
    # Calculate initial alpha and beta parameters using method of moments 
    # estimates as the starting point for MLE optimization.
    common_factor = (mu * (1.0 - mu) / var) - 1.0
    a_init = mu * common_factor
    b_init = (1.0 - mu) * common_factor
    
    # A function to calculate the negative log-likelihood of the observed data under 
    # the Beta-Binomial distribution with parameters a and b.
    def beta_neg_log_likelihood(params):
        a, b = params
        if a <= 0 or b <= 0:
            return 1e10  # Return a large number for invalid parameters
        
        # betabinom.pmf(k_obs, n_obs, a, b) gives the probability mass function of 
        # observing k_obs successes in n_obs trials under a Beta-Binomial distribution 
        # with parameters a and b. We take the log of this probability and sum it across 
        # all observations to get the total log-likelihood. Since we want to maximize the 
        # likelihood, we minimize the negative log-likelihood.
        log_probs = betabinom.logpmf(k_obs, n_obs, a, b)

        # If any log probability is -inf (which can happen if the parameters are not suitable 
        # for the data), we return a large number to penalize those parameters.
        log_probs = np.nan_to_num(log_probs, nan=-1e10, posinf=1e10, neginf=-1e10)
        return -np.sum(log_probs)
        
    # Use the initial estimates as the starting point for optimization
    try:
        bounds = [(0.01, 100.0), (0.01, 1000.0)]
        # Use L-BFGS-B algorithm for optimization, which is suitable for 
        # problems with bound constraints.
        result = minimize(
            beta_neg_log_likelihood, 
            x0=[a_init, b_init],  # Starting point for optimization
            method='L-BFGS-B',  
            bounds=bounds
        ) # 对于 Beta-Binomial 最大似然估计，L-BFGS-B 是目前最佳的优化算法之一，
          # 因为它能够处理参数的边界约束（alpha 和 beta 必须为正）。相比之下，
          # Nelder-Mead 不支持边界约束，可能会导致参数更新到无效的区域，从而产生非
          # 物理意义的结果。
        
        if result.success:
            a_mle, b_mle = result.x
            sys.stderr.write(f"- MLE optimization successful: alpha={a_mle:.6f}, beta={b_mle:.6f}\n")
        else:
            sys.stderr.write(f"- MLE optimization failed: {result.message}. Using method of "
                             f"moments estimates: alpha={a_init:.6f}, beta={b_init:.6f}\n")
            a_mle, b_mle = a_init, b_init
        
        return a_mle, b_mle
            
    except Exception as e:
        sys.stderr.write(f"[WARNING] MLE encountered error ({e}). Using default estimates: "
                         f"alpha={a_init:.6f}, beta={b_init:.6f}\n")
        return a_init, b_init


def qc(input_vcf_path, output_vcf_path, args):
    """Parse VCF file using pysam to extract variant information and store original records."""
    # Only get samples id from input vcffile
    _vcf = pysam.VariantFile(input_vcf_path, 'r')
    samples = list(_vcf.header.samples)
    _vcf.close()
    
    vafs, k_obs, n_obs = [], [], []  # Collect VAFs, supporting read counts, and total read depths across all variants for multi-sample KL divergence calculation
    results = []
    pos_kl_div = {} # Store multi-sample KL divergence for each position
    for variant in variant_generator(input_vcf_path):
        # print(f"{variant['chrom']}:{variant['pos']} {variant['ref']}->{variant['alts']} - QUAL: {variant['qual']} FILTER: {variant['filter']} INFO: {variant['info']}")
        ref_len = len(variant['ref'])
        if any(pos in GLOBAL_BLACKLIST_SITES_SET for pos in range(variant['pos'], variant['pos'] + ref_len)):
            sys.stderr.write(f"[INFO] Variant at {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alts']} "
                             f"is located in blacklisted regions. Skipping.\n")
            continue
        
        # Maximum ALT alleles per site
        if len(variant['alts']) > args.max_alt_alleles:
            sys.stderr.write(f"[INFO] Variant at {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alts']} "
                             f"has more than {args.max_alt_alleles} ALT alleles. Skipping.\n")
            continue

        ref_alts = [variant.get('ref', '')] + list(variant.get('alts', []))

        # Collect sample information for this variant by parsing FORMAT fields in VCF
        # A list of dictionaries, each dictionary contains all relevant information of 
        # a sample for this variant, such as GT, AF, DP, AD, HQ, FS, SOR, SB, etc.
        samples_info = collect_samples_info(variant, samples) 

        # Estimate background error rates from suitable samples for this variant
        background_error_rates = fetch_background_error_rates(samples_info, args)
        if len(background_error_rates) == 0:
            sys.stderr.write(f"[WARNING] All heterozygous samples. No background "
                             f"error rate calculated for this position: "
                             f"{variant['chrom']}:{variant['pos']}. Ignore.\n")
            continue

        # Calculate mean background error rate
        p_error = np.mean(background_error_rates)
        
        # Estimate background noise distribution using normal (homozygous) samples error rates
        # The prefix "q_" is used to distinguish the parameters of the background noise distribution (Beta) 
        # from the parameters of the true mutation distribution (Beta) in the Bayesian filter, which will be 
        # updated iteratively later.
        q_alpha, q_beta, bin_edges = estimate_background_noise(background_error_rates)
        sys.stderr.write(
            f"Processing variant at {variant['chrom']}:{variant['pos']} background error "
            f"rate: {p_error} - Parameters of Beta Distribution for background noise of "
            f"this position: alpha={q_alpha:.6f} beta={q_beta:.6f}\n"
        )
        
        # The prefix "h1_" is used to denote the parameters of the true mutation distribution (Beta) in the Bayesian filter,
        # which will be updated iteratively later based on the observed data and the current mutation calls. 
        # The initial values can be set to Beta(1,1) or estimated from the data using the method of moments or MLE. 
        # Here we use MLE to get a more informed initial estimate for the true mutation distribution. 
        # h1_alpha, h1_beta = estimate_true_mutation_params_betabinom(samples_info, p_error, args)
        # sys.stderr.write(f"Estimated Beta distribution parameters for true mutations at "
        #                  f"{variant['chrom']}:{variant['pos']} : alpha={h1_alpha} beta={h1_beta}\n")
        
        vaf, k, n = ff(samples_info, p_error, args)
        vafs.extend(vaf)
        k_obs.extend(k)
        n_obs.extend(n)
        
        # QC for each sample
        all_available_vaf_obs = []
        for sample in samples_info:
            sample_name = sample['sample']
            gt  = sample['GT']
            dp  = sample['DP']
            vaf = sample['AF']
            hqs = sample['HQ']
            ads = sample['AD']
            srf = sample['SRF']
            
            # Pre-filtering based on DP and AQ thresholds
            is_pre_filtered = False
            # if (dp is None) or dp < ploidy * args.DP_to_ploidy_threshold:
            #     sys.stderr.write(f"[WARNING] Sample {sample_name} at {variant['chrom']}:{variant['pos']} failed "
            #                      f"DP/Ploidy threshold ({dp}/{ploidy} < {args.DP_to_ploidy_threshold}). Skipping.\n")
            #     is_pre_filtered = True
            if ((dp is None) 
                or any(f is None for f in vaf) 
                or any(a is None for a in ads) 
                or any(h is None for h in hqs)):
                sys.stderr.write(
                    f"[WARNING] Missing DP, AF, AD or HQ for sample {sample_name} "
                    f"at {variant['chrom']}:{variant['pos']}. Skipping.\n"
                )
                is_pre_filtered = True
            else :
                # Calculate KL divergence for each ALT allele observed in this sample and 
                # use the average KL divergence across all ALTs as the final KL divergence 
                # for this sample at this position, which will be used in the Bayesian filter 
                # to compute the posterior probability of true mutation.
                for v_obs, a_obs in zip(vaf, ads):
                    if (v_obs is None) or (a_obs is None):
                        sys.stderr.write(f"[WARNING] Missing VAF or AD for sample {sample_name} "
                                         f"at {variant['chrom']}:{variant['pos']}. Skipping.\n")
                        break
                    
                    all_available_vaf_obs.append(v_obs)

            results.append({
                'sample': sample_name,
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': [ref_alts[g_idx] if (g_idx is not None) else '.' for g_idx in gt] if gt else ['.'],
                
                # Keep original GT for reference, but the final mutation 
                # call will be determined by 'is_mutation' field.
                'gt': gt,
                'ploidy': len(gt),
                
                'vaf': vaf,
                'HQ': hqs,
                'A': ads,
                'D': dp,
                'srf': srf,
                
                'q_alpha': q_alpha,  # Store the background noise distribution parameters for this position, which will be used in the Bayesian filter to compute the posterior probability of true mutation for each sample at this position. These parameters will be updated iteratively later based on the observed data and the current mutation calls.
                'q_beta': q_beta,    # Store the background noise distribution parameters for this position, which will be used in the Bayesian filter to compute the posterior probability of true mutation for each sample at this position. These parameters will be updated iteratively later based on the observed data and the current mutation calls.
                'p_error': p_error,

                'posterior': 0.0,      # update pp in each iteration of the iterative_beta_fit_and_call function, which will be used to update mutation calls and re-estimate the true mutation distribution parameters until convergence.
                'is_mutation': False,  # Initial mutation call based on threshold, updated iteratively later in the iterative_beta_fit_and_call function. This is to avoid making hard mutation calls before we have a stable estimate of the true mutation distribution parameters.
                
                'pre_filtered': is_pre_filtered
            })
        
        # Calculate multi-sample KL divergence for the position
        kl_div_multi = calculate_kl_divergence_multi(all_available_vaf_obs, q_alpha, q_beta, bin_edges)
        pos_kl_div[(variant['chrom'], variant['pos'])] = kl_div_multi if np.isfinite(kl_div_multi) else 1000

    # Iteratively estimate beta distribution parameters from called mutations and update mutation calls until convergence. 
    # For multi-allelic sites, treat each ALT as an independent event and use the product of posteriors as the final posterior for the site.
    results, alpha_hist, beta_hist, diff_hist = iterative_beta_fit_and_call(
        results, np.array(vafs), np.array(k_obs), np.array(n_obs),
        hq_threshold=args.HQ_threshold, pi=args.pi, 
        threshold=args.threshold
    )
    sys.stderr.write(f"\nTotal variants processed: {len(results)}, Global Beta distribution for "
                     f"true variants: alpha={alpha_hist[-1]:.6f} beta={beta_hist[-1]:.6f}\n")
    
    sample_variant_count = {
        'all': [0 for _ in samples], 
        'mutations': [0 for _ in samples]
    }
    sample_index = {s:i for i, s in enumerate(samples)}
    
    vaf_true_list = []
    vaf_false_list = []
    for r in results:
        if any(alt != '.' and alt != r['ref'] for alt in r['alt']):  # raw variant exists
            sample_variant_count['all'][sample_index[r['sample']]] += 1
            
        if r['is_mutation']:
            for g, v in zip(r['gt'], r['vaf']):
                if g is not None and g > 0 and v is not None and 0 < v < 1:
                    vaf_true_list.append(v)
                
            # if any ALT allele is observed in the original GT, count this sample as having 
            # a mutation at this position
            if any(gt for gt in r['gt'] if gt is not None and gt > 0):  
                sample_variant_count['mutations'][sample_index[r['sample']]] += 1
        else:
            vaf_false_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])

        # r['kl_divergence_single'] = np.mean(r['kl_divergence_single']) if r['kl_divergence_single'] else 0
        r.pop('vaf', None)  # Remove VAF from final results
        r.pop('A', None)    # Remove A from final results
        r.pop('D', None)    # Remove D from final results
          
    write_vcf(input_vcf_path, output_vcf_path, results, pos_kl_div)
    return results, sample_variant_count, vaf_true_list, vaf_false_list, alpha_hist, beta_hist, diff_hist


def __old_iterative_beta_fit_and_call(
    results, vafs, k_obs, n_obs,
    hq_threshold=20,
    pi=5e-8 * 16569,
    threshold=0.9,
    max_iter=50,
    tol=1e-5
):
    """
    Iteratively estimate beta distribution parameters from called mutations and update mutation calls
    until convergence. For multi-allelic sites, treat each ALT as an independent event and use the
    product of posteriors as the final posterior for the site.
    """
    alpha_hist, beta_hist, diff_hist = [1.0], [1.0], []
    prev_is_mut = [r['is_mutation'] for r in results]
    
    # 1. Collect VAFs from currently called mutations
    high_VAF_indices = []
    low_VAF_indices = []
    for i, vaf in enumerate(vafs):
        if vaf is not None and 0 < vaf < 1:
            if vaf > 0.9:
                high_VAF_indices.append(i)
            else:
                low_VAF_indices.append(i)
    
    sampling_num = min(len(high_VAF_indices), len(low_VAF_indices))      
    for i in range(max_iter):
        if sampling_num > 0:
            h_indices_sampled = np.random.choice(high_VAF_indices, size=sampling_num, replace=False)
            l_indices_sampled = np.random.choice(low_VAF_indices, size=sampling_num, replace=False)
            
            vafs_sampled  = list(vafs[h_indices_sampled]) + list(vafs[l_indices_sampled])
            k_obs_sampled = list(k_obs[h_indices_sampled]) + list(k_obs[l_indices_sampled])
            n_obs_sampled = list(n_obs[h_indices_sampled]) + list(n_obs[l_indices_sampled])
        else:
            vafs_sampled  = list(vafs)
            k_obs_sampled = list(k_obs)
            n_obs_sampled = list(n_obs)

        # 2. Fit Beta distribution parameters using MLE
        alpha_h1, beta_h1 = estimate_params_betabinom(vafs_sampled, k_obs_sampled, n_obs_sampled)
        sys.stderr.write(f"- Iteration {i+1}: {len(vafs_sampled)} VAFs (first-5: {vafs_sampled[:5]}) "
                         f"collected for Beta fitting. Estimated Beta distribution parameters for "
                         f"true mutations: alpha={alpha_h1:.6f}, beta={beta_h1:.6f}\n")

        # Ensure alpha_h1 and beta_h1 are valid
        alpha_hist.append(alpha_h1)
        beta_hist.append(beta_h1)
        
        # 3. Recalculate posterior and is_mutation for each record
        for r in results:
            if r['pre_filtered']:
                continue  # Skip pre-filtered samples
            
            # For multi-allelic sites, calculate posterior for each ALT and take the product
            new_pps = []
            r['new_gt'] = r['gt']  # Initialize new GT to original GT, will update based on new mutation calls if needed
            # for a_obs, kl_div, srf, hq in zip(r['A'], r['kl_divergence_single'], r['srf'], r['HQ']):
            for a_obs, srf, hq in zip(r['A'], r['srf'], r['HQ']):
                new_pp = bayesian_filter(
                    A=a_obs,
                    D=r['D'],
                    srf=srf,
                    hq=hq,
                    hq_threshold=hq_threshold,
                    
                    q_alpha=r['q_alpha'],
                    q_beta=r['q_beta'],
                    
                    alpha_h1=alpha_h1,  # Use updated Beta parameters
                    beta_h1=beta_h1,    # Use updated Beta parameters
                    
                    # Use the same prior probability for true mutationas before, 
                    # which can be adjusted based on the expected mutation rate and 
                    # the size of the genomic region being analyzed.
                    pi=pi
                )
                new_pps.append(new_pp)

            r['posterior'] = np.max(new_pps) if new_pps else 0  # For multi-allelic sites, take the maximum posterior across all ALTs as the final posterior for the site. This is a more conservative approach than taking the product, as it assumes that the evidence from each ALT is not independent and that the strongest evidence should dominate the mutation call.
            r['is_mutation'] = r['posterior'] > threshold       # Update mutation call based on new posterior
            r['new_gt'] = [g if (pp > threshold and g is not None) else None for pp, g in zip(new_pps, r['gt'])]  # Update GT to only include alleles that are called as mutations based on the new posterior (optional, can keep original GT for reference)
            
        # 4. Check for convergence
        curr_is_mut = [r['is_mutation'] for r in results]
        if prev_is_mut is not None:
            diff = np.mean([a != b for a, b in zip(prev_is_mut, curr_is_mut)])
            diff_hist.append(diff)
            if diff < tol:
                break
            
        prev_is_mut = curr_is_mut.copy()

    # After convergence, update the final GT calls based on 
    for r in results:
        if r['pre_filtered']:
            continue  # Skip pre-filtered samples
            
        r['gt'] = r['new_gt']  # Update GT to the new GT based on final mutation calls
        r.pop('new_gt', None)  # Remove the temporary 'new_gt' field from the final results to clean up the output.

    return results, alpha_hist, beta_hist, diff_hist


def iterative_beta_fit_and_call(
    results, vafs, k_obs, n_obs,
    hq_threshold=20,
    pi=5e-8 * 16569,
    threshold=0.9
):
    """
    Estimate beta distribution parameters from called mutations and update mutation calls
    until convergence. For multi-allelic sites, treat each ALT as an independent event and 
    use the maximum of posteriors as the final posterior for the site.
    """
    alpha_hist, beta_hist, diff_hist = [1.0], [1.0], []
    
    # 1. Collect VAFs from currently called mutations
    high_VAF_indices = []
    low_VAF_indices = []
    for i, vaf in enumerate(vafs):
        if vaf is not None and 0 < vaf < 1:
            if vaf > 0.95:
                high_VAF_indices.append(i)
            else:
                low_VAF_indices.append(i)
    
    sampling_num = min(len(high_VAF_indices), len(low_VAF_indices))
    by_step = max(len(high_VAF_indices), len(low_VAF_indices)) // sampling_num if (
        sampling_num > 0
    ) else 1
    if len(high_VAF_indices) > len(low_VAF_indices):
        high_VAF_indices = high_VAF_indices[::by_step]
    else:
        low_VAF_indices = low_VAF_indices[::by_step]
    
    if sampling_num > 0:
        vafs_sampled  = list(vafs[low_VAF_indices]) + list(vafs[high_VAF_indices])
        k_obs_sampled = list(k_obs[low_VAF_indices]) + list(k_obs[high_VAF_indices])
        n_obs_sampled = list(n_obs[low_VAF_indices]) + list(n_obs[high_VAF_indices])
    else:
        vafs_sampled  = list(vafs)
        k_obs_sampled = list(k_obs)
        n_obs_sampled = list(n_obs)
    
    # 2. Fit Beta distribution parameters using MLE
    alpha_h1, beta_h1 = estimate_params_betabinom(vafs_sampled, k_obs_sampled, n_obs_sampled)
    sys.stderr.write(f"Training dataset: {len(vafs_sampled)} VAFs (first-5: {vafs_sampled[:5]}) "
                     f"collected for Beta fitting. Estimated Beta distribution parameters for "
                     f"true mutations: alpha={alpha_h1:.6f}, beta={beta_h1:.6f}\n")

    # Ensure alpha_h1 and beta_h1 are valid
    alpha_hist.append(alpha_h1)
    beta_hist.append(beta_h1)
    
    # 3. Recalculate posterior and is_mutation for each record
    for r in results:
        if r['pre_filtered']:
            continue  # Skip pre-filtered samples
        
        # For multi-allelic sites, calculate posterior for each ALT and take the product
        new_pps = []
        r['new_gt'] = r['gt']  # Initialize new GT to original GT, will update based on new mutation calls if needed
        # for a_obs, kl_div, srf, hq in zip(r['A'], r['kl_divergence_single'], r['srf'], r['HQ']):
        for a_obs, srf, hq in zip(r['A'], r['srf'], r['HQ']):
            new_pp = bayesian_filter(
                A=a_obs,
                D=r['D'],
                srf=srf,
                hq=hq,
                hq_threshold=hq_threshold,
                
                q_alpha=r['q_alpha'],
                q_beta=r['q_beta'],
                
                alpha_h1=alpha_h1,  # Use updated Beta parameters
                beta_h1=beta_h1,    # Use updated Beta parameters
                
                # Use the same prior probability for true mutationas before, 
                # which can be adjusted based on the expected mutation rate and 
                # the size of the genomic region being analyzed.
                pi=pi
            )
            new_pps.append(new_pp)

        # r['posterior'] = np.prod(new_pps)  # For multi-allelic sites, take the product of posteriors for all ALTs as the final posterior for the site. This is a simple way to combine evidence from multiple ALTs, but it assumes that the evidence from each ALT is independent, which may not always be the case. More sophisticated methods could be used to combine evidence from multiple ALTs if needed.
        # r['posterior'] = np.mean(new_pps) if new_pps else 0  # For multi-allelic sites, take the average posterior across all ALTs as the final posterior for the site. This is a more balanced approach than taking the product, as it assumes that the evidence from each ALT contributes equally to the mutation call without assuming independence or dominance.
        r['posterior'] = np.max(new_pps) if new_pps else 0  # For multi-allelic sites, take the maximum posterior across all ALTs as the final posterior for the site. This is a more conservative approach than taking the product, as it assumes that the evidence from each ALT is not independent and that the strongest evidence should dominate the mutation call.
        r['is_mutation'] = r['posterior'] > threshold       # Update mutation call based on new posterior
        r['new_gt'] = [g if (g is not None and pp > threshold) else None for pp, g in zip(new_pps, r['gt'])]  # Update GT to only include alleles that are called as mutations based on the new posterior (optional, can keep original GT for reference)

    # After convergence, update the final GT calls based on 
    for r in results:
        if r['pre_filtered']:
            continue  # Skip pre-filtered samples
            
        r['gt'] = r['new_gt']  # Update GT to the new GT based on final mutation calls
        r.pop('new_gt', None)  # Remove the temporary 'new_gt' field from the final results to clean up the output.

    return results, alpha_hist, beta_hist, diff_hist


def calculate_kl_divergence_single(v_obs, a, b):
    """Calculate KL divergence for a single sample VAF."""
    if v_obs is None or v_obs <= 0 or v_obs >= 1:
        return 0  # Invalid VAF value

    q_v = beta.pdf(v_obs, a, b)
    sys.stderr.write(f" - Q({v_obs:.4f}|{a:.4f}, {b:.4f}) = {q_v}\n")
    return -np.log(q_v) if q_v > 0 else np.inf


def calculate_kl_divergence_multi(vafs, q_alpha, q_beta, bin_edges):
    """Calculate KL divergence for multi-sample VAF distribution."""
    if len(vafs) == 0:
        return 0
    
    px_at_bins, _ = np.histogram(vafs, bins=bin_edges, density=True)  # p_x
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    qx_at_bins = beta.pdf(bin_centers, q_alpha, q_beta)
    qx_at_bins = np.where(qx_at_bins > 0, qx_at_bins, 1e-10)  # Avoid division by zero

    delta_bin_width = bin_edges[1] - bin_edges[0]
    kl_div = np.sum(px_at_bins * np.log(px_at_bins / qx_at_bins + 1e-16) * delta_bin_width)

    return kl_div if kl_div > 0 else 0


# Bayesian filter to compute posterior probability of a true mutation using log-likelihoods, 
# with improved handling of edge cases and numerical stability.
# 2026-04-15: Updated Bayesian filter with improved handling of edge cases and numerical stability
def bayesian_filter(
    A,
    D,
    srf,
    hq,
    hq_threshold,
    q_alpha=1.0,  
    q_beta=1.0,  
    alpha_h1=1.0,
    beta_h1=1.0,
    pi=5e-8 * 16569
):
    """Compute posterior probability for H1 (true mutation) vs H0 (background error).

    H0: observed alternate counts come from a Beta-Binomial distribution with
        Beta(q_alpha, q_beta).
    H1: observed alternate counts come from a Beta-Binomial distribution with
        Beta(alpha_h1, beta_h1).

# 这段注释可以删掉了，因为我现在不使用 KL divergence 来调整后验概率了。
# `kl_div` is the KL divergence of the observed VAF distribution from the background noise distribution, 
# which serves as an additional piece of evidence to distinguish true mutations from background errors. 
# A higher KL divergence indicates that the observed VAF distribution is more different from the expected 
# background noise distribution, which can increase the likelihood of H1 (true mutation) being the correct 
# model. The `kl_weight` parameter allows us to adjust how much influence this KL divergence evidence has on 
# the final posterior probability calculation. By incorporating KL divergence into the Bayesian filter, we 
# can leverage the information about how well the observed data fits the expected background noise model to 
# make more informed decisions about whether a variant is likely to be a true mutation or just a result of 
# sequencing errors.
    
    Args:
        A: int, observed alternate allele count
        D: int, total depth
        srf: float, strand ratio factor
        hq: float, allele quality
        hq_threshold: float, HQ threshold
        q_alpha: float, alpha parameter for background Beta distribution under H0
        q_beta: float, beta parameter for background Beta distribution under H0
        alpha_h1: float, alpha parameter for true-mutation Beta distribution under H1
        beta_h1: float, beta parameter for true-mutation Beta distribution under H1
        pi: float, prior probability of mutation (H1)
    
    """
    if A is None or D is None:
        return 0.0

    try:
        A = int(A)
        D = int(D)
    except (TypeError, ValueError):
        return 0.0

    if D <= 0 or A < 0 or A > D:
        return 0.0

    hq_k=10

    # Sigmoid penalty functions for SRF and HQ to smoothly penalize low-quality evidence for mutations, 
    # which helps to reduce false positives by down-weighting the likelihood of H1 when quality metrics 
    # are poor.
    penalty_srf = 1 - srf
    
    # HQ高于阈值，免除惩罚；低于阈值，增加惩罚.
    penalty_hq = 1.0 / (1.0 + np.exp(hq_k * (hq_threshold - hq))) if hq < hq_threshold else 1.0
    log_penalty = np.log(penalty_hq * penalty_srf + 1e-12)  # Add small value to avoid log(0)

    # Log-likelihood under H0
    # The log_L0 is calculated using the binomial distribution PMF, which directly models the probability 
    # of observing A successes (alternate allele reads) in D trials (total reads) given a success probability 
    # of p_error (background error rate). This is a more direct and appropriate way to calculate the likelihood 
    # under H0 for count data, as it accounts for the discrete nature of the data and the specific error model 
    # we are using. The other commented-out functions are more general forms that could be used in different 
    # contexts, but for our specific case of modeling read counts with a known error rate, the binomial PMF is 
    # the most suitable choice. The math formula for log_L0 using the binomial PMF is:
    # log_L0 = log(P(X=A | D, p_error)) 
    #        = log(binom.pmf(A, D, p_error)) 
    #        = log(D choose A) + A*log(p_error) + (D - A)*log(1 - p_error)
    # 这是按照之前的版本计算的log_L0，在H0下直接建模了一个binomial分布来描述背景错误率。
    # 后面我要改成 beta-binomial 来更好地建模背景噪音分布，所以这个版本的log_L0就不再适用了。
    # p_error = np.clip(p_error, 1e-12, 1.0 - 1e-12)
    # log_L0 = binom.logpmf(A, D, p_error)

    # Log-likelihood under H0 (now using Beta-Binomial with background Beta parameters)
    # The log_L0 is calculated using the Beta-Binomial distribution PMF, which models the probability of 
    # observing A successes in D trials given a success probability that follows a Beta distribution with 
    # parameters q_alpha and q_beta. This is appropriate for H0 because we are modeling the background noise 
    # distribution as a Beta-Binomial, which accounts for overdispersion in the data that can arise from 
    # biological variability and technical noise. The math formula for log_L0 using the Beta-Binomial PMF is:
    # log_L0 = log(P(X=A | D, q_alpha, q_beta)) 
    #        = log(betabinom.pmf(A, D, q_alpha, q_beta)) 
    #        = log(Gamma(D+1)) - log(Gamma(A+1)) - log(Gamma(D-A+1)) + 
    #          log(Beta(A + q_alpha, D - A + q_beta)) - log(Beta(q_alpha, q_beta))
    
    # log_L0 = betabinom.logpmf(A, D, q_alpha, q_beta)  # 这是完整调用，但我不需要完整调用，因为我只需要后面两项来计算对数似然比，前面那三项在计算对数似然比时会相互抵消掉。
    # without the constant term log(Gamma(D+1)) - log(Gamma(A+1)) - log(Gamma(D-A+1)) which will cancel out in the posterior odds ratio anyway.
    log_L0 = betaln(A + q_alpha, D - A + q_beta) - betaln(q_alpha, q_beta)  
    
    # Log-likelihood under H1
    # The log_L1 is calculated using the Beta-Binomial distribution PMF, which models the probability of 
    # observing A successes in D trials given a success probability that follows a Beta distribution with 
    # parameters alpha_h1 and beta_h1. This is appropriate for H1 because we are modeling the true mutation 
    # distribution as a Beta-Binomial, which accounts for overdispersion in the data that can arise from 
    # biological variability and technical noise. The math formula for log_L1 using the Beta-Binomial PMF is:
    # log_L1 = log(P(X=A | D, alpha_h1, beta_h1)) 
    #        = log(betabinom.pmf(A, D, alpha_h1, beta_h1)) 
    #        = log(Gamma(D+1)) - log(Gamma(A+1)) - log(Gamma(D-A+1)) + 
    #          log(Beta(A + alpha_h1, D - A + beta_h1)) - log(Beta(alpha_h1, beta_h1))
    
    # log_L1 = betabinom.logpmf(A, D, alpha_h1, beta_h1)  # 这是完整调用，但我不需要完整调用，因为我只需要后面两项来计算对数似然比，前面那三项在计算对数似然比时会相互抵消掉。
    # without the constant term log(Gamma(D+1)) - log(Gamma(A+1)) - log(Gamma(D-A+1)) which will cancel out in the posterior odds ratio anyway.
    log_L1 = betaln(A + alpha_h1, D - A + beta_h1) - betaln(alpha_h1, beta_h1)  
    
    # Incorporate the fused penalty into the H1 to further penalize low-quality evidence 
    # for mutations, which helps to reduce false positives by down-weighting the posterior 
    # odds in favor of H1 when quality metrics are poor.
    log_L1 += log_penalty
    
    # Log prior odds of H1 vs H0, calculated from the prior probability pi of mutation (H1). 
    # The log prior odds is the logarithm of the ratio of the prior probabilities of H1 and H0, 
    # which is log(pi / (1 - pi)). This term represents our prior belief about how likely a true 
    # mutation is compared to a background error before observing the data. A higher pi will 
    # increase the log prior odds in favor of H1, while a lower pi will decrease it. 
    pi = np.clip(pi, 1e-12, 1.0 - 1e-12)
    log_prior_odds = np.log(pi) - np.log1p(-pi)

    # Log-posterior odds
    log_posterior_odds = log_prior_odds + log_L1 - log_L0

    # Posterior probability: 该样本/该等位基因在当前深度和背景错误率下，属于真实突变（H1）的概率
    if np.isfinite(log_posterior_odds):
        posterior_h1 = 1.0 / (1.0 + np.exp(-log_posterior_odds))
    else:
        posterior_h1 = 1.0 if log_posterior_odds > 0 else 0.0

    sys.stderr.write(
        f">> Bayesian filter: A={A}, D={D}, srf={srf:.6f}, hq={hq}, "
        f"q_alpha={q_alpha:.6f}, q_beta={q_beta:.6f}, "
        f"alpha_h1={alpha_h1:.6f}, beta_h1={beta_h1:.6f}, pi={pi:.6f}, "
        f"log_prior_odds={log_prior_odds:.6f}, log_penalty={log_penalty:.6f}, "
        f"log_posterior_odds={log_posterior_odds:.6f}, "
        f"posterior_h1={posterior_h1:.6f}\n"
    )
    
    return float(np.clip(posterior_h1, 0.0, 1.0))


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="Perform mtDNA variant quality control analysis from VCF file.")
    parser.add_argument('--vcf', required=True, help="Input VCF file path")
    parser.add_argument('--output', required=True, help="Output CSV file path")
    parser.add_argument('--output-vcf', required=True, help="Output quality-controlled VCF file path")
    
    # For maximum ALT alleles per site
    parser.add_argument('--max-alt-alleles', type=int, default=2, help="Maximum number of ALT alleles per site to consider. Default is 2")
    
    # For per-sample pre-filtering parameters
    parser.add_argument('--DP-to-ploidy-threshold', type=int, default=100, 
                        help="Parameter for training. Minimum depth to ploidy (DP/ploidy) threshold for considering variants "
                             "for each sample. Default is 100")
    parser.add_argument('--HQ-threshold', type=int, default=20, 
                        help="Parameter for training. Minimum variant quality (HQ) threshold for considering "
                             "variants for each sample. Default is 20")
    
    # parser.add_argument('--FS-threshold', type=float, default=100.0, 
    #                     help="Maximum Fisher Strand (FS) value threshold for considering "
    #                          "variants for each sample. Default is 100.0")
    # parser.add_argument('--SOR-threshold', type=float, default=5.0, 
    #                     help="Maximum Strand Odds Ratio (SOR) value threshold for considering "
    #                          "variants for each sample. Default is 5.0")
    
    # For Bayesian filter parameters
    parser.add_argument('--bins', type=int, default=100, help="Number of histogram bins. Default is 100")
    parser.add_argument('--pi', type=float, default=5e-8 * 16569, help="Prior probability of mutation. Default is 5e-8 * 16569")
    parser.add_argument('--threshold', type=float, default=0.9, help="Posterior probability threshold for calling mutations. Default is 0.9")
    args = parser.parse_args()
    
    try:
        # Parse VCF file
        (variants, 
         sample_variant_count, 
         vaf_true_list, 
         vaf_false_list, 
         alpha_hist, 
         beta_hist, 
         diff_hist) = qc(args.vcf, args.output_vcf, args)
        print(f'diff_hist: {diff_hist}\nalpha_hist: {alpha_hist}\nbeta_hist: {beta_hist}')
        
        # Save results to CSV
        df_results = pd.DataFrame(variants).explode('alt')
        # Filter out records where ref==alt and keep only mutations, then drop the 'gt' column if present
        df_results = df_results[
            (df_results['ref'] != df_results['alt']) & 
            df_results['is_mutation']
        ].drop(columns=['gt', 'pre_filtered'], errors='ignore')
        df_results.to_csv(args.output, sep='\t', index=False)
        
        # Plot convergence of beta parameters and mutation calls
        plot_beta_fit_convergence(alpha_hist, beta_hist, diff_hist, 
                                  save_path=args.output.split('.')[0] + 
                                  '_beta_fit_convergence.png')
        
        # plot Beta distribution fit and VAFs
        plot_beta_fit_and_vaf(alpha_hist, beta_hist, vaf_true_list, vaf_false_list, 
                              save_path=args.output.split('.')[0] + '_beta_fit_vaf.png')
        print(f"Quality control analysis completed. Results saved to {args.output_vcf}")
        
        # Plot variant counts per sample
        plot_variant_counts(sample_variant_count, save_path=args.output.split('.')[0] + '_variant_counts_per_sample.png')
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    