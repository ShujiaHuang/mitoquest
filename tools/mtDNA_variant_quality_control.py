
"""Perform mtDNA variant quality control analysis from VCF file.
This script parses a VCF file, estimates background noise from normal samples,
computes KL divergence for variant allele frequencies, and applies a Bayesian filter
to determine the posterior probability of true mutations.

Returns:
    It generates a filtered VCF file and a CSV report of the results.

Usage:
    python mtDNA_variant_quality_control.py --vcf input.vcf --output output.csv --output_vcf filtered_output.vcf
    [--bins 100] [--lambda_kl 0.1] [--pi 5e-8] [--threshold 0.9] 
    
Author: Shujia Huang
Date: 2025-07-10
"""

import argparse
import sys
import os
import pysam

import numpy as np
import pandas as pd

from scipy.stats import beta
from scipy.special import betaln
from typing import Dict, Generator, Optional
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend suitable for cluster/terminal
import matplotlib.pyplot as plt


GLOBAL_BLACKLIST_SITES_SET = set()
_GLOBAL_BLACKLISTED_REGIONS = [
    # blacklisted regions (chr, start, end)
    ('chrM', 299, 317),
    ('chrM', 511, 525),
    ('chrM', 564, 571),
    ('chrM', 952, 955),
    ('chrM', 3106, 3108),
    ('chrM', 5895, 5899),
    ('chrM', 8268, 8279),
    ('chrM', 13645, 13650),
    ('chrM', 16180, 16187),
]
for chrom, start, end in _GLOBAL_BLACKLISTED_REGIONS:
    for pos in range(start, end + 1):
        GLOBAL_BLACKLIST_SITES_SET.add((chrom, pos))

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
                'alt': record.alts,
                'qual': record.qual,
                'filter': list(record.filter.keys()),
                'info': dict(record.info),
                'format': record.format,
            }
            
            # Add samples information if present
            if record.samples:
                # Get all sample data: sample => record.samples[sample]
                variant_info['samples'] = {sample: record.samples[sample]  
                                           for sample in record.samples}

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
        # vcf_in.header.info['KL_DIV'] = ('1', 'Float', 'Kullback-Leibler divergence of multi-sample VAF distribution from background noise distribution')
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
                if qc_status.get(chrom_pos, False):
                    record.filter.add('PASS')
                elif chrom_pos in GLOBAL_BLACKLIST_SITES_SET:
                    record.filter.add('BLACKLISTED_SITE')
                else:
                    record.filter.add('LOW_QUALITY')
                
                record.qual = pos_kl_div.get(chrom_pos, 0)
                vcf_out.write(record)
        
        # If output VCF is compressed, index it
        if output_vcf_path.endswith('.gz'):
            pysam.tabix_index(output_vcf_path, preset='vcf', force=True)


def qc(input_vcf_path, output_vcf_path, args):
    """Parse VCF file using pysam to extract variant information and store original records."""
    # Only get samples id from input vcffile
    _vcf = pysam.VariantFile(input_vcf_path, 'r')
    samples = list(_vcf.header.samples)
    _vcf.close()
    
    results = []
    fs_threshold  = 1000
    sor_threshold = 1000
    pos_kl_div = {} # Store multi-sample KL divergence for each position
    for variant in variant_generator(input_vcf_path):
        chrom_pos = (variant['chrom'], variant['pos'])
        if chrom_pos in GLOBAL_BLACKLIST_SITES_SET:
            sys.stderr.write(f"[INFO] Variant at {variant['chrom']}:{variant['pos']} "
                             f"is located in blacklisted regions. Skipping.\n")
            continue
        
        # Maximum ALT alleles per site
        if len(variant['alt']) > args.max_alt_alleles:
            sys.stderr.write(f"[INFO] Variant at {variant['chrom']}:{variant['pos']} has "
                             f"more than {args.max_alt_alleles} ALT alleles. Skipping.\n")
            continue

        # Process each variant as needed
        ref_alts = [variant.get('ref', '')] + list(variant.get('alt', []))
        sample_gts = []
        vaf_obs = []
        background_error_rates = []
        for sample in samples:
            gt = list(variant['samples'][sample].get('GT'))
            hf = list(variant['samples'][sample].get('HF'))

            sample_gts.append(gt)
            vaf_obs.append(hf)
            ploidy = len(gt)
            
            # Pre-filtering based on DP and HQ thresholds
            if variant['samples'][sample].get('DP', 0) / ploidy < args.DP_to_ploidy_threshold:
                continue
            
            sample_hq = variant['samples'][sample].get('HQ', (0))
            if any(hq < args.HQ_threshold for hq in sample_hq):
                continue
            
            sample_fs = variant['samples'][sample].get('FS')
            if (sample_fs is not None) and any(fs > fs_threshold for fs in sample_fs):
                continue

            sample_sor = variant['samples'][sample].get('SOR')
            if (sample_sor is not None) and any(sor > sor_threshold for sor in sample_sor):
                continue

            if (ploidy == 1) and (hf is not None) and all(h is not None for h in hf): 
                # Only use non mutated (homozygous) samples to estimate background error rate.
                background_error_rates.append(1.0 - sum(hf))
        
        if len(background_error_rates) == 0:
            sys.stderr.write(f"[WARNING] All heterozygous samples. No background error "
                             f"rate calculated for this position: "
                             f"{variant['chrom']}:{variant['pos']}. Ignore.\n")
            continue

        # Estimate background noise distribution using normal (homozygous) samples error rates
        q_alpha, q_beta, bin_edges = estimate_background_noise(background_error_rates)

        # Calculate mean background error rate
        p_error = np.mean(background_error_rates)

        # QC for each sample
        sys.stderr.write(f"Processing variant at {variant['chrom']}:{variant['pos']}. "
                         f"Background Error Rate: {p_error}. \n")
        for sample, gt, vaf in zip(samples, sample_gts, vaf_obs):
            sys.stderr.write(f" - {sample}, VAF: {vaf}\n")
            
            # Calculate KL divergence for each sample for each ALT
            kl_div_singles = []
            dp   = variant['samples'][sample].get('DP')
            hqs  = variant['samples'][sample].get('HQ')   # A tuple
            ads  = variant['samples'][sample].get('AD')   # A tuple
            fss  = variant['samples'][sample].get('FS')   # A tuple
            sors = variant['samples'][sample].get('SOR')  # A tuple
            
            # Pre-filtering based on DP and HQ thresholds
            is_pre_filtered = False
            ploidy = len(gt)
            if dp / ploidy < args.DP_to_ploidy_threshold:
                sys.stderr.write(f"[WARNING] Sample {sample} at {variant['chrom']}:{variant['pos']} failed DP/Ploidy "
                                 f"threshold ({dp}/{ploidy} < {args.DP_to_ploidy_threshold}). Skipping.\n")
                pp = 0.0
                is_pre_filtered = True
            elif any(hq < args.HQ_threshold for hq in hqs):
                sys.stderr.write(f"[WARNING] Sample {sample} at {variant['chrom']}:{variant['pos']} "
                                 f"failed HQ threshold ({hqs} < {args.HQ_threshold}). Skipping.\n")
                pp = 0.0
                is_pre_filtered = True
            elif (fss is not None) and any(fs > fs_threshold for fs in fss):
                sys.stderr.write(f"[WARNING] Sample {sample} at {variant['chrom']}:{variant['pos']} "
                                 f"failed FS threshold ({fss} > {fs_threshold}). Skipping.\n")
                pp = 0.0
                is_pre_filtered = True
            elif (sors is not None) and any(sor > sor_threshold for sor in sors):
                sys.stderr.write(f"[WARNING] Sample {sample} at {variant['chrom']}:{variant['pos']} "
                                 f"failed SOR threshold ({sors} > {sor_threshold}). Skipping.\n")
                pp = 0.0
                is_pre_filtered = True
            else :
                # Calculate posterior probability for each sample
                pp = 1.0 # Initialize posterior probability
                for v_obs, a_obs in zip(vaf, ads):
                    if (v_obs is None) or (a_obs is None) or (dp is None):
                        pp = 0
                        sys.stderr.write(f"[WARNING] Missing VAF, AD, or DP for sample {sample} "
                                         f"at {variant['chrom']}:{variant['pos']}. Skipping.\n")
                        continue
                    
                    kl_div = calculate_kl_divergence_single(v_obs, q_alpha, q_beta)
                    kl_div_singles.append(kl_div)
                    
                    # Bayesian filter to compute posterior probability of true mutation
                    posterior = bayesian_filter(
                        A=a_obs,
                        D=dp,
                        p_error=p_error,
                        alpha_h1=1,  # Use initial Beta(1,1), will be updated later in ``iterative_beta_fit_and_call``
                        beta_h1=1,   # Use initial Beta(1,1), will be updated later in ``iterative_beta_fit_and_call``
                        lambda_kl=args.lambda_kl,
                        kl_div=kl_div,
                        pi=args.pi
                    )
                    pp *= posterior

            results.append({
                'sample': sample,
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': [ref_alts[g_idx] if (g_idx is not None) else '.' for g_idx in gt] if gt else ['.'],
                'gt': gt if gt and (not is_pre_filtered and pp > args.threshold) else [None],  # A GT tuple
                
                'vaf': vaf,
                'A': ads,
                'D': dp,
                'p_error': p_error,
                
                'ploidy': len(gt),
                'kl_divergence_single': kl_div_singles,
                'posterior': pp,
                'is_mutation': pp > args.threshold,  # Initial mutation call based on threshold
                'pre_filtered': is_pre_filtered
            })
        
        # Calculate multi-sample KL divergence for the position
        kl_div_multi = calculate_kl_divergence_multi(vaf_obs, q_alpha, q_beta, bin_edges)
        pos_kl_div[(variant['chrom'], variant['pos'])] = kl_div_multi if np.isfinite(kl_div_multi) else 10000

    results, alpha_h1, beta_h1, alpha_hist, beta_hist, diff_hist = iterative_beta_fit_and_call(
        results, lambda_kl=args.lambda_kl, pi=args.pi, threshold=args.threshold
    )

    sample_variant_count = {
        'all': [0 for _ in samples], 
        'mutations': [0 for _ in samples]
    }
    sample_index = {s: i for i, s in enumerate(samples)}
    vaf_true_list = []
    vaf_false_list = []
    for r in results:
        if any(alt != '.' and alt != r['ref'] for alt in r['alt']):  # raw variant exists
            sample_variant_count['all'][sample_index[r['sample']]] += 1
            
        if r['is_mutation']:
            vaf_true_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
            if any(gt for gt in r['gt']):  # gt is not None and gt != 0
                sample_variant_count['mutations'][sample_index[r['sample']]] += 1
        else:
            vaf_false_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])

        r['kl_divergence_single'] = np.mean(r['kl_divergence_single']) if r['kl_divergence_single'] else 0
        r.pop('vaf', None)  # Remove VAF from final results
        r.pop('A', None)    # Remove A from final results
        r.pop('D', None)    # Remove D from final results

    print(f"Total variants processed: {len(results)}, Beta distribution for true variants: Beta({alpha_h1}, {beta_h1}).\n")
    write_vcf(input_vcf_path, output_vcf_path, results, pos_kl_div)

    return results, sample_variant_count, vaf_true_list, vaf_false_list, alpha_hist, beta_hist, diff_hist


def iterative_beta_fit_and_call(results, lambda_kl=0.1, pi=5e-8 * 16569, threshold=0.9, max_iter=50, tol=1e-5):
    """
    Iteratively estimate beta distribution parameters from called mutations and update mutation calls
    until convergence. For multi-allelic sites, treat each ALT as an independent event and use the
    product of posteriors as the final posterior for the site.
    """
    safe_alpha_h1, safe_beta_h1 = 1, 1  # Default values for Beta distribution parameters
    alpha_hist, beta_hist, diff_hist = [1], [1], []
    prev_is_mut = [r['is_mutation'] for r in results]
    for i in range(max_iter):
        # 1. Collect VAFs from currently called mutations
        het_var_list = []
        hom_var_list = []
        for r in results:
            if r['is_mutation'] and (not r['pre_filtered']) and (r.get('vaf') is not None):
                if r['ploidy'] > 1:
                    het_var_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
                else:
                    hom_var_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
                
        by_step = len(hom_var_list) // len(het_var_list) if len(hom_var_list) > len(het_var_list) else 1
        vaf_list = het_var_list + hom_var_list[::by_step]  # Use same number of het and hom variants, by downsampling hom variants by step
        print(f"- Iteration {i+1}: {len(vaf_list)} VAFs collected for Beta fitting.")

        # 2. Fit Beta distribution parameters using MLE
        alpha_h1, beta_h1 = estimate_beta_params_mle(vaf_list)

        # Ensure alpha_h1 and beta_h1 are valid
        safe_alpha_h1 = alpha_h1 if alpha_h1 is not None else 1
        safe_beta_h1 = beta_h1 if beta_h1 is not None else 1
        alpha_hist.append(safe_alpha_h1)
        beta_hist.append(safe_beta_h1)
        
        # 3. Recalculate posterior and is_mutation for each record
        for r in results:
            if r['pre_filtered']:
                continue  # Skip pre-filtered samples
            
            # For multi-allelic sites, calculate posterior for each ALT and take the product
            new_pps = []
            for a_obs, kl_div in zip(r['A'], r['kl_divergence_single']):
                new_pp = bayesian_filter(
                    A=a_obs,
                    D=r['D'],
                    p_error=r['p_error'],
                    alpha_h1=safe_alpha_h1,  # Use updated Beta parameters
                    beta_h1=safe_beta_h1,    # Use updated Beta parameters
                    lambda_kl=lambda_kl,
                    kl_div=kl_div,
                    pi=pi,
                )
                new_pps.append(new_pp)

            r['posterior'] = np.prod(new_pps)  # Update posterior as product of individual posteriors
            r['is_mutation'] = r['posterior'] > threshold  # Update mutation call based on new posterior
            r['gt'] = r['gt'] if r['is_mutation'] else [None]
            
        # 4. Check for convergence
        curr_is_mut = [r['is_mutation'] for r in results]
        if prev_is_mut is not None:
            diff = np.mean([a != b for a, b in zip(prev_is_mut, curr_is_mut)])
            diff_hist.append(diff)
            if diff < tol:
                break
            
        prev_is_mut = curr_is_mut.copy()

    return results, safe_alpha_h1, safe_beta_h1, alpha_hist, beta_hist, diff_hist


def estimate_background_noise(background_error_rates, bins=100):
    """Estimate background noise distribution using normal samples (GT='0', '1', '2', ...)."""
    
    # Fit Beta distribution
    if len(background_error_rates) == 0 or np.all(background_error_rates == 0):
        return 1e-4, 1e-8, np.linspace(0, 1, bins + 1)  #, np.zeros(bins), np.linspace(0, 1, bins + 1)
    
    hist, bin_edges = np.histogram(background_error_rates, bins=bins, density=True, range=(0, 1))
    # a, b = parametrize_vaf_distribution(background_error_rates)
    a, b = estimate_beta_params_mle(background_error_rates)

    return a, b, bin_edges  # return alpha and beta parameters for Beta distribution


def estimate_beta_params_mle(vaf_list):
    """Fit Beta distribution parameters using Maximum Likelihood Estimation (MLE)."""
    # Remove None and extreme values to avoid fitting errors
    vaf_arr = np.array([v for v in vaf_list if v is not None and 0 < v < 1])
    if len(vaf_arr) < 2:
        return 1, 1  # Return default if not enough data
        
    # Use scipy's beta.fit with fixed loc=0, scale=1
    a, b, loc, scale = beta.fit(vaf_arr, floc=0, fscale=1)
    return a, b


def parametrize_vaf_distribution(vaf_samples):
    """Fit a Beta distribution to the VAF samples."""
    if not vaf_samples:
        return 1e-4, 1e-8

    mean_vaf = np.mean(vaf_samples)
    var_vaf = np.var(vaf_samples) if len(vaf_samples) > 1 else mean_vaf * (1 - mean_vaf)
    if var_vaf == 0:
        var_vaf = 1e-6  # Prevent division by zero
        
    w = mean_vaf * (1 - mean_vaf) / var_vaf - 1

    # `alpha` and `beta` parameters for Beta distribution
    a = max(mean_vaf * w, 1e-4)  # Ensure positive values
    b = max((1 - mean_vaf) * w, 1e-8)
    
    return a, b


def calculate_kl_divergence_single(v_obs, a, b):
    """Calculate KL divergence for a single sample VAF."""
    if v_obs is None or v_obs <= 0 or v_obs >= 1:
        return 0  # Invalid VAF value

    q_v = beta.pdf(v_obs, a, b)
    sys.stderr.write(f" - Q({v_obs:.4f}|{a:.4f}, {b:.4f}) = {q_v}\n")
    return -np.log(q_v) if q_v > 0 else np.inf


def calculate_kl_divergence_multi(vafs, q_alpha, q_beta, bin_edges):
    """Calculate KL divergence for multi-sample VAF distribution."""
    vafs_flatten = []
    for vaf in vafs:
        vafs_flatten.extend([x for x in vaf if x is not None])

    if len(vafs_flatten) == 0:
        return 0
        
    px_at_bins, _ = np.histogram(vafs_flatten, bins=bin_edges, density=True)  # p_x
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    qx_at_bins = beta.pdf(bin_centers, q_alpha, q_beta)
    qx_at_bins = np.where(qx_at_bins > 0, qx_at_bins, 1e-10)  # Avoid division by zero

    delta_bin_width = bin_edges[1] - bin_edges[0]
    kl_div = np.sum(px_at_bins * np.log(px_at_bins / qx_at_bins + 1e-16) * delta_bin_width)

    return kl_div if kl_div > 0 else 0

# Alternative implementation using histogram for KL divergence calculation (Not be used for final version)
# def calculate_kl_divergence_multi(vafs, noise_alpha, noise_beta, bins=100):
#     """Calculate KL divergence for multi-sample VAF distribution."""
#     if len(vafs) == 0:
#         return 0
    
#     hist_p, bin_edges_p = np.histogram(vafs, bins=bins, density=True, range=(0, 1))
#     kl_div = 0
#     for i in range(len(hist_p)):
#         p_x = hist_p[i] if hist_p[i] > 0 else 1e-10
#         bin_center = (bin_edges_p[i] + bin_edges_p[i+1]) / 2
#         q_x = beta.pdf(bin_center, noise_alpha, noise_beta) if beta.pdf(bin_center, noise_alpha, noise_beta) > 0 else 1e-10
#         kl_div += p_x * np.log(p_x / q_x + 1e-16)
    
#     return kl_div if kl_div > 0 else 0


def bayesian_filter(A, D, p_error, alpha_h1=1, beta_h1=1, lambda_kl=0.1, 
                    kl_div=None, pi=5e-8 * 16569):
    """Bayesian filter to compute posterior probability of a true mutation using log-likelihoods.
    """
    if kl_div == np.inf or p_error == 0.0:
        return 1 # If KL divergence is infinite, return 1 (indicating mutation)
    
    # Log-likelihood under H0 (without binomial coefficient)
    log_L0 = A * np.log(p_error) + (D - A) * np.log(1 - p_error)
    
    # Log-likelihood under H1 (without binomial coefficient)
    log_L1 = betaln(A + alpha_h1, D - A + beta_h1) - betaln(alpha_h1, beta_h1)
     
    # Adjust log_L1 with KL divergence (which is a measure of how much the observed distribution diverges from the expected background noise)
    kl_div = kl_div if kl_div is not None else 0
    log_L1_adjusted = log_L1 + lambda_kl * kl_div

    # Log-posterior odds
    log_odds = log_L1_adjusted - log_L0 + np.log(pi) - np.log(1.0 - pi)
    
    # Posterior probability
    posterior_h1 = 1 / (1 + np.exp(-log_odds))
    
    return posterior_h1


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
                        help="Minimum depth (DP/ploidy) threshold for considering variants "
                             "for each sample. Default is 100")
    parser.add_argument('--HQ-threshold', type=int, default=20, 
                        help="Minimum base quality (HQ) threshold for considering "
                             "variants for each sample. Default is 20")
    
    # For Bayesian filter parameters
    parser.add_argument('--bins', type=int, default=100, help="Number of histogram bins. Default is 100")
    parser.add_argument('--lambda-kl', type=float, default=0.1, help="Weight for KL divergence. Default is 0.1")
    parser.add_argument('--pi', type=float, default=5e-8 * 16569, help="Prior probability of mutation. Default is 5e-8 * 16569")
    parser.add_argument('--threshold', type=float, default=0.9, help="Posterior probability threshold for calling mutations. Default is 0.9")
    
    args = parser.parse_args()
    try:
        # Parse VCF file
        (variants, sample_variant_count, vaf_true_list, vaf_false_list, 
         alpha_hist, beta_hist, diff_hist) = qc(args.vcf, args.output_vcf, args)
        
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
                                  save_path=args.output.split('.')[0] + '_beta_fit_convergence.png')
        print(f'diff_hist: {diff_hist}\nalpha_hist: {alpha_hist}\nbeta_hist: {beta_hist}')
        
        # plot Beta distribution fit and VAFs
        plot_beta_fit_and_vaf(alpha_hist, beta_hist, vaf_true_list, vaf_false_list, 
                              save_path=args.output.split('.')[0] + '_beta_fit_vaf.png')
        print(f"Quality control analysis completed. Results saved to {args.output_vcf}")
        
        # Plot variant counts per sample
        plot_variant_counts(sample_variant_count, 
                            save_path=args.output.split('.')[0] + '_variant_counts_per_sample.png')
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    