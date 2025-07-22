
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
    # Build lookup: {(chrom, pos, sample): (pp, is_mutation)}
    result_lookup = {}
    for r in results:
        key = (r['chrom'], r['pos'], r['sample'])
        result_lookup[key] = (r['posterior'], int(r['is_mutation']))

    df_results = pd.DataFrame(results)
    qc_status = df_results.groupby(['chrom', 'pos'])['is_mutation'].max().reset_index().set_index(['chrom', 'pos']).to_dict()['is_mutation']
    with pysam.VariantFile(input_vcf_path) as vcf_in:

        # Add new FORMAT fields to header
        vcf_in.header.add_line('##FORMAT=<ID=PP,Number=1,Type=Float,Description="Posterior probability of true mutation">')
        vcf_in.header.add_line('##FORMAT=<ID=IS_MUT,Number=1,Type=String,Description="Mutation call (True=mutation, False=not)">')
        vcf_in.header.add_line('##FILTER=<ID=QC_FAIL,Description="Failed mtDNA QC filter">')
        
        with pysam.VariantFile(output_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                for sample in record.samples:
                    key = (record.chrom, record.pos, sample)
                    if key in result_lookup:
                        pp, is_mut = result_lookup[key]
                        record.samples[sample].update({'PP': float(pp), 'IS_MUT': 'True' if is_mut else 'False'})

                # Fix: use record.filter.clear() and record.filter.add() to update FILTER.
                record.filter.clear()
                chrom_pos = (record.chrom, record.pos)
                if qc_status.get(chrom_pos, False):
                    record.filter.add('PASS')
                else:
                    record.filter.add('QC_FAIL')
                
                record.qual = pos_kl_div.get(chrom_pos, 0)
                vcf_out.write(record)
        
        # If output VCF is compressed, index it
        if output_vcf_path.endswith('.gz'):
            pysam.tabix_index(output_vcf_path, preset='vcf', force=True)


def qc(input_vcf_path, output_vcf_path, bins=100, lambda_kl=0.1, pi=5e-8 * 16569, threshold=0.9):
    """Parse VCF file using pysam to extract variant information and store original records."""
    # Only get samples id from input vcffile
    _vcf = pysam.VariantFile(input_vcf_path, 'r')
    samples = list(_vcf.header.samples)
    _vcf.close()
    
    results = []
    pos_kl_div = {} # Store multi-sample KL divergence for each position
    for variant in variant_generator(input_vcf_path):
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
            if (ploidy == 1) and (hf is not None) and all(h is not None for h in hf): 
                # Only use non mutated (homozygous) samples to 
                # estimate background error rate.
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
        sys.stderr.write(f"Processing variant at {variant['chrom']}:{variant['pos']}. Background Error Rate: {p_error}. \n")
        for sample, gt, vaf in zip(samples, sample_gts, vaf_obs):
            sys.stderr.write(f" - {sample}, VAF: {vaf}\n")
            
            # Calculate KL divergence for each sample for each ALT
            dp = variant['samples'][sample].get('DP')
            ad = list(variant['samples'][sample].get('AD'))
            
            pp = 1.0 # Initialize posterior probability
            kl_div_singles = []
            for v_obs, a_obs in zip(vaf, ad):
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
                    alpha_h1=1,
                    beta_h1=1,
                    lambda_kl=lambda_kl,
                    kl_div=kl_div,
                    pi=pi
                )
                pp *= posterior

            results.append({
                'sample': sample,
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': [ref_alts[g_idx] if g_idx is not None else '.' for g_idx in gt] if gt else ['.'],
                
                'vaf': vaf,
                'A': ad,
                'D': dp,
                'p_error': p_error,
                
                'ploidy': len(gt),
                'kl_divergence_single': kl_div_singles,
                'posterior': pp,
                'is_mutation': pp > threshold
            })
        
        # Calculate multi-sample KL divergence for the position
        kl_div_multi = calculate_kl_divergence_multi(vaf_obs, q_alpha, q_beta, bin_edges)
        pos_kl_div[(variant['chrom'], variant['pos'])] = kl_div_multi if np.isfinite(kl_div_multi) else 10000

    results, alpha_h1, beta_h1, alpha_hist, beta_hist, diff_hist = iterative_beta_fit_and_call(results, lambda_kl=lambda_kl, pi=pi, threshold=threshold)
    for r in results:
        r.pop('vaf', None)  # Remove VAF from final results
        r.pop('A', None)    # Remove A from final results
        r.pop('D', None)    # Remove D from final results
        r['kl_divergence_single'] = np.mean(r['kl_divergence_single']) if r['kl_divergence_single'] else 0

    print(f"Total variants processed: {len(results)}, Beta distribution for true variants: Beta({alpha_h1}, {beta_h1}).\n")
    write_vcf(input_vcf_path, output_vcf_path, results, pos_kl_div)
    
    return results, alpha_hist, beta_hist, diff_hist


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
            if r['is_mutation'] and r.get('vaf') is not None:
                if r['ploidy'] > 1:
                    het_var_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
                else:
                    hom_var_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
                
                # if isinstance(r['vaf'], list):
                #     vaf_list.extend([v for v in r['vaf'] if v is not None and 0 < v < 1])
                # else:
                #     vaf_list.append(r['vaf'])
                
        # vaf_list = []
        # hom_var_list = np.array(hom_var_list)
        # np.random.shuffle(hom_var_list) # Random shuffling
        # n = min(len(het_var_list), len(hom_var_list)) 
        step = len(hom_var_list) // len(het_var_list) if len(hom_var_list) > len(het_var_list) else 1
        vaf_list = het_var_list + hom_var_list[::step]  # Use same number of het and hom variants
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
            # For multi-allelic sites, calculate posterior for each ALT and take the product
            new_pps = []
            for a_obs, kl_div in zip(r['A'], r['kl_divergence_single']):
                new_pp = bayesian_filter(
                    A=a_obs,
                    D=r['D'],
                    p_error=r['p_error'],
                    alpha_h1=safe_alpha_h1,
                    beta_h1=safe_beta_h1,
                    lambda_kl=lambda_kl,
                    kl_div=kl_div,
                    pi=pi,
                )
                new_pps.append(new_pp)

            r['posterior'] = np.prod(new_pps)
            r['is_mutation'] = r['posterior'] > threshold
            
        # 4. Check for convergence
        curr_is_mut = [r['is_mutation'] for r in results]
        if prev_is_mut is not None:
            diff = np.mean([a != b for a, b in zip(prev_is_mut, curr_is_mut)])
            diff_hist.append(diff)
            if diff < tol:
                break
            
        prev_is_mut = curr_is_mut.copy()

    return results, safe_alpha_h1, safe_beta_h1, alpha_hist, beta_hist, diff_hist

def estimate_beta_params_mle(vaf_list):
    """Fit Beta distribution parameters using Maximum Likelihood Estimation (MLE)."""
    # Remove None and extreme values to avoid fitting errors
    vaf_arr = np.array([v for v in vaf_list if v is not None and 0 < v < 1])
    if len(vaf_arr) < 2:
        return 1, 1  # Return default if not enough data
        
    # Use scipy's beta.fit with fixed loc=0, scale=1
    a, b, loc, scale = beta.fit(vaf_arr, floc=0, fscale=1)
    return a, b


def estimate_background_noise(background_error_rates, bins=100):
    """Estimate background noise distribution using normal samples (GT='0', '1', '2', ...)."""
    
    # Fit Beta distribution
    if len(background_error_rates) == 0 or np.all(background_error_rates == 0):
        return 1e-4, 1e-8, np.linspace(0, 1, bins + 1)  #, np.zeros(bins), np.linspace(0, 1, bins + 1)
    
    hist, bin_edges = np.histogram(background_error_rates, bins=bins, density=True, range=(0, 1))
    a, b = parametrize_vaf_distribution(background_error_rates)

    return a, b, bin_edges  # return alpha and beta parameters for Beta distribution


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
    parser.add_argument('--bins', type=int, default=100, help="Number of histogram bins")
    parser.add_argument('--lambda-kl', type=float, default=0.1, help="Weight for KL divergence")
    parser.add_argument('--pi', type=float, default=5e-8 * 16569, help="Prior probability of mutation")
    parser.add_argument('--threshold', type=float, default=0.9, help="Posterior probability threshold for calling mutations")
    
    args = parser.parse_args()
    try:
        # Parse VCF file
        variants, alpha_hist, beta_hist, diff_hist = qc(args.vcf, args.output_vcf, args.bins, args.lambda_kl, args.pi, args.threshold)
        # Plot convergence of beta parameters and mutation calls
        plot_beta_fit_convergence(alpha_hist, beta_hist, diff_hist, save_path=args.output.split('.')[0] + '_beta_fit_convergence.png')
        print(f'diff_hist: {diff_hist}\nalpha_hist: {alpha_hist}\nbeta_hist: {beta_hist}')

        # Save results to CSV
        df_results = pd.DataFrame(variants).explode('alt')
        df_results = df_results[df_results['ref'] != df_results['alt']]
        df_results.to_csv(args.output, sep='\t', index=False)
        print(f"Quality control analysis completed. Results saved to {args.output_vcf}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    