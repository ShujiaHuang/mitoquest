#!/usr/bin/env python3
"""Detect NUMT bias by using VAF of mtDNA variants and copy number of mtDNA (mtCN).

Authors: Shujia Huang
Date: 2026-04-24
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats


def numt_match(df1, df2):
    # 1. Sort
    df1 = df1.sort_values('Pos')
    df2 = df2[df2['is_NUMT']].sort_values('Pos')  # 只挑 NUMT 风险位点
    
    # 2. region
    df1_start = df1['Pos'].values
    df1_end   = df1_start + df1['REF'].astype(str).str.len().values - 1
    
    df2_pos = df2['Pos'].values
    df2_numt = df2['is_NUMT'].values
    if len(df2_pos) == 0:
        df1['is_NUMT'] = False
        return df1

    # Creat an initial array as the size as df1 and set to False value
    numt_record = np.full(len(df1), False)
    
    idx, df2_size = 0, len(df2_pos)
    for i in range(len(df1)):
        is_overlap = False    
        for j in range(idx, df2_size):
            if df1_start[i] > df2_pos[j]: continue
            if df1_end[i] < df2_pos[j]: break

            idx = j
            is_overlap = True
            break
            
        if is_overlap:
            numt_record[i] = df2_numt[idx]
    
    df1['is_NUMT'] = numt_record
    return df1


def detect_numt_artifacts_by_copynumber(df, 
                                        pos_col='Pos',
                                        ref_col='REF',
                                        alt_col='ALT',
                                        vaf_col='VAF',
                                        copynum_col='copynum', 
                                        min_samples=15, 
                                        pval_threshold=0.05, 
                                        corr_threshold=-0.4):
    """
    通过 VAF 与真实线粒体拷贝数 (Copy Number) 的负相关性，标记 NUMT bias。
    兼容队列数据中存在混合测序深度 (如 7X, 15X, 30X 混合) 的场景。
    
    参数:
    df: 
    min_samples: 最少样本数。
    pval_threshold: Spearman 负相关的显著性 P 值阈值。
    corr_threshold: 负相关阈值，越接近 -1 越代表是典型的 NUMT。
    
    返回:
    带有 'is_NUMT' 的 df 结果
    """
    required_cols = {pos_col, ref_col, alt_col, vaf_col, copynum_col}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")
    
    # 过滤掉零数据
    working_df = df[df[copynum_col] > 0].copy()
    
    # 将 CN 转换到对数空间，使 NUMT ~ 1/(CN+1) 的非线性关系趋于线性，提升灵敏度
    # 该关系的数学推导在 PPT 中记录
    working_df['log_CN'] = np.log10((working_df[copynum_col]+1))
    numt_results = []
    
    # 2. 按位点和突变类型逐一进行分组
    grouped = working_df.groupby(['Pos', 'ALT'])
    for (locus, allele), group in grouped:
        n_samples = len(group)
        if n_samples < min_samples:
            numt_results.append({
                'Pos': locus, 
                'ALT': allele,
                'NUMT_Corr': np.nan, 
                'NUMT_Pval': np.nan, 
                'is_NUMT': False
            })
            continue
            
        # 3. 使用 Spearman 而不是 Pearson，更适合处理离群样本
        result = stats.spearmanr(group['log_CN'], group['VAF'])
        corr = float(result.correlation)
        pval = float(result.pvalue)
        
        # 4. 是否显著且足够负相关
        is_numt = False
        if (corr < corr_threshold) and (pval < pval_threshold):
            is_numt = True
            
        numt_results.append({
            'Pos': locus, 
            'ALT': allele,
            'NUMT_Corr': corr, 
            'NUMT_Pval': pval, 
            'is_NUMT': is_numt
        })
        
    numt_df = pd.DataFrame(numt_results)
    final_df = numt_match(working_df, numt_df[numt_df['is_NUMT']])
    detected_numts = final_df[final_df['is_NUMT']][['Pos', 'ALT']].drop_duplicates()
    if not detected_numts.empty:
        print(f"Successfully detected {len(detected_numts)} NUMT artifacts driven by Copy Number.")
        
    return final_df, numt_df


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Detect NUMT bias by using VAF of mtDNA variants and copy number of mtDNA (mtCN).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input TSV file (contains: Pos, REF,ALT,VAF, and copy-number columns)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output tidy TSV file",
    )
    parser.add_argument(
        '-P', '--pos', type=str, default='Pos', 
        help="Name of POS column. Default: Pos"
    )
    parser.add_argument(
        '-R', '--ref', type=str, default='REF', 
        help="Name of REF column. Default: REF"
    )
    parser.add_argument(
        '-A', '--alt', type=str, default='ALT', 
        help="Name of ALT column. Default: ALT"
    )
    parser.add_argument(
        '-V', '--vaf', type=str, default='VAF', 
        help="Name of VAF column. Default: VAF"
    )
    parser.add_argument(
        '-C', '--copynum', type=str, default='VAF', 
        help="Name of MT copynumber column. Default: copynum"
    )
    parser.add_argument(
        '--min-samples', type=int, default=15, 
        help="minimum sample size for detecting the rise of NUMT. Default: 15"
    )
    parser.add_argument(
        '--pvalue-threshold', type=float, default=0.05, 
        help="Pvalue threshold for NUMT detection. Default: 0.05"
    )
    parser.add_argument(
        '--corr-threshold', type=float, default=-0.4, 
        help="Spearmanr Correlation threshold for NUMT detection. Default: -0.4"
    )
    
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input, sep='\t')
    final_df, numt_df = detect_numt_artifacts_by_copynumber(
        df, 
        pos_col=args.pos,
        ref_col=args.ref,
        alt_col=args.alt,
        vaf_col=args.vaf,
        copynum_col=args.copynum,
        min_samples=args.min_samples, 
        pval_threshold=args.pvalue_threshold, 
        corr_threshold=args.corr_threshold
    )
    final_df.to_csv(args.output, sep='\t', index=False)
    

if __name__ == "__main__":
    main()
