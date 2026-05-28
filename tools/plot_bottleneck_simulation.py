#!/usr/bin/env python3
"""Plot the deCODE-style "Figure 5" bottleneck-parameter simulation.

Reproduces the per-maternal-VAF-bin "Observed and simulated means for
bottleneck parameter, b" panel from Helgason et al. 2024 Cell using the
TSV emitted by `mitoquest ne-estimate --bin-simulation`.

The TSV stores, per equal-width maternal-VAF bin:

    obs_var       =  empirical mean of (p_c - p_m)^2          (raw drift)
    obs_var_corr  =  empirical mean of (d_i - s_i)            (sampling-corrected)
    obs_F         =  empirical mean of (d_i - s_i) / [p_m(1 - p_m)]
                                                              (= 1 - b per bin)
    obs_F_se      =  standard error of obs_F in the bin

and the theoretical curves p_m(1 - p_m) / Ne and 1/Ne at the fitted Ne
plus its 95% profile-likelihood CI bounds.  This script overlays the
observed bin means (with error bars) on those theoretical curves and
optionally annotates the Wonnapinij/Kimura cross-check Ne.

Usage:
    python plot_bottleneck_simulation.py --input cohort.bin_sim.tsv \
                                         --output cohort.bottleneck.png

Author: Shujia Huang
Date: 2026-05-28
"""
import argparse
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend suitable for cluster/terminal
import matplotlib.pyplot as plt


def parse_header_metadata(tsv_path):
    """Read the leading commented metadata of the bin-simulation TSV.

    Returns a dict with keys such as 'fitted_ne', 'fitted_ne_ci_low',
    'fitted_ne_ci_high', 'kimura_ne', 'kimura_b', 'model', 'n_pairs_used'.
    Missing keys are absent from the returned dict.
    """
    meta = {}
    with open(tsv_path, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            stripped = line.strip()
            if '=' not in stripped:
                continue
            # The very first line is the verbatim command; ignore it here
            # because there is no key=value pair worth keeping.
            if stripped.startswith('#mitoquest_ne_estimate_command='):
                continue
            key, _, value = stripped.lstrip('#').partition('=')
            key = key.strip()
            value = value.strip()
            try:
                meta[key] = float(value)
            except ValueError:
                meta[key] = value
    return meta


def load_bin_simulation(tsv_path):
    """Load the per-bin observed-vs-simulated drift table.

    Returns:
        meta:  dict of header metadata (see parse_header_metadata).
        df:    DataFrame with one row per bin.
    """
    meta = parse_header_metadata(tsv_path)
    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    return meta, df


def plot_drift_panel(ax, df, meta):
    """Draw the deCODE-style raw drift panel (observed vs simulated)."""
    p = np.linspace(0.0, 1.0, 201)
    ne_pt = meta.get('fitted_ne')
    ne_lo = meta.get('fitted_ne_ci_low')
    ne_hi = meta.get('fitted_ne_ci_high')

    # Observed bin means; use the standard error of mean F_i as an
    # approximate error bar on (p_m(1 - p_m)) * F_i.
    pbar = df['bin_center'].values
    obs = df['obs_var'].values
    obs_err = df['obs_F_se'].values * pbar * (1.0 - pbar)
    sizes = 30 + 4 * df['n_pairs'].values  # marker scales with N

    # Theoretical curves at fitted Ne and CI bounds.  Smaller Ne ->
    # higher drift, so the upper ribbon edge corresponds to Ne_low and
    # the lower edge to Ne_high.
    if ne_pt and ne_pt > 0:
        ax.plot(p, p * (1.0 - p) / ne_pt, color='#d62728', lw=2.0,
                label=f'Simulated  p(1-p)/Ne   [Ne_MLE = {ne_pt:.2f}]')
    if ne_lo and ne_hi and ne_lo > 0 and ne_hi > 0:
        upper = p * (1.0 - p) / ne_lo
        lower = p * (1.0 - p) / ne_hi
        ax.fill_between(p, lower, upper, color='#d62728', alpha=0.15,
                        label=f'95% CI  [{ne_lo:.2f}, {ne_hi:.2f}]')

    # Optional Wonnapinij/Kimura overlay.
    ne_kim = meta.get('kimura_ne')
    if ne_kim and ne_kim > 0:
        ax.plot(p, p * (1.0 - p) / ne_kim, color='#1f77b4', lw=1.5,
                ls='--', label=f'Kimura  p(1-p)/Ne   [Ne_Kimura = {ne_kim:.2f}]')

    ax.scatter(pbar, obs, s=sizes, color='black', zorder=4,
               label='Observed bin mean (size ~ N pairs)')
    ax.errorbar(pbar, obs, yerr=obs_err, fmt='none', ecolor='black',
                elinewidth=1.0, capsize=3, zorder=3)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel('Maternal heteroplasmy   $p_m$  (bin center)')
    ax.set_ylabel(r'Mean  $(p_c - p_m)^2$  per bin')
    ax.set_title('Observed vs simulated drift   $E[(p_c - p_m)^2] = p_m(1-p_m) / N_e$')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(loc='upper right', fontsize=9)


def plot_F_panel(ax, df, meta):
    """Draw the per-bin (1 - b) panel: F_i = (d - s) / [p(1-p)] vs p_m."""
    pbar = df['bin_center'].values
    obsF = df['obs_F'].values
    err = df['obs_F_se'].values
    sizes = 30 + 4 * df['n_pairs'].values

    ne_pt = meta.get('fitted_ne')
    ne_lo = meta.get('fitted_ne_ci_low')
    ne_hi = meta.get('fitted_ne_ci_high')
    ne_kim = meta.get('kimura_ne')

    if ne_pt and ne_pt > 0:
        ax.axhline(1.0 / ne_pt, color='#d62728', lw=2.0,
                   label=f'Simulated  1/Ne   [Ne_MLE = {ne_pt:.2f}]')
    if ne_lo and ne_hi and ne_lo > 0 and ne_hi > 0:
        ax.axhspan(1.0 / ne_hi, 1.0 / ne_lo, color='#d62728', alpha=0.15,
                   label=f'95% CI  [{ne_lo:.2f}, {ne_hi:.2f}]')
    if ne_kim and ne_kim > 0:
        ax.axhline(1.0 / ne_kim, color='#1f77b4', lw=1.5, ls='--',
                   label=f'Kimura  1/Ne   [Ne_Kimura = {ne_kim:.2f}]')

    ax.scatter(pbar, obsF, s=sizes, color='black', zorder=4,
               label='Observed bin mean F_i (size ~ N pairs)')
    ax.errorbar(pbar, obsF, yerr=err, fmt='none', ecolor='black',
                elinewidth=1.0, capsize=3, zorder=3)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel('Maternal heteroplasmy   $p_m$  (bin center)')
    ax.set_ylabel(r'$F_i = (d_i - s_i) / [p_m(1-p_m)]$   ($= 1 - b$)')
    ax.set_title('Per-bin bottleneck parameter   $1 - b = 1 / N_e$')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(loc='upper right', fontsize=9)


def make_plot(tsv_path, output_path, dpi=300, figsize=(13, 5.2), title=None):
    """Render the two-panel deCODE-style figure to `output_path`."""
    meta, df = load_bin_simulation(tsv_path)
    if df.empty:
        raise ValueError(f'[plot] No bin rows found in {tsv_path}')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    plot_drift_panel(ax1, df, meta)
    plot_F_panel(ax2, df, meta)

    n_pairs = int(meta.get('n_pairs_used', df['n_pairs'].sum()))
    model = meta.get('model', 'continuous')
    suptitle = title or (
        f'mitoquest ne-estimate  --  bottleneck simulation '
        f'({n_pairs} pairs, model={model})'
    )
    fig.suptitle(suptitle, y=1.02, fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'[plot] Wrote {output_path}', file=sys.stderr)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot the deCODE-style observed-vs-simulated bottleneck "
                    "drift panel from `mitoquest ne-estimate --bin-simulation`.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Bin-simulation TSV produced by `mitoquest ne-estimate '
             '--bin-simulation`.',
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output figure path (PNG or PDF, inferred from extension).',
    )
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Output DPI for raster formats.',
    )
    parser.add_argument(
        '--title', type=str, default=None,
        help='Optional figure suptitle override.',
    )
    parser.add_argument(
        '--figsize', type=str, default='13,5.2',
        help='Figure size in inches as "width,height".',
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        w, h = (float(x) for x in args.figsize.split(','))
    except ValueError:
        raise SystemExit(f'[plot] Invalid --figsize: {args.figsize!r} '
                         '(expected "width,height").')
    make_plot(args.input, args.output, dpi=args.dpi,
              figsize=(w, h), title=args.title)


if __name__ == "__main__":
    main()
