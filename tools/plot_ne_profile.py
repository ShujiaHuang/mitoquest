#!/usr/bin/env python3
"""Plot the deCODE-style Ne-profile comparison between MMLE and Kimura.

The deCODE 2024 *Cell* paper picked the mtDNA bottleneck size Ne ~~ 3 by
scanning candidate Ne values, simulating the Kimura allele-frequency-
change distribution at each one, and choosing the Ne that best fits the
observed distribution across 137 variants in 53,041 mother-child pairs.

This script reproduces that diagnostic for the data fed to
`mitoquest ne-estimate`.  The companion C++ command writes a TSV that
scores every candidate Ne under both estimators in the program:

    mmle_log_lik(Ne) -- global marginal log-likelihood under the configured
                        model (continuous Beta-diffusion or discrete
                        Beta-Binomial).  Maximised at the fitted Ne_MMLE.

    kimura_ssr(Ne)   -- Sigma_i ((d_i - s_i) - p_m_i (1 - p_m_i) / Ne)^2.
                        Per-pair least-squares fit of the one-generation
                        Wright-Fisher prediction.  Minimised at the analytic
                        Ne_kimura_ssr  =  Sigma w^2 / Sigma rw.

The output figure has two panels:

    Left   --  -2(LL - LL_max) vs Ne, with the chi^2 95% CI threshold (3.84).
               The fitted Ne_MMLE and its 95% profile CI are annotated.

    Right  --  Kimura per-pair SSR vs Ne, with the analytic Ne_kimura_ssr
               and the moment-based Ne_kimura (1 / (1 - b)) annotated, plus
               a deCODE Ne = 3 reference line.

Usage:
    python plot_ne_profile.py --input cohort.ne_profile.tsv \
                              --output cohort.ne_profile.png

Author: Shujia Huang
Date: 2026-05-30
"""
import argparse
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend suitable for cluster/terminal
import matplotlib.pyplot as plt


def parse_header_metadata(tsv_path: str) -> dict:
    """Parse the leading `#key=value` header lines into a dict.

    Lines that do not match the `#key=value` pattern are kept under a
    catch-all `__raw_header__` key so the caller can still inspect e.g. the
    `#mitoquest_ne_estimate_command=...` provenance line.
    """
    meta = {"__raw_header__": []}
    with open(tsv_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            stripped = line.lstrip("#").rstrip("\n")
            meta["__raw_header__"].append(stripped)
            if "=" in stripped:
                key, _, value = stripped.partition("=")
                meta[key.strip()] = value.strip()
    return meta


def load_ne_profile(tsv_path: str):
    """Return (metadata_dict, profile_dataframe).

    Drops rows where mmle_log_lik is non-finite (Ne = 1 under the
    continuous model is a degenerate point: see ne_estimate.cpp), so the
    plot does not blow up on the y-axis.
    """
    meta = parse_header_metadata(tsv_path)
    df = pd.read_csv(tsv_path, sep="\t", comment="#")
    # v1.8.6 renamed mle_* columns to mmle_*; accept the legacy names too.
    if "mmle_log_lik" not in df.columns and "mle_log_lik" in df.columns:
        df = df.rename(columns={"mle_log_lik":   "mmle_log_lik",
                                "mle_delta_2ll": "mmle_delta_2ll"})
    df = df[np.isfinite(df["mmle_log_lik"])].reset_index(drop=True)
    return meta, df


def _meta_float(meta: dict, key: str) -> float:
    """Parse a numeric metadata field; return NaN when missing/unparseable."""
    if key not in meta:
        return float("nan")
    try:
        return float(meta[key])
    except (TypeError, ValueError):
        return float("nan")


def plot_mmle_panel(ax, df: pd.DataFrame, meta: dict) -> None:
    """Left panel: -2(LL - LL_max) vs Ne under the MMLE."""
    ax.plot(df["ne_candidate"], df["mmle_delta_2ll"],
            color="#c0392b", linewidth=2.0, label="MMLE  -2 (LL - LL_max)")

    # 95% profile-likelihood threshold: chi2_{1, 0.95} = 3.841.
    ax.axhline(3.841, color="grey", linestyle="--", linewidth=1.0,
               label=r"$\chi^2_{1,\,0.95}=3.841$")

    # v1.8.6 renamed metadata keys; fall back to the v1.8.5 names.
    ne_mmle = _meta_float(meta, "fitted_ne_mmle")
    if not np.isfinite(ne_mmle):
        ne_mmle = _meta_float(meta, "fitted_ne_mle")
    ne_lo   = _meta_float(meta, "fitted_ne_mmle_ci_low")
    if not np.isfinite(ne_lo):
        ne_lo = _meta_float(meta, "fitted_ne_mle_ci_low")
    ne_hi   = _meta_float(meta, "fitted_ne_mmle_ci_high")
    if not np.isfinite(ne_hi):
        ne_hi = _meta_float(meta, "fitted_ne_mle_ci_high")
    if np.isfinite(ne_mmle):
        ax.axvline(ne_mmle, color="#c0392b", linestyle="-", linewidth=1.5,
                   label=f"Fitted Ne_MMLE = {ne_mmle:.2f}")
    if np.isfinite(ne_lo) and np.isfinite(ne_hi):
        ax.axvspan(ne_lo, ne_hi, color="#c0392b", alpha=0.12,
                   label=f"95% CI [{ne_lo:.2f}, {ne_hi:.2f}]")

    ax.set_xlabel("Candidate Ne")
    ax.set_ylabel(r"$-2\,(\log L - \log L_{\max})$")
    ax.set_title(f"MMLE marginal log-likelihood profile (model={meta.get('model', '?')})")
    # Cap y-axis so the curve and 3.841 line are both visible; values
    # outside the cap usually correspond to far-off-optimum Ne.
    finite_y = df["mmle_delta_2ll"][np.isfinite(df["mmle_delta_2ll"])]
    if len(finite_y) > 0:
        y_cap = float(np.clip(np.nanpercentile(finite_y, 99), 20.0, 200.0))
        ax.set_ylim(-1.0, y_cap)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)


def plot_kimura_panel(ax, df: pd.DataFrame, meta: dict) -> None:
    """Right panel: Kimura per-pair SSR vs Ne (normalised)."""
    ax.plot(df["ne_candidate"], df["kimura_norm_ssr"],
            color="#2c3e50", linewidth=2.0,
            label="Kimura SSR / SSR_min")

    # Analytic and grid-best Kimura Ne under SSR.
    ne_kim_ssr  = _meta_float(meta, "best_ne_kimura_analytic")
    ne_kim_grid = _meta_float(meta, "best_ne_kimura_on_grid")
    ne_kim_mom  = _meta_float(meta, "kimura_ne_moments")
    ne_kim_lo   = _meta_float(meta, "kimura_ne_ci_low")
    ne_kim_hi   = _meta_float(meta, "kimura_ne_ci_high")

    if np.isfinite(ne_kim_ssr):
        ax.axvline(ne_kim_ssr, color="#2c3e50", linestyle="-", linewidth=1.5,
                   label=f"Best Ne_Kimura_SSR = {ne_kim_ssr:.2f}")
    elif np.isfinite(ne_kim_grid):
        ax.axvline(ne_kim_grid, color="#2c3e50", linestyle="-", linewidth=1.5,
                   label=f"Best Ne_Kimura_SSR = {ne_kim_grid:.2f}")

    if np.isfinite(ne_kim_mom):
        ax.axvline(ne_kim_mom, color="#3498db", linestyle="--", linewidth=1.3,
                   label=f"Wonnapinij Ne_Kimura = {ne_kim_mom:.2f}")
    if np.isfinite(ne_kim_lo) and np.isfinite(ne_kim_hi):
        ax.axvspan(ne_kim_lo, ne_kim_hi, color="#3498db", alpha=0.12,
                   label=f"Wonnapinij 95% CI [{ne_kim_lo:.2f}, {ne_kim_hi:.2f}]")

    # deCODE 2024 Cell reference (Ne = 3).
    ax.axvline(3.0, color="#27ae60", linestyle=":", linewidth=1.3,
               label="deCODE reference Ne = 3")

    ax.set_xlabel("Candidate Ne")
    ax.set_ylabel("Kimura SSR / SSR_min")
    ax.set_title("Kimura per-pair SSR profile")
    # Most users care about the bowl near the minimum; cap at 5x for clarity.
    finite_y = df["kimura_norm_ssr"][np.isfinite(df["kimura_norm_ssr"])]
    if len(finite_y) > 0:
        y_cap = float(min(np.nanpercentile(finite_y, 99), 10.0))
        y_cap = max(y_cap, 3.0)
        ax.set_ylim(0.95, y_cap)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)


def make_plot(tsv_path: str, output_path: str,
              dpi: int, figsize: tuple, title: str) -> None:
    meta, df = load_ne_profile(tsv_path)
    if df.empty:
        raise RuntimeError(
            "[plot_ne_profile] No usable rows in profile (all mmle_log_lik are -inf?)")

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=figsize)
    plot_mmle_panel(ax_left, df, meta)
    plot_kimura_panel(ax_right, df, meta)

    if title:
        fig.suptitle(title, fontsize=12)

    n_pairs = meta.get("n_pairs_used", "?")
    model   = meta.get("model", "?")
    ne_mmle_disp = meta.get("fitted_ne_mmle", meta.get("fitted_ne_mle", "?"))
    fig.text(0.5, 0.005,
             f"n_pairs_used = {n_pairs}, model = {model},  "
             f"Ne_MMLE = {ne_mmle_disp},  "
             f"Ne_Kimura_SSR = {meta.get('best_ne_kimura_analytic', '?')},  "
             f"Ne_Kimura_moments = {meta.get('kimura_ne_moments', '?')}",
             ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the MMLE vs Kimura Ne-profile diagnostic emitted by "
                    "`mitoquest ne-estimate --ne-profile`.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
                        help="Ne-profile TSV produced by `mitoquest ne-estimate "
                             "--ne-profile`.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output figure path (.png/.pdf/.svg).")
    parser.add_argument("--dpi", type=int, default=150,
                        help="Output DPI for raster formats.")
    parser.add_argument("--title", default="",
                        help="Optional super-title (e.g. cohort name).")
    parser.add_argument("--figsize", default="13,5.2",
                        help="Figure size 'W,H' in inches.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        w, h = (float(x) for x in args.figsize.split(","))
    except ValueError:
        sys.exit(f"[plot_ne_profile] Invalid --figsize: {args.figsize!r} "
                 "(expected 'W,H')")
    make_plot(args.input, args.output, args.dpi, (w, h), args.title)


if __name__ == "__main__":
    main()
