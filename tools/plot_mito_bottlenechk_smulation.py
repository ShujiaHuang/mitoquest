#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MitoQuest multi-generation bottleneck plot
==========================================
Reads the JSON output of ``mitoquest ne-estimate`` (see ``src/ne_estimate.cpp``)
and renders the multi-generation bottleneck retention parameter ``b`` across
pair types of varying genealogical distance.  The Continuous Beta-diffusion
MMLE and the optional Kimura / Wonnapinij cross-check are overlaid with their
95 % confidence-interval bands.

Theory (consistent with ``ne_estimate.cpp``)
--------------------------------------------
``ne_estimate.cpp`` models a single mother-child transmission as a haploid
Wright-Fisher / Beta-diffusion step::

    p_child | p_mother  ~  Beta( p_m * (Ne - 1),  (1 - p_m) * (Ne - 1) )
    Var(p_child | p_mother)  =  p_m * (1 - p_m) / Ne

The Wonnapinij retention parameter is therefore::

    1 - b  =  Var(p_child - p_mother) / [p_m * (1 - p_m)]  =  1 / Ne
       b   =  1 - 1 / Ne                                  (single transmission)

Extending neutrally to ``t`` transmissions between two samples (label
``"m-n"`` with ``t = m + n``) gives the haploid Wright-Fisher accumulation::

    b_t  =  ( 1 - 1 / Ne ) ** t

For mother-child (``t = 1``) this reproduces the value reported in the JSON
field ``Kimura_Cross_Check.b`` exactly, so the plotted curve and the cpp
output are on the same scale.

Both the Continuous Beta-diffusion model and the Kimura / Wonnapinij model
share this first-moment retention prediction; they differ only in the
fitted ``Ne``.

Dependencies: numpy, seaborn, matplotlib.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns


# ---------------------------------------------------------------------------
# Pair-type catalogue
# ---------------------------------------------------------------------------
# Label "m-n" -> (m, n) where m and n are the number of transmissions on the
# two sides of the common ancestor.  Total drift distance is t = m + n.
PAIR_TRANSMISSIONS: Dict[str, Tuple[int, int]] = {
    "0-1": (0, 1),  # mother - child           (t = 1)
    "0-2": (0, 2),  # grandmother - grandchild (t = 2)
    "1-1": (1, 1),  # full siblings            (t = 2)
    "2-2": (2, 2),  # first cousins            (t = 4)
    "3-3": (3, 3),  # second cousins           (t = 6)
    "4-4": (4, 4),
    "5-5": (5, 5),
    "6-6": (6, 6),
    "7-7": (7, 7),
}


# ---------------------------------------------------------------------------
# Core model
# ---------------------------------------------------------------------------

def retention_b(ne: float, t: int) -> float:
    """Haploid Wright-Fisher retention parameter after ``t`` transmissions.

    Matches ``ne_estimate.cpp`` continuous Beta-diffusion model at ``t = 1``:
        b = 1 - 1 / Ne
    Multi-generation extension under neutral drift:
        b_t = (1 - 1 / Ne) ** t
    """
    if t <= 0:
        return 1.0
    if ne <= 1.0:
        # Ne = 1 means complete fixation in one transmission, so b = 0.
        return 0.0
    return float(np.clip((1.0 - 1.0 / ne) ** t, 0.0, 1.0))


def compute_curve(ne: float, pair_types: List[str]) -> List[float]:
    return [retention_b(ne, sum(PAIR_TRANSMISSIONS[pt])) for pt in pair_types]


# ---------------------------------------------------------------------------
# JSON parsing
# ---------------------------------------------------------------------------

def _safe_num(v) -> Optional[float]:
    """Convert a JSON value to a finite float.

    ``ne_estimate.cpp`` emits non-finite numbers as JSON strings such as
    ``"Infinity"`` or ``"NaN"``; this helper returns ``None`` for those
    so the plotting code can hide the affected curve / band.
    """
    if v is None:
        return None
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v) if np.isfinite(v) else None
    if isinstance(v, str):
        try:
            f = float(v)
        except ValueError:
            return None
        return f if np.isfinite(f) else None
    return None


def load_ne_json(filepath: str) -> dict:
    """Parse ne-estimate JSON output, skipping leading '#' provenance lines."""
    if not os.path.exists(filepath):
        print(f"[ERROR] File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    with open(filepath, "r") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]

    try:
        data = json.loads("".join(lines))
    except json.JSONDecodeError as e:
        print(f"[ERROR] JSON parse failed: {e}", file=sys.stderr)
        sys.exit(1)

    for key in ("Ne", "CI_95_Low", "CI_95_High"):
        if key not in data:
            print(f"[ERROR] Missing required JSON field: '{key}'", file=sys.stderr)
            sys.exit(1)

    print(
        f"[INFO] Continuous MMLE: Ne = {data['Ne']:.4f}  "
        f"(95% CI [{data['CI_95_Low']:.4f}, {data['CI_95_High']:.4f}])"
    )

    if "Kimura_Cross_Check" in data:
        k = data["Kimura_Cross_Check"]
        ne_k = _safe_num(k.get("Ne_Kimura"))
        ne_lo = _safe_num(k.get("Ne_Kimura_CI_95_Low"))
        ne_hi = _safe_num(k.get("Ne_Kimura_CI_95_High"))
        if ne_k is None:
            print("[INFO] Kimura cross-check: Ne_Kimura is non-finite; skipping band.")
        else:
            ci_str = (
                f" (95% CI [{ne_lo:.4f}, {ne_hi:.4f}])"
                if ne_lo is not None and ne_hi is not None
                else ""
            )
            print(f"[INFO] Kimura cross-check: Ne_Kimura = {ne_k:.4f}{ci_str}")
    else:
        print("[INFO] Kimura cross-check not present in JSON "
              "(re-run with `--cross-check kimura` to add it).")

    return data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_bottleneck(json_data: dict, output_path: str) -> None:
    """Render the Continuous MMLE / Kimura cross-check overlay."""
    sns.set_theme(style="whitegrid", context="talk", font_scale=0.9)

    pair_types = list(PAIR_TRANSMISSIONS.keys())
    x_pos = np.arange(len(pair_types))

    fig, ax = plt.subplots(figsize=(13, 7.5))

    # --- Reference grid: theoretical b vs pair-type curves at canonical Ne. -
    grid_nes = [1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100]
    grey_palette = sns.color_palette("Greys", n_colors=len(grid_nes) + 2)
    for i, ne in enumerate(grid_nes):
        y_vals = compute_curve(ne, pair_types)
        ax.plot(
            x_pos, y_vals,
            color=grey_palette[i + 1], alpha=0.55, linewidth=0.8, zorder=1,
        )
        ax.text(
            len(pair_types) - 0.05, y_vals[-1], f"$N_e$={ne:g}",
            fontsize=7, color="gray", va="center", alpha=0.75,
        )

    # --- Continuous MMLE: matches the JSON "Ne" field. ----------------------
    main_ne = float(json_data["Ne"])
    ci_lo = float(json_data["CI_95_Low"])
    ci_hi = float(json_data["CI_95_High"])
    main_y = compute_curve(main_ne, pair_types)
    # Larger Ne -> larger retention b, so map ci_high -> upper edge of band.
    band_lower = compute_curve(ci_lo, pair_types)
    band_upper = compute_curve(ci_hi, pair_types)

    color_main = "#B22222"
    ax.fill_between(x_pos, band_lower, band_upper,
                    color=color_main, alpha=0.15, zorder=5)
    ax.plot(
        x_pos, main_y,
        color=color_main, linewidth=2.8, marker="o", markersize=6, zorder=10,
        label=(f"Continuous MMLE  ($N_e$={main_ne:.2f}, "
               f"95% CI [{ci_lo:.2f}, {ci_hi:.2f}])"),
    )

    # --- Kimura cross-check (optional). -------------------------------------
    kimura = json_data.get("Kimura_Cross_Check")
    if kimura is not None:
        k_ne = _safe_num(kimura.get("Ne_Kimura"))
        k_lo = _safe_num(kimura.get("Ne_Kimura_CI_95_Low"))
        k_hi = _safe_num(kimura.get("Ne_Kimura_CI_95_High"))

        if k_ne is not None:
            color_kimura = "#1E6091"
            k_y = compute_curve(k_ne, pair_types)

            label = f"Kimura cross-check  ($N_e$={k_ne:.2f}"
            if k_lo is not None and k_hi is not None:
                label += f", 95% CI [{k_lo:.2f}, {k_hi:.2f}])"
                k_band_lo = compute_curve(k_lo, pair_types)
                k_band_hi = compute_curve(k_hi, pair_types)
                ax.fill_between(x_pos, k_band_lo, k_band_hi,
                                color=color_kimura, alpha=0.12, zorder=4)
            else:
                label += ")"

            ax.plot(
                x_pos, k_y,
                color=color_kimura, linewidth=2.8, linestyle="--",
                marker="s", markersize=6, zorder=9, label=label,
            )

    # --- Cosmetic ----------------------------------------------------------
    ax.set_xticks(x_pos)
    ax.set_xticklabels(pair_types, fontsize=12)
    ax.set_xlabel("Type of pair (transmission distance $m$-$n$)", fontsize=14)
    ax.set_ylabel(r"Bottleneck retention $b = (1 - 1/N_e)^{\,t}$", fontsize=14)

    model_name = json_data.get("Model", "continuous")
    pairs_used = json_data.get("Pairs_Used", "?")
    vaf_min = json_data.get("Min_VAF", "?")
    vaf_max = json_data.get("Max_VAF", "?")
    ax.set_title(
        "Mitochondrial bottleneck: Continuous MMLE vs Kimura cross-check\n"
        f"{pairs_used} pairs | maternal VAF in [{vaf_min}, {vaf_max}] | "
        f"model = {model_name}",
        fontsize=14, fontweight="bold",
    )
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
    ax.set_xlim(-0.3, len(pair_types) - 0.3)
    ax.set_ylim(0.0, 1.02)
    ax.legend(fontsize=11, loc="lower left", frameon=True, edgecolor="gray")
    sns.despine(ax=ax)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"\n[SUCCESS] Figure saved to: {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot multi-generation mtDNA bottleneck retention from the JSON "
            "produced by `mitoquest ne-estimate`."
        ),
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to the ne-estimate JSON output file.",
    )
    parser.add_argument(
        "--output", "-o", default="mito_bottleneck.png",
        help="Output figure path (PNG / PDF / SVG).  Default: %(default)s",
    )
    args = parser.parse_args()

    json_data = load_ne_json(args.input)
    plot_bottleneck(json_data, args.output)


if __name__ == "__main__":
    main()
