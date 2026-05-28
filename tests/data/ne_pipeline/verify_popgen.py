#!/usr/bin/env python3
"""
verify_popgen.py — Verify that `synthesize.py` satisfies Wright-Fisher
one-generation AND multi-generation moments before any read-level
noise is added.

This complements `validate_estimators.py`:

  * `validate_estimators.py` tests the *estimators* on simulated cohorts.
  * `verify_popgen.py`        tests the *simulator itself* against theory,
                              upstream of read sampling and any estimator.

The check
---------

Wright-Fisher predicts that, for a single generation through an
effective population of size Ne, the post-bottleneck heteroplasmy
`p_c` given pre-bottleneck `p_m` satisfies:

    E[p_c | p_m]       =  p_m                       (martingale)
    Var(p_c | p_m)     =  p_m * (1 - p_m) / Ne      (drift variance)

For `g` independent generations of drift starting from `p_m = p_0`,
iterating the Wright-Fisher transition gives the standard result:

    E[p_g | p_0]       =  p_0                       (still martingale)
    Var(p_g | p_0)     =  p_0 * (1 - p_0) * F_g
        where  F_g     =  1 - (1 - 1/Ne)^g.

The quantity `F_g` is the cumulative drift coefficient.  For g=1 it
recovers `1/Ne`; for g -> infinity it tends to 1 (full fixation).
Defining the dimensionless drift quantity

    F_i := (p_g_i - p_0_i)^2 / (p_0_i * (1 - p_0_i))

we must therefore see

    E[F_i]  =  F_g  =  1 - (1 - 1/Ne)^g

This is exactly what the Wonnapinij/Kimura `b` parameter measures:
`b = 1 - F_g`, and the *apparent* per-generation Ne returned by an
estimator that assumes a single generation is

    Ne_apparent  =  1 / F_g  =  1 / (1 - (1 - 1/Ne)^g)

so multi-generation pedigrees deflate the apparent Ne (a known
limitation of single-generation estimators).

Both bottleneck samplers in `synthesize.py` are designed to satisfy
one-generation Wright-Fisher exactly in expectation, and iterating
them `g` times is the standard way to simulate `g` generations of
independent drift:

  * discrete:   k ~ Multinomial(round(Ne), p_in); p_out = k / Ne.
                Bernoulli sum, Var(p_out|p_in) = p_in*(1-p_in)/Ne.
  * continuous: p_out ~ Dirichlet(p_in * (Ne - 1)).
                Kimura diffusion limit; for biallelic this reduces to
                Beta(p_alt*(Ne-1), (1-p_alt)*(Ne-1)).
                Direct calculation:
                  alpha + beta            = Ne - 1
                  Var(Beta(alpha, beta))  = alpha*beta / ((alpha+beta)^2 * (alpha+beta+1))
                                          = p_in(1-p_in)(Ne-1)^2 / ((Ne-1)^2 * Ne)
                                          = p_in(1-p_in) / Ne.

If the simulator is correct, the empirical mean of `F_i` across many
(p_0, p_g) draws must converge to `F_g`, with relative error scaling
as `1/sqrt(N_samples)`.  This script checks that for a panel of
(model, Ne, g) triples and prints a comparison table.

Usage:
    python3 verify_popgen.py                  # default panel
    python3 verify_popgen.py --n-samples 500000

Author: Shujia Huang (hshujia@qq.com)
Date:   2026-05-28
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import random
import sys
from pathlib import Path
from typing import List, Tuple

HERE = Path(__file__).resolve().parent

# Import synthesize.py as a module so we can reuse its gamma_sample()
# without copy-pasting it (keeps this check honest).
_spec = importlib.util.spec_from_file_location("syn", HERE / "synthesize.py")
_syn  = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_syn)


def sample_pc_discrete(p_m: float, true_ne: int, n_generations: int,
                       rng: random.Random) -> float:
    """Iterated one-generation Wright-Fisher: at each step,
    k ~ Binomial(Ne, p_curr), p_curr = k/Ne.  Returns p_g."""
    p = p_m
    for _ in range(n_generations):
        cnt_alt = 0
        for _i in range(true_ne):
            if rng.random() < p:
                cnt_alt += 1
        p = cnt_alt / true_ne
    return p


def sample_pc_continuous(p_m: float, true_ne: float, n_generations: int,
                         rng: random.Random) -> float:
    """Iterated Kimura diffusion: at each step,
    p_curr ~ Beta(p_curr*(Ne-1), (1-p_curr)*(Ne-1)).  Returns p_g."""
    ne1 = float(true_ne) - 1.0
    p = p_m
    for _ in range(n_generations):
        xa = _syn.gamma_sample(max(p       * ne1, 1e-12), rng)
        xb = _syn.gamma_sample(max((1 - p) * ne1, 1e-12), rng)
        s = xa + xb
        p = (xa / s) if s > 0 else 0.0
    return p


def measure(model: str, true_ne: float, n_generations: int,
            n_samples: int,
            vaf_low: float, vaf_high: float,
            seed: int) -> Tuple[float, float, float]:
    """Returns (mean_drift, mean_F, se_F).

    F_i := (p_g - p_0)^2 / (p_0 * (1 - p_0))
    """
    rng = random.Random(seed)
    sum_drift = 0.0
    sum_F     = 0.0
    sum_F2    = 0.0
    n = 0
    ne_int = max(1, int(round(true_ne)))
    for _ in range(n_samples):
        p_0 = rng.uniform(vaf_low, vaf_high)
        if model == "discrete":
            p_g = sample_pc_discrete(p_0, ne_int, n_generations, rng)
        else:
            p_g = sample_pc_continuous(p_0, true_ne, n_generations, rng)
        d = p_g - p_0
        F = (d * d) / (p_0 * (1.0 - p_0))
        sum_drift += d
        sum_F     += F
        sum_F2    += F * F
        n += 1
    mean_drift = sum_drift / n
    mean_F     = sum_F / n
    var_F      = max(0.0, sum_F2 / n - mean_F * mean_F)
    se_F       = math.sqrt(var_F / n)
    return mean_drift, mean_F, se_F


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-samples", type=int, default=200_000,
                   help="Number of (p_0, p_g) draws per row [200000].")
    p.add_argument("--vaf-low",   type=float, default=0.15)
    p.add_argument("--vaf-high",  type=float, default=0.85)
    p.add_argument("--seed",      type=int,   default=20260528)
    args = p.parse_args()

    # Panel: (model, true_ne, n_generations).
    panel: List[Tuple[str, float, int]] = [
        # ---- one-generation baseline (g = 1, F_g = 1/Ne) ---------------
        ("discrete",   3,    1),
        ("discrete",   5,    1),
        ("discrete",   10,   1),
        ("continuous", 3,    1),
        ("continuous", 5,    1),
        ("continuous", 7.5,  1),
        ("continuous", 30,   1),
        # ---- multi-generation drift accumulation -----------------------
        # F_g = 1 - (1 - 1/Ne)^g; the apparent Ne under a single-
        # generation estimator is 1 / F_g.  These rows verify the
        # iterated bottleneck matches theory across pedigree depths.
        ("continuous", 10,   2),    # F_g = 0.190, Ne_app ~ 5.26
        ("continuous", 10,   5),    # F_g = 0.410, Ne_app ~ 2.44
        ("continuous", 10,  10),    # F_g = 0.651, Ne_app ~ 1.54
        ("discrete",   10,   2),    # same F_g, discrete sampler
        ("discrete",   10,   5),
        ("continuous", 30,   5),    # F_g = 0.156, Ne_app ~ 6.39
        ("continuous", 30,  10),    # F_g = 0.288, Ne_app ~ 3.47
    ]

    print(f"# Wright-Fisher one- and multi-generation moment check")
    print(f"# Theory: E[F] = F_g = 1 - (1 - 1/Ne)^g")
    print(f"# n_samples = {args.n_samples}, "
          f"p_0 ~ Uniform[{args.vaf_low}, {args.vaf_high}], "
          f"seed = {args.seed}")
    print()
    header = ("Model",
              "True_Ne",
              "g",
              "E[F]_emp", "E[F]_se",
              "F_g_theory",
              "Ne_app",
              "F_relerr",
              "Mean_drift")
    fmt = "{:<11} {:>8} {:>3} {:>12} {:>10} {:>11} {:>9} {:>10} {:>12}"
    print(fmt.format(*header))
    print("-" * 112)

    all_ok = True
    for model, ne, g in panel:
        mean_drift, mean_F, se_F = measure(
            model, ne, g, args.n_samples,
            args.vaf_low, args.vaf_high, args.seed)
        F_g_theory = 1.0 - (1.0 - 1.0 / ne) ** g
        ne_app     = 1.0 / mean_F if mean_F > 0 else float("inf")
        F_relerr   = (mean_F - F_g_theory) / F_g_theory

        print(fmt.format(
            model,
            f"{ne:.2f}",
            f"{g}",
            f"{mean_F:.6f}",
            f"±{se_F:.5f}",
            f"{F_g_theory:.6f}",
            f"{ne_app:.3f}",
            f"{F_relerr:+.3%}",
            f"{mean_drift:+.5f}",
        ))

        # Sanity: the empirical mean of F must be within 4 SE of theory.
        z = abs(mean_F - F_g_theory) / max(se_F, 1e-12)
        if z > 4.0:
            print(f"  !! FAIL: empirical E[F] is {z:.2f} SE away from F_g theory")
            all_ok = False

    print()
    if all_ok:
        print("[OK] All rows within 4 SE of Wright-Fisher prediction.")
        print("     Both bottleneck samplers correctly reproduce one- and")
        print("     multi-generation drift accumulation F_g = 1 - (1 - 1/Ne)^g.")
        sys.exit(0)
    else:
        print("[FAIL] At least one row deviates from theory; investigate above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
