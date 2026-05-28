#!/usr/bin/env python3
"""
validate_estimators.py
======================

Three-estimator validation harness for `mitoquest ne-estimate`.

This script generates a battery of synthetic mother-child cohorts that
target the *known statistical strengths and weaknesses* of each
available Ne estimator, then runs every estimator on every cohort and
emits a single comparison table.

Estimators under test
---------------------
  * Continuous MLE  (`--model continuous`, default since v1.8.2):
      p_child | p_m  ~  Beta(p_m*(Ne-1), (1-p_m)*(Ne-1))
      c_alt   | p_c  ~  Binomial(c_dp, p_c)
      Real-valued Ne via two-phase golden-section search (since v1.8.3).

  * Discrete MLE    (`--model discrete`):
      k     ~ BetaBinomial(Ne, m_alt+1, m_ref+1)
      c_alt ~ Binomial(c_dp, k/Ne)
      Integer Ne grid; can be upward-biased on very deep mtDNA reads.

  * Kimura cross-check (`--cross-check kimura`):
      Wonnapinij/Kimura method-of-moments,
      V = sum(d_i - s_i) / sum(w_i),  Ne = 1 / V.
      Plug-in sampling-error correction; closed-form, no MLE; subject
      to Jensen's-inequality bias and outlier inflation.

Scenario design
---------------
Each scenario is chosen to highlight one specific behaviour:

  S1  baseline_discrete       --bottleneck-model discrete,    true_ne=5
                              MLE_discrete should win or tie.
  S2  baseline_continuous     --bottleneck-model continuous,  true_ne=5
                              MLE_continuous should win or tie.
  S3  fractional_ne           --bottleneck-model continuous,  true_ne=7.5
                              Only MLE_continuous returns a fractional
                              Ne; MLE_discrete must round to 7 or 8.
  S4  small_cohort            n_pairs=20  (vs default 200)
                              CIs widen; Kimura noise dominates.
  S5  low_depth               m_dp = c_dp = 100 (vs default 2000)
                              Sampling correction is unreliable; Kimura
                              tends upward, both MLEs handle it better.
  S6  outliers_5pct           5% of (site, pair) cells forced to extreme
                              drift -> Kimura inflates strongly while
                              MLEs degrade more gracefully.
  S7  large_ne                true_ne = 30 -> all three should converge.

Usage
-----
  bash run_demo.sh                        # quick smoke test (one scenario)
  python3 validate_estimators.py          # run all 7 scenarios
  python3 validate_estimators.py --scenarios baseline_continuous outliers_5pct
  python3 validate_estimators.py --keep-cohorts  # keep generated VCFs for inspection

Outputs
-------
  scenarios/<name>/cohort.vcf, cohort.fam, cohort.transmission_pairs.tsv
  scenarios/<name>/ne.continuous.json
  scenarios/<name>/ne.discrete.json
  validation_results.md              <-- summary comparison table
  validation_results.tsv             <-- machine-readable copy
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

HERE      = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent.parent
DEFAULT_BIN = REPO_ROOT / "bin" / "mitoquest"


# ---------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------

@dataclass
class Scenario:
    """A single (cohort recipe, expected behaviour) pair."""
    name:               str
    description:        str           # one-line summary printed in the table
    true_ne:            float
    n_pairs:            int           = 200
    n_sites:            int           = 20
    m_dp:               int           = 2000
    c_dp:               int           = 2000
    vaf_low:            float         = 0.15
    vaf_high:           float         = 0.85
    bottleneck_model:   str           = "continuous"
    n_multiallelic:     int           = 0
    missing_gt_rate:    float         = 0.0
    outlier_frac:       float         = 0.0
    n_generations:      int           = 1
    seed:               int           = 20260528
    expected_winner:    str           = "MLE_continuous"  # purely informational

    def synthesize_args(self) -> List[str]:
        return [
            "--n-pairs",         str(self.n_pairs),
            "--n-sites",         str(self.n_sites),
            "--true-ne",         f"{self.true_ne}",
            "--m-dp",            str(self.m_dp),
            "--c-dp",            str(self.c_dp),
            "--vaf-low",         f"{self.vaf_low}",
            "--vaf-high",        f"{self.vaf_high}",
            "--bottleneck-model", self.bottleneck_model,
            "--n-multiallelic",  str(self.n_multiallelic),
            "--missing-gt-rate", f"{self.missing_gt_rate}",
            "--outlier-frac",    f"{self.outlier_frac}",
            "--n-generations",   str(self.n_generations),
            "--seed",            str(self.seed),
        ]

    @property
    def expected_apparent_ne(self) -> float:
        """Wright-Fisher prediction of the *apparent* per-generation Ne
        seen by a single-generation estimator: 1 / (1 - (1-1/Ne)^g).

        For g=1 this reduces to the true Ne.  For g>=2 it is strictly
        smaller than the true Ne (drift accumulates).
        """
        F_g = 1.0 - (1.0 - 1.0 / self.true_ne) ** self.n_generations
        return 1.0 / F_g if F_g > 0 else float("inf")


SCENARIOS: List[Scenario] = [
    Scenario(name="S1_baseline_discrete",
             description="discrete bottleneck, integer true_ne=5 (matches discrete MLE)",
             true_ne=5,
             bottleneck_model="discrete",
             expected_winner="MLE_discrete"),

    Scenario(name="S2_baseline_continuous",
             description="continuous Beta-diffusion, true_ne=5 (matches continuous MLE)",
             true_ne=5,
             bottleneck_model="continuous",
             expected_winner="MLE_continuous"),

    Scenario(name="S3_fractional_ne",
             description="fractional true_ne=7.5 (only continuous MLE can resolve it)",
             true_ne=7.5,
             bottleneck_model="continuous",
             expected_winner="MLE_continuous"),

    Scenario(name="S4_small_cohort",
             description="n_pairs=20 (small sample, wide CIs everywhere)",
             true_ne=5,
             n_pairs=20,
             bottleneck_model="continuous",
             expected_winner="MLE_continuous"),

    Scenario(name="S5_low_depth",
             description="depth=100 (sampling correction unreliable, Kimura biased up)",
             true_ne=5,
             m_dp=100, c_dp=100,
             bottleneck_model="continuous",
             expected_winner="MLE_continuous"),

    Scenario(name="S6_outliers_5pct",
             description="5% high-drift outliers (NUMTs/errors) inflate Kimura",
             true_ne=5,
             outlier_frac=0.05,
             bottleneck_model="continuous",
             expected_winner="MLE_continuous"),

    Scenario(name="S7_large_ne",
             description="true_ne=30 (all three estimators should converge)",
             true_ne=30,
             bottleneck_model="continuous",
             expected_winner="all converge"),

    # ------------------------------------------------------------------
    # Multi-generation pedigree scenarios.
    # All single-generation estimators implicitly assume mother-child;
    # feeding them g-generation pairs makes them estimate the *apparent*
    # Ne_app = 1 / (1 - (1 - 1/Ne)^g), which is biased downward.
    # These scenarios quantify that bias.
    # ------------------------------------------------------------------
    Scenario(name="S8_g2_ne10",
             description="g=2, true_ne=10 (apparent Ne ~5.26 due to 2-gen drift)",
             true_ne=10,
             n_generations=2,
             bottleneck_model="continuous",
             expected_winner="MLE_cont matches 1/F_g"),

    Scenario(name="S9_g5_ne10",
             description="g=5, true_ne=10 (apparent Ne ~2.44, large pedigree bias)",
             true_ne=10,
             n_generations=5,
             bottleneck_model="continuous",
             expected_winner="MLE_cont matches 1/F_g"),

    Scenario(name="S10_g10_ne30",
             description="g=10, true_ne=30 (apparent Ne ~3.48, severe drift accumulation)",
             true_ne=30,
             n_generations=10,
             bottleneck_model="continuous",
             expected_winner="MLE_cont matches 1/F_g"),

    Scenario(name="S11_g3_ne10_discrete",
             description="g=3, true_ne=10, discrete sampler (apparent Ne ~3.7)",
             true_ne=10,
             n_generations=3,
             bottleneck_model="discrete",
             expected_winner="all match 1/F_g"),
]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

@dataclass
class EstimatorResult:
    ne:                   Optional[float] = None
    ci_low:               Optional[float] = None
    ci_high:              Optional[float] = None
    pairs_used:           Optional[int]   = None
    max_loglik:           Optional[float] = None
    kimura_ne:            Optional[float] = None
    kimura_ci_low:        Optional[float] = None
    kimura_ci_high:       Optional[float] = None
    kimura_b:             Optional[float] = None
    kimura_n_informative: Optional[int]   = None


def _strip_comments(text: str) -> str:
    """ne-estimate JSON output can be prefixed by a `#mitoquest_...` line.
    Strip any leading lines that start with `#` before json.loads."""
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("#")
    )


def _read_json(path: Path) -> dict:
    return json.loads(_strip_comments(path.read_text()))


def parse_ne_json(json_path: Path) -> EstimatorResult:
    j = _read_json(json_path)
    r = EstimatorResult(
        ne          = j.get("Ne"),
        ci_low      = j.get("CI_95_Low"),
        ci_high     = j.get("CI_95_High"),
        pairs_used  = j.get("Pairs_Used"),
        max_loglik  = j.get("Max_LogLik"),
    )
    k = j.get("Kimura_Cross_Check") or {}
    r.kimura_ne            = k.get("Ne_Kimura")
    r.kimura_ci_low        = k.get("Ne_Kimura_CI_95_Low")
    r.kimura_ci_high       = k.get("Ne_Kimura_CI_95_High")
    r.kimura_b             = k.get("b")
    r.kimura_n_informative = k.get("N_Informative")
    return r


def run(cmd: List[str], log_to: Optional[Path] = None) -> None:
    """Run `cmd` and abort on non-zero exit."""
    if log_to is not None:
        with log_to.open("w") as fh:
            res = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    else:
        res = subprocess.run(cmd)
    if res.returncode != 0:
        sys.stderr.write(f"\n[validate] Command failed (rc={res.returncode}):\n")
        sys.stderr.write("           " + " ".join(cmd) + "\n")
        if log_to is not None:
            sys.stderr.write(f"           Log: {log_to}\n")
        sys.exit(res.returncode)


def synthesize_cohort(synth_py: Path, scenario: Scenario,
                      out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["python3", str(synth_py),
           "--out-dir", str(out_dir),
           *scenario.synthesize_args()]
    run(cmd, log_to=out_dir / "synthesize.log")


def run_estimator(mitoquest: Path, tsv: Path, model: str,
                  out_json: Path, threads: int = 4) -> EstimatorResult:
    cmd = [
        str(mitoquest), "ne-estimate",
        "-i", str(tsv),
        "--model", model,
        "--min-vaf", "0.10", "--max-vaf", "0.90",
        "--min-ne", "1", "--max-ne", "200",
        "--cross-check", "kimura",
        "--kimura-bootstrap", "500",
        "-t", str(threads),
        "-o", str(out_json),
    ]
    run(cmd, log_to=out_json.with_suffix(".log"))
    return parse_ne_json(out_json)


# ---------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------

def fmt_ci(low, high) -> str:
    if low is None or high is None:
        return "n/a"
    return f"{low:.2f}-{high:.2f}"


def fmt_ne(x, decimals: int = 2) -> str:
    if x is None:
        return "n/a"
    return f"{x:.{decimals}f}"


def collect_rows(results: Dict[str, Dict[str, EstimatorResult]],
                 scenarios: List[Scenario]) -> List[List[str]]:
    """Build the long-form row list for both markdown and TSV output."""
    rows: List[List[str]] = []
    for sc in scenarios:
        per = results.get(sc.name) or {}
        cont = per.get("continuous") or EstimatorResult()
        disc = per.get("discrete")   or EstimatorResult()
        # Kimura values come out of either run; prefer the continuous run
        # because its informative-pair set is identical (same VAF gates).
        kim  = cont if cont.kimura_ne is not None else disc

        rows.append([
            sc.name,
            fmt_ne(sc.true_ne, 2),
            str(sc.n_generations),
            fmt_ne(sc.expected_apparent_ne, 2),
            str(sc.n_pairs),
            f"{sc.m_dp}/{sc.c_dp}",
            sc.bottleneck_model,
            f"{sc.outlier_frac:.2f}",
            fmt_ne(cont.ne),     fmt_ci(cont.ci_low,  cont.ci_high),
            fmt_ne(disc.ne, 0),  fmt_ci(disc.ci_low,  disc.ci_high),
            fmt_ne(kim.kimura_ne),
            fmt_ci(kim.kimura_ci_low, kim.kimura_ci_high),
            str(cont.pairs_used or disc.pairs_used or 0),
            sc.expected_winner,
        ])
    return rows


HEADER = [
    "Scenario", "True_Ne", "g", "Ne_app_thy",
    "N_pairs", "m_dp/c_dp",
    "Sim_model", "Outliers",
    "MLE_cont", "CI_cont",
    "MLE_disc", "CI_disc",
    "Ne_Kimura", "CI_Kimura",
    "Pairs_Used", "Expected",
]


def write_markdown(rows: List[List[str]], out_path: Path,
                   binary_path: Path) -> None:
    sep_line = "| " + " | ".join([":---:"] * len(HEADER)) + " |"
    lines = [
        "# Three-estimator validation results",
        "",
        f"Binary: `{binary_path}`",
        "",
        "Each row reports the same `cohort.transmission_pairs.tsv` fed to",
        "`mitoquest ne-estimate` three times: once with `--model continuous`",
        "(MLE_cont), once with `--model discrete` (MLE_disc), with",
        "`--cross-check kimura` enabled in both runs (Ne_Kimura).",
        "",
        "Column `g` is the simulated pedigree depth (number of generations",
        "between mother and child). Column `Ne_app_thy` is the Wright-Fisher",
        "prediction of the *apparent* per-generation Ne that any single-",
        "generation estimator should report: `1 / (1 - (1 - 1/Ne)^g)`.",
        "For `g=1` this equals the true Ne; for `g>1` it is strictly smaller.",
        "",
        "| " + " | ".join(HEADER) + " |",
        sep_line,
    ]
    for r in rows:
        lines.append("| " + " | ".join(r) + " |")
    lines.append("")
    lines.append("Per-scenario raw JSON files: `scenarios/<name>/ne.continuous.json` "
                 "and `scenarios/<name>/ne.discrete.json`.")
    lines.append("")
    out_path.write_text("\n".join(lines) + "\n")


def write_tsv(rows: List[List[str]], out_path: Path) -> None:
    out_path.write_text("\t".join(HEADER) + "\n"
                        + "\n".join("\t".join(r) for r in rows) + "\n")


def print_table(rows: List[List[str]]) -> None:
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(HEADER)]
    sep = "  ".join("-" * w for w in widths)
    line = lambda cells: "  ".join(c.ljust(widths[i]) for i, c in enumerate(cells))
    print(line(HEADER))
    print(sep)
    for r in rows:
        print(line(r))


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mitoquest", default=str(DEFAULT_BIN),
                   help=f"Path to the `mitoquest` binary [{DEFAULT_BIN}].")
    p.add_argument("--scenarios", nargs="+", default=None,
                   help="Run only the named scenarios (default: all).")
    p.add_argument("--out-dir",   default=str(HERE / "scenarios"),
                   help="Where to write per-scenario cohorts and JSON output.")
    p.add_argument("--report-md", default=str(HERE / "validation_results.md"))
    p.add_argument("--report-tsv", default=str(HERE / "validation_results.tsv"))
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-cohorts", action="store_true",
                   help="Keep the generated cohort.{vcf,fam,tsv} files.")
    p.add_argument("--clean", action="store_true",
                   help="Remove --out-dir before running.")
    args = p.parse_args()

    mitoquest = Path(args.mitoquest)
    if not mitoquest.exists() or not os.access(mitoquest, os.X_OK):
        sys.stderr.write(f"[validate] mitoquest binary not found / not executable: {mitoquest}\n")
        sys.stderr.write(f"           build it first with: cmake --build build -j8\n")
        sys.exit(1)

    synth_py = HERE / "synthesize.py"
    out_dir  = Path(args.out_dir)
    if args.clean and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    selected: List[Scenario]
    if args.scenarios:
        names = set(args.scenarios)
        selected = [s for s in SCENARIOS if s.name in names]
        unknown = names - {s.name for s in SCENARIOS}
        if unknown:
            sys.stderr.write(f"[validate] Unknown scenarios: {', '.join(sorted(unknown))}\n")
            sys.exit(1)
    else:
        selected = SCENARIOS

    results: Dict[str, Dict[str, EstimatorResult]] = {}

    for sc in selected:
        print(f"\n=== Scenario {sc.name} ===")
        print(f"    {sc.description}")
        print(f"    true_ne={sc.true_ne}  g={sc.n_generations}  "
              f"Ne_app(theory)={sc.expected_apparent_ne:.2f}  "
              f"n_pairs={sc.n_pairs}  depth={sc.m_dp}/{sc.c_dp}  "
              f"sim={sc.bottleneck_model}  outliers={sc.outlier_frac}")

        sc_dir = out_dir / sc.name
        synthesize_cohort(synth_py, sc, sc_dir)
        tsv = sc_dir / "cohort.transmission_pairs.tsv"

        results[sc.name] = {}
        for model in ("continuous", "discrete"):
            json_out = sc_dir / f"ne.{model}.json"
            r = run_estimator(mitoquest, tsv, model, json_out,
                              threads=args.threads)
            results[sc.name][model] = r
            print(f"    --model {model:<10}  Ne={fmt_ne(r.ne)}  "
                  f"CI={fmt_ci(r.ci_low, r.ci_high)}  "
                  f"pairs={r.pairs_used}  "
                  f"Ne_Kimura={fmt_ne(r.kimura_ne)}")

        if not args.keep_cohorts:
            # Drop the bulky VCF / TSV but keep the JSON output for review.
            for f in ("cohort.vcf", "cohort.fam",
                      "cohort.transmission_pairs.tsv",
                      "synthesize.log"):
                fp = sc_dir / f
                if fp.exists(): fp.unlink()

    rows = collect_rows(results, selected)
    print()
    print_table(rows)

    write_markdown(rows, Path(args.report_md), mitoquest)
    write_tsv(rows, Path(args.report_tsv))

    print(f"\n[validate] Wrote {args.report_md}")
    print(f"[validate] Wrote {args.report_tsv}")


if __name__ == "__main__":
    main()
