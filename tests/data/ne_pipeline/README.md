# `tests/data/ne_pipeline/` — Ne-estimate validation suite

This directory ships:

1. A small, deterministic mother-child mtDNA cohort that lets you
   verify `mitoquest trans-prep` and `mitoquest ne-estimate` run
   smoothly end-to-end (the **smoke test**, see `run_demo.sh`).
2. A reproducible **three-estimator validation harness** that
   simulates a battery of one- and multi-generation scenarios designed
   to expose the statistical strengths and weaknesses of each Ne
   estimator (`validate_estimators.py`).
3. A **simulator self-test** (`verify_popgen.py`) that empirically
   verifies the synthetic cohorts satisfy Wright-Fisher one- and
   multi-generation drift moments to within sampling noise — i.e.
   that the data fed to the estimator validation actually obeys the
   population-genetics model the estimators are designed for.

---

## 1. Files

| File                              | Purpose                                                                                           |
| --------------------------------- | ------------------------------------------------------------------------------------------------- |
| `cohort.vcf`                      | Default smoke-test multi-sample VCF (40 samples × 12 mtDNA SNV sites).                            |
| `cohort.fam`                      | PLINK FAM describing the 20 mother-child trios (fathers are placeholders).                        |
| `cohort.transmission_pairs.tsv`   | Pre-baked TSV in the schema emitted by `mitoquest trans-prep`.                                    |
| `synthesize.py`                   | Deterministic simulator (no numpy). Supports **discrete** and **continuous** bottleneck models, optional outlier injection, and **multi-generation** transmission (`--n-generations`). |
| `verify_popgen.py`                | Empirical Wright-Fisher check: draws (p₀, p_g) pairs from each sampler and compares E[F] to the theoretical drift coefficient F_g = 1 − (1 − 1/Ne)ᵍ. |
| `run_demo.sh`                     | One-shot end-to-end smoke test: `trans-prep` → `ne-estimate` (with Kimura cross-check).           |
| `validate_estimators.py`          | Three-estimator validation harness (11 scenarios incl. multi-gen). Generates per-scenario cohorts under `scenarios/<name>/`, runs `ne-estimate` once with `--model continuous` and once with `--model discrete` (both with `--cross-check kimura`), and writes a comparison table. |
| `validation_results.md`           | Markdown comparison table emitted by `validate_estimators.py`.                                    |
| `validation_results.tsv`          | Machine-readable copy of the same table.                                                          |
| `scenarios/<name>/ne.*.json`      | Raw per-estimator JSON output for every scenario.                                                 |
| `ne_diagnostic.py`                | Stand-alone diagnostic that recomputes Kimura on an existing TSV (legacy, kept for parity).       |

---

## 2. The three estimators under test

`mitoquest ne-estimate` exposes two MLE models and one cross-check
estimator. Each is appropriate for a different statistical regime:

### 2.1 Continuous MLE — `--model continuous` (default since v1.8.2)

```txt
p_child | p_mother  ~  Beta( p_m·(Ne−1),  (1−p_m)·(Ne−1) )
c_alt   | p_child   ~  Binomial( c_dp, p_child )
```

Wright-Fisher diffusion limit: post-bottleneck heteroplasmy lives on
the continuous interval `[0, 1]`. Since v1.8.3 the optimizer returns
**real-valued Ne** (two-phase golden-section search, 0.01 precision)
with profile-likelihood 95% CI. Recommended default for mtDNA.

### 2.2 Discrete MLE — `--model discrete`

```txt
k     ~ BetaBinomial( Ne, m_alt+1, m_ref+1 )
c_alt ~ Binomial( c_dp, k / Ne )
```

Hard discrete bottleneck: the child's heteroplasmy is restricted to
the grid `{0, 1/Ne, 2/Ne, …, 1}`. Searches integer Ne; can be
**upward-biased on deep mtDNA reads** when the true biology is
diffusion-like (the grid cannot represent intermediate values, so the
optimizer has to inflate Ne to densify the grid).

### 2.3 Kimura cross-check — `--cross-check kimura`

Wonnapinij/Kimura method-of-moments:

```txt
V  =  Σ (d_i  −  s_i)  /  Σ w_i           (weighted mean drift)
b  =  1  −  V
Ne =  1 / (1 − b) = 1 / V
```

* `d_i = (p_c_i − p_m_i)²`     (observed drift)
* `s_i = p_m·(1−p_m)/m_dp + p_c·(1−p_c)/c_dp`   (sampling correction)
* `w_i = p_m·(1−p_m)`         (information weight)

Closed-form, no optimization. Subject to:

* **Jensen's-inequality bias**: `1/V` is convex, so `E[1/V̂] > 1/E[V̂]`
  → Kimura tends to over-estimate Ne when V̂ is noisy.
* **Plug-in correction noise**: replaces population frequencies with
  noisy point estimates → variance leaks into V.
* **Outlier sensitivity**: a few extreme drifts (NUMTs, mosaicism,
  errors) inflate the mean drift; mitigated by `--kimura-trim 0.10`.

Use Kimura as a **diagnostic cross-check** alongside the MLE, not as
the primary estimate.

---

## 3. Population-genetics foundation of the simulator

Before running the estimators we have to be confident the synthetic
cohorts actually behave like a real Wright-Fisher mtDNA bottleneck.
This section spells out the model and how `synthesize.py` implements it.

### 3.1 Wright-Fisher one-generation drift

For a heteroplasmic site in generation `t`, let `p_t ∈ [0, 1]` be the
frequency of the ALT allele. A single transmission through an
effective population of size `Ne` follows Wright-Fisher:

```txt
k_{t+1} | p_t   ~  Binomial(Ne, p_t)
p_{t+1}         =  k_{t+1} / Ne
```

This gives the two foundational moments:

```txt
E [ p_{t+1} | p_t ]    =  p_t                       (martingale)
Var( p_{t+1} | p_t )   =  p_t (1 − p_t) / Ne        (drift variance)
```

The variance scales as **1/Ne**: tighter bottlenecks (small Ne) drive
larger jumps in heteroplasmy.

For multi-allelic sites (REF + n_alts ≥ 2 ALTs) the natural
generalisation is the Multinomial bottleneck on the category vector
`p_t = (p_REF, p_ALT1, …, p_ALTn)`, which preserves both moments
component-wise and reproduces the Multinomial covariance structure.

### 3.2 Kimura diffusion limit (continuous bottleneck)

Letting `Ne → ∞` while keeping `t/Ne` fixed gives the Kimura diffusion.
For a **single generation** this collapses to a closed-form Beta
distribution at the same first two moments:

```txt
p_{t+1} | p_t  ~  Beta( α = p_t·(Ne−1),  β = (1−p_t)·(Ne−1) )
```

Direct calculation from the Beta variance formula confirms the moments:

```txt
α + β               =  Ne − 1
E[Beta(α, β)]       =  α / (α + β)              =  p_t                 ✓
Var(Beta(α, β))     =  α·β / [ (α+β)² · (α+β+1) ]
                    =  p_t(1−p_t)·(Ne−1)² / [ (Ne−1)² · Ne ]
                    =  p_t (1 − p_t) / Ne                              ✓
```

So at the per-generation level the discrete Multinomial sampler and
the continuous Beta sampler are **statistically equivalent in the
first two moments**; they differ only in higher moments (the discrete
sampler has heavier tails / discreteness on the grid `k/Ne`).

For multi-allelic sites the analogous continuous law is Dirichlet:

```txt
(p_{t+1,REF}, p_{t+1,ALT1}, …, p_{t+1,ALTn})
        ~  Dirichlet( (Ne−1) · (p_REF, p_ALT1, …, p_ALTn) )
```

Each marginal still has the variance `p_t(1−p_t)/Ne`.

### 3.3 Multi-generation drift accumulation

Iterating the transition kernel `g` times (independent generations,
same Ne per step) gives the standard Wright-Fisher result:

```txt
E [ p_g | p_0 ]    =  p_0
Var( p_g | p_0 )   =  p_0 (1 − p_0) · F_g

where     F_g  =  1 − (1 − 1/Ne)^g       (cumulative drift coefficient)
```

`F_g` is exactly the Kimura **`b`-deficit** (with `b_Kimura = 1 − F_g`).
Limits:

| g           | F_g                         | Interpretation                           |
| :---:       | :---                         | :---                                     |
| `g = 1`     | `1/Ne`                       | one generation of drift                  |
| `g → ∞`     | `1`                          | full fixation                            |
| `Ne · g → ∞`| `1 − exp(−g/Ne)` (continuous)| diffusion-time approximation             |

**Key consequence for estimators.** All three estimators in
§2 fit a **single-generation** model. If you feed them pairs that
actually span `g` generations, they recover the *apparent*
per-generation Ne implied by the cumulative drift:

```txt
Ne_apparent  =  1 / F_g  =  1 / [ 1 − (1 − 1/Ne)^g ]   ≤   Ne
```

so multi-generation pedigrees **deflate** the reported Ne.
Numerically:

| True Ne | g=1   | g=2  | g=3  | g=5  | g=10 |
| :---:   | :---: | :---:| :---:| :---:| :---:|
| 10      | 10.00 | 5.26 | 3.69 | 2.44 | 1.54 |
| 30      | 30.00 | 15.25| 10.25| 6.39 | 3.48 |

This is why pedigree depth matters: a study reporting "Ne ≈ 3" from
grandmother→grandchild pairs may be entirely consistent with a true
per-generation Ne ≈ 10.

### 3.4 What `synthesize.py` actually does (step by step)

For each `(site, mother-child pair)` cell, the simulator executes:

```txt
1. Pick mother-side genotype counts:
     k_mom_alt ~ Binomial(m_dp, p_mom_alt)             # read sampling
     k_mom_ref =  m_dp − k_mom_alt
     p_m       =  k_mom_alt / m_dp                     # observed maternal VAF
     (Multi-allelic: per-ALT counts via Multinomial(m_dp, (p_REF, p_ALTs)))

2. Apply the bottleneck g times (default g = 1):
     For step = 1 .. n_generations:
        if --bottleneck-model discrete:
           counts ~ Multinomial(round(Ne), p_curr)
           p_next =  counts / round(Ne)
        if --bottleneck-model continuous:
           p_next ~  Dirichlet( (Ne − 1) · p_curr )
        p_curr := p_next
     p_child_alt = p_curr_alt

3. Optional outlier injection (--outlier-frac):
     With probability outlier_frac, replace p_child_alt with 0 or 1
     (50/50). Simulates NUMTs / mosaicism / heavy genotyping errors.

4. Child-side read sampling:
     k_child_alt ~ Binomial(c_dp, p_child_alt)         # read sampling
     k_child_ref =  c_dp − k_child_alt
     (Multi-allelic: Multinomial(c_dp, (p_REF, p_ALTs)))

5. Optional GT='.' truncation (--missing-gt-rate):
     With probability missing_gt_rate, mask the genotype as missing.
```

The bottleneck step is the only place the population-genetics
parameter `Ne` enters; everything else is read-level noise.
**Verifying that step 2 matches Wright-Fisher in expectation is the
job of `verify_popgen.py`** — see §4.

---

## 4. Verifying the simulator against theory (`verify_popgen.py`)

We don't just *claim* the bottleneck step is Wright-Fisher; we
**prove** it empirically. `verify_popgen.py` draws hundreds of
thousands of `(p_0, p_g)` pairs directly from the bottleneck sampler
(no read noise, no estimator) and compares the empirical mean of

```txt
F_i  :=  (p_g_i − p_0_i)²  /  [ p_0_i · (1 − p_0_i) ]
```

against the theoretical drift coefficient `F_g = 1 − (1 − 1/Ne)ᵍ`
(§3.3). If our sampler is correct, `mean(F_i) → F_g` with statistical
error `~ 1/√N`.

### 4.1 How to run

```bash
python3 tests/data/ne_pipeline/verify_popgen.py
# Optionally: --n-samples 500000 --vaf-low 0.1 --vaf-high 0.9 --seed 42
```

The script imports `synthesize.gamma_sample` directly so the same
Beta-via-Gamma sampler used to generate cohorts is the one being
tested (no copy-paste).

### 4.2 Reference output

```txt
# Wright-Fisher one- and multi-generation moment check
# Theory: E[F] = F_g = 1 − (1 − 1/Ne)^g
# n_samples = 200000, p_0 ~ Uniform[0.15, 0.85], seed = 20260528

Model        True_Ne   g     E[F]_emp    E[F]_se   F_g_theory   Ne_app   F_relerr   Mean_drift
discrete        3.00   1     0.333296   ±0.00095   0.333333      3.000   −0.011%   −0.00036
discrete        5.00   1     0.198935   ±0.00060   0.200000      5.027   −0.533%   +0.00010
discrete       10.00   1     0.099524   ±0.00031   0.100000     10.048   −0.476%   +0.00042
continuous      3.00   1     0.334348   ±0.00096   0.333333      2.991   +0.304%   +0.00022
continuous      5.00   1     0.201055   ±0.00061   0.200000      4.974   +0.528%   −0.00015
continuous      7.50   1     0.133382   ±0.00041   0.133333      7.497   +0.037%   −0.00007
continuous     30.00   1     0.033335   ±0.00011   0.033333     29.999   +0.004%   −0.00005
continuous     10.00   2     0.189419   ±0.00057   0.190000      5.279   −0.306%   −0.00013
continuous     10.00   5     0.408916   ±0.00115   0.409510      2.445   −0.145%   +0.00037
continuous     10.00  10     0.650750   ±0.00166   0.651322      1.537   −0.088%   +0.00019
discrete       10.00   2     0.188813   ±0.00056   0.190000      5.296   −0.625%   −0.00030
discrete       10.00   5     0.409515   ±0.00115   0.409510      2.442   +0.001%   +0.00021
continuous     30.00   5     0.155985   ±0.00047   0.155920      6.411   +0.042%   −0.00023
continuous     30.00  10     0.287798   ±0.00084   0.287529      3.475   +0.094%   +0.00026

[OK] All rows within 4 SE of Wright-Fisher prediction.
```

### 4.3 Reading the table

* **`E[F]_emp`** — empirical mean of `F_i` across `n_samples` draws.
* **`E[F]_se`** — standard error of that mean (≈ `σ_F / √N`).
* **`F_g_theory`** — predicted drift coefficient `1 − (1 − 1/Ne)ᵍ`.
* **`F_relerr`** — relative error `(emp − thy) / thy`. We expect
  `|F_relerr|` ≲ a few `1/√N` ≈ 0.2–0.7% for `n_samples = 200,000`.
* **`Ne_app`** — apparent Ne implied by the empirical drift,
  `1 / E[F]_emp`. For `g = 1` this should match the true Ne; for
  `g > 1` it equals `1 / F_g` (the deflated apparent Ne).
* **`Mean_drift`** — empirical `E[p_g − p_0]`. The martingale
  property predicts this is exactly 0; we observe ≲ 5e-4 (consistent
  with `1/√N` noise on a bounded variable).

### 4.4 What this proves

* **One-generation correctness.** Rows `g = 1` confirm that both
  samplers reproduce `Var(p_g | p_0) = p_0(1−p_0)/Ne` to within
  statistical noise — the moments the discrete and continuous MLEs
  rely on are correct *by construction* in the simulator.
* **Multi-generation correctness.** Rows `g > 1` confirm that
  *iterating* the per-step sampler reproduces the cumulative
  Wright-Fisher drift `F_g`. This is non-trivial: it would fail if,
  e.g., the continuous sampler accidentally re-injected the original
  `p_0` instead of carrying state forward, or if numerical underflow
  killed mass at the boundaries.
* **Discrete and continuous match each other (g > 1).** At `g = 2,
  Ne = 10` both models give `E[F] ≈ 0.189` (theory: 0.190); at
  `g = 5, Ne = 10` both give `E[F] ≈ 0.409` (theory: 0.4095). They
  diverge only in higher moments (the discrete sampler is
  fingerprinted by its grid, see §6.2 below).
* **Martingale property.** All `Mean_drift` values are within a few
  ×10⁻⁴ of zero, confirming `E[p_g − p_0] = 0` to ~6 σ.

If any row in the panel exceeded **4 SE** from theory the script
exits non-zero and lists the offending row; this can be wired into
CI to catch sampler regressions.

---

## 5. Quick smoke test

Build `mitoquest` first (`cmake --build build -j8` from the repo root),
then run the one-shot demo:

```bash
bash tests/data/ne_pipeline/run_demo.sh
```

Expected (numbers are deterministic):

```txt
==> trans-prep:  cohort.vcf + cohort.fam -> cohort.trans_prep.run.tsv
[trans-prep] Processed 12 variant records.
[trans-prep] Wrote 280 rows (...)

==> ne-estimate (with --cross-check kimura):
[ne-estimate] Optimal Ne ≈ 2–3 (95% CI), max logL ≈ -1625
[ne-estimate] Kimura cross-check: b ≈ 0.79, Ne_Kimura ≈ 4.7
```

The default smoke-test cohort uses the discrete bottleneck simulator
(historical default), so this single result is *not* a fair
comparison — see §6 for the full validation suite.

---

## 6. Three-estimator validation suite

```bash
# From the repo root, after building `mitoquest`:
python3 tests/data/ne_pipeline/validate_estimators.py --clean
```

The harness generates **11 scenarios** — 7 single-generation
scenarios that probe model misspecification and noise robustness,
plus 4 multi-generation scenarios that quantify pedigree-depth bias.
Cohorts are written to `scenarios/<name>/`; raw JSON results are
kept; the bulky VCF/FAM/TSV files are removed by default (pass
`--keep-cohorts` to preserve them).

### 6.1 Scenarios

**Single-generation scenarios** (`g = 1`):

| Scenario | Sim_model | True_Ne | n_pairs | depth | Outliers | What it tests |
| :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **S1** baseline_discrete | discrete | 5 | 200 | 2000 | 0% | Hard bottleneck (matches `--model discrete`) — discrete MLE should win, continuous MLE is misspecified. |
| **S2** baseline_continuous | continuous | 5 | 200 | 2000 | 0% | Beta-diffusion bottleneck (matches `--model continuous`) — continuous MLE should win, discrete MLE is misspecified. |
| **S3** fractional_ne | continuous | **7.5** | 200 | 2000 | 0% | Non-integer truth — only continuous MLE can return a fractional Ne. |
| **S4** small_cohort | continuous | 5 | **20** | 2000 | 0% | Few pairs → CIs widen everywhere; tests stability with small N. |
| **S5** low_depth | continuous | 5 | 200 | **100** | 0% | Low sequencing depth → sampling correction unreliable. |
| **S6** outliers_5pct | continuous | 5 | 200 | 2000 | **5%** | NUMT/mosaicism/error simulation — extreme-drift contamination. |
| **S7** large_ne | continuous | **30** | 200 | 2000 | 0% | Loose bottleneck → all three estimators should converge. |

**Multi-generation scenarios** (`g > 1`, theoretical apparent Ne is the deflated `1/F_g`):

| Scenario | Sim_model | True_Ne | g | Ne_app (theory) | What it tests |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **S8** g2_ne10  | continuous | 10 | 2  | **5.26** | Mild pedigree depth (e.g. grandmother→grandchild). |
| **S9** g5_ne10  | continuous | 10 | 5  | **2.44** | Deeper pedigree — drift halves the apparent Ne. |
| **S10** g10_ne30| continuous | 30 | 10 | **3.48** | Severe accumulation: a "loose" Ne=30 looks like a tight Ne~3.5. |
| **S11** g3_ne10_discrete | discrete | 10 | 3 | **3.69** | Multi-gen *discrete* sampler — keeps the `k/Ne` grid fingerprint. |

### 6.2 Reference results

Reproducible verbatim by `python3 validate_estimators.py --clean`:

| Scenario | True_Ne | g | Ne_app (thy) | MLE_cont (CI)         | MLE_disc (CI)       | Ne_Kimura (CI)        |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| S1 baseline_discrete    | 5    | 1  | 5.00  | 2.55 (2.49–2.61)     | **5.00** (5–5)     | 4.99 (4.76–5.23)    |
| S2 baseline_continuous  | 5    | 1  | 5.00  | **5.12** (4.97–5.28) | 31 (31–31)         | 5.06 (4.86–5.25)    |
| S3 fractional_ne        | 7.5  | 1  | 7.50  | **7.48** (7.23–7.74) | 31 (31–31)         | 7.45 (7.16–7.72)    |
| S4 small_cohort         | 5    | 1  | 5.00  | **5.24** (4.75–5.78) | 36 (36–36)         | 5.40 (4.76–6.16)    |
| S5 low_depth            | 5    | 1  | 5.00  | **4.72** (4.57–4.87) | 11 (11–11)         | 4.73 (4.54–4.93)    |
| S6 outliers_5pct        | 5    | 1  | 5.00  | 3.19 (3.11–3.27)     | 29 (29–29)         | 3.91 (3.75–4.11)    |
| S7 large_ne             | 30   | 1  | 30.00 | **29.27** (28.06–30.51) | 37 (37–37)      | 29.71 (28.46–31.06) |
| S8 g2_ne10              | 10   | 2  | 5.26  | 4.91 (4.76–5.06)     | 32 (32–32)         | **5.30** (5.07–5.52) |
| S9 g5_ne10              | 10   | 5  | 2.44  | 1.90 (1.87–1.93)     | 27 (27–27)         | **2.36** (2.27–2.43) |
| S10 g10_ne30            | 30   | 10 | 3.48  | 2.66 (2.60–2.72)     | 27 (27–27)         | **3.50** (3.38–3.62) |
| S11 g3_ne10_discrete    | 10   | 3  | 3.69  | 2.25 (2.20–2.29)     | **10** (10–10)     | 3.81 (3.66–3.97)    |

Bold = closest to the relevant ground truth (true Ne for `g=1`,
apparent Ne_app=`1/F_g` for `g>1`, or — for S11 — true per-generation
Ne, see §6.4).

### 6.3 Single-generation interpretation (S1–S7)

#### Discrete MLE catastrophically inflates on continuous data (S2–S7)

In every scenario where the simulator emits Beta-diffusion
heteroplasmy, `--model discrete` returns Ne much higher than truth
(often 25–35 for true Ne = 5–7.5). The discrete grid `{k/Ne}` cannot
represent intermediate `p_child` values that the diffusion produces
freely; the optimizer responds by pushing Ne up so the grid becomes
dense enough to fit the observations. **The discrete model is only
appropriate when the bottleneck is genuinely a hard sampling event
(virus-passage experiments, animal breeding bottlenecks).** mtDNA
inheritance is much closer to a continuous diffusion process.

#### Continuous MLE under-estimates on hard discrete data (S1)

When the simulator does pull `round(true_ne)` integer copies (S1,
matching `--model discrete`), continuous MLE returns Ne ≈ 2.55 against
truth Ne = 5. The Beta-diffusion likelihood places mass on
intermediate `p_child` values, but the data clusters at `{k/Ne}`
points — explaining the data with a Beta requires a smaller Ne to
concentrate mass near the modes. Symmetric to the discrete-on-
continuous failure above.

#### Continuous MLE and Kimura agree closely on clean continuous data

In S2, S3, S4, S5, S7 the continuous MLE point estimate is within
~3% of Kimura's. This is the regime in which both estimators are
correctly specified; their CIs overlap heavily. **Tight agreement
between MLE_cont and Kimura is therefore a strong sanity signal.**

#### Outliers pull continuous MLE down faster than Kimura (S6)

With 5% extreme-drift contamination the gap reverses: MLE_cont = 3.19,
Ne_Kimura = 3.91. Why does the MLE move faster than the
method-of-moments estimator?

* Under Beta-diffusion, the probability that `p_child` lands at exactly
  0 or 1 vanishes for any moderate Ne. Even a handful of forced-
  boundary outliers among 4,000 sites carry enormous likelihood weight,
  so the MLE shrinks Ne until the Beta puts non-negligible mass near
  0 and 1.
* Kimura is a *weighted average* of squared drifts; one extreme drift
  contributes at most `(p_m)²` ≤ 1 to the numerator, capped by the
  weight `w_i = p_m(1−p_m)` ≤ 0.25. The pull is bounded.

**Practical consequence on real data:** a noticeable gap with
**MLE_cont < Ne_Kimura** is a fingerprint of high-drift outlier
contamination (NUMTs, mosaicism, sequencing errors). When you see
this, run with `--kimura-trim 0.10 --top-drift-k 20` to inspect and
optionally drop the worst pairs.

#### Loose bottlenecks are easy (S7)

At true Ne = 30, the Beta diffusion is tight (low variance) and
`d_i` is small. MLE_cont, MLE_disc-rounded-up, and Kimura all land
within 20% of truth. When biological intuition says Ne is loose,
disagreement among the three estimators should be small.

### 6.4 Multi-generation interpretation (S8–S11)

These scenarios feed `g`-generation pairs to estimators that all
assume `g = 1`. Wright-Fisher (§3.3) predicts they should report the
*apparent* `Ne_app = 1 / [1 − (1 − 1/Ne)ᵍ]`, which is strictly
smaller than the true per-generation Ne.

#### Kimura tracks the apparent Ne almost exactly (S8–S10)

In every continuous multi-gen scenario, Ne_Kimura matches the
theoretical apparent Ne within statistical noise:

| Scenario | Ne_app theory | Ne_Kimura observed | Relative error |
| :---:    | :---:         | :---:              | :---:          |
| S8       | 5.26          | 5.30               | +0.8%          |
| S9       | 2.44          | 2.36               | −3.3%          |
| S10      | 3.48          | 3.50               | +0.6%          |
| S11      | 3.69          | 3.81               | +3.3%          |

This is exactly what Kimura is designed to estimate: the
**variance-implied Ne**. Its method-of-moments form `Ne = 1/V` is
dimensionally identical to `1/F_g`, with no shape assumption beyond
"variance is `p(1−p)·F`". So when the data carry `g` generations
worth of drift, Kimura faithfully reports `1/F_g`.

#### Continuous MLE further undershoots `Ne_app` (S8–S10)

For the continuous-bottleneck multi-gen scenarios:

| Scenario | Ne_app theory | MLE_cont observed | Undershoot |
| :---:    | :---:         | :---:             | :---:      |
| S8       | 5.26          | 4.91              | −7%        |
| S9       | 2.44          | 1.90              | −22%       |
| S10      | 3.48          | 2.66              | −24%       |

Why does the MLE go *below* the apparent Ne? Iterating the Beta
diffusion `g` times yields a distribution with **heavier tails than
a single-step Beta** at the same variance — extra excess kurtosis
accumulates with `g`. The continuous MLE is fitting a single Beta to
this heavier-tailed reality and interprets the heavy tails as
evidence of even *more* drift, so it pushes Ne further down. This
is consistent with the §6.3 outlier discussion: anything that puts
extra weight near the boundaries pulls MLE_cont down faster than it
pulls Kimura down.

**Practical takeaway.** When `g > 1` is suspected, **prefer
Ne_Kimura over MLE_cont** as the primary estimate of the apparent Ne.
A consistent gap with `MLE_cont < Ne_Kimura` is therefore *also* a
fingerprint of pedigree-depth bias (in addition to the outlier
fingerprint of §6.3); the two are distinguishable by inspecting the
empirical drift distribution and the cohort metadata.

#### Discrete MLE on iterated discrete data: surprising integer recovery (S11)

S11 is the most surprising row. With `g = 3` and `true_ne = 10` under
the **discrete** sampler:

* Theoretical apparent Ne (variance-based) = `1 / F_3 ≈ 3.69`.
* Ne_Kimura recovers the variance-based 3.81 (correct).
* MLE_cont = 2.25 (further down, consistent with the multi-gen
  Beta-misspecification above).
* **MLE_disc = 10**, exactly matching the *true per-generation Ne*.

This is not a coincidence. The discrete sampler iterated `g` times
still produces post-bottleneck VAFs on the integer grid `{0, 1/Ne,
2/Ne, …, 1}` (each step rounds to multiples of `1/Ne`). The discrete
MLE doesn't just measure variance — it also picks up this **grid
fingerprint** in `p_child`, which the BetaBinomial-Binomial
likelihood scores as decisive evidence for Ne = 10 even though the
data also exhibit `g = 3` worth of accumulated variance. So:

* **For real mtDNA (continuous diffusion), MLE_disc is misspecified
  and the grid fingerprint is absent** — see S2, S3, S5, S7 where
  MLE_disc inflates to densify the grid against intermediate values.
* **For genuine hard-bottleneck data with iterated transmission**,
  MLE_disc can in principle recover the per-generation Ne
  *independently* of pedigree depth, while MLE_cont and Kimura
  recover only the variance-based apparent Ne. This is a niche but
  real distinction.

### 6.5 Decision rules

| Situation | Recommended action |
| :--- | :--- |
| mtDNA mother-child cohort, deep sequencing, `g = 1` | **`--model continuous --cross-check kimura`** (default). Trust MLE_cont; use Kimura agreement as a sanity check. |
| `MLE_cont` and `Ne_Kimura` agree within ~10% | High confidence — report MLE_cont. |
| `MLE_cont < Ne_Kimura` by 15–25% | Possible **outlier contamination** *or* **pedigree depth `g > 1`**. Re-run with `--kimura-trim 0.10 --top-drift-k 20`; inspect the per-pair generation gaps in your FAM. If outliers are the cause the gap shrinks after trimming; if pedigree depth is the cause Kimura still tracks `1/F_g`. |
| `MLE_cont < Ne_Kimura` by >25% with no obvious outliers | Suspect **deep pedigrees** (`g ≥ 3`). Treat **Ne_Kimura as the apparent Ne**; if you can estimate `g` from the FAM, recover the per-generation Ne via `Ne = 1 / [1 − (1 − F_g)^{1/g}]` where `F_g = 1 / Ne_Kimura`. |
| `MLE_cont > Ne_Kimura` by >20% | Possible heavy missingness or boundary clipping; check `Pairs_Used` and `CI_*_Clipped` in JSON. |
| Hard sampling bottleneck (virus passage, animal breeding) | `--model discrete` is more appropriate; continuous would under-estimate. **Bonus**: under iterated hard bottlenecks, discrete MLE can recover the per-generation Ne even when `g > 1` (S11). |
| Loose bottleneck (Ne ≳ 30) and `g = 1` | All three estimators converge — any of them is fine; continuous is still preferred for mtDNA. |
| Low depth (`m_dp` or `c_dp` < 200) | Continuous MLE remains accurate; discrete MLE biases up; Kimura's sampling correction becomes noisy. |
| Few transmission pairs (n < 30) | Use continuous MLE and report wide CIs honestly. Kimura bootstrap CI is a useful agreement check. |

### 6.6 When discrete MLE is actually the right choice

Reading §6.2 it is tempting to conclude that **discrete MLE is
useless** — it loses on 9 of 11 scenarios and only matches truth on
S1 and S11. That conclusion is wrong, but it deserves a careful
answer because the validation panel is **deliberately
mtDNA-flavoured**: 9 of 11 scenarios use the continuous Beta-diffusion
bottleneck because that is how human mtDNA actually behaves. On any
continuous-bottleneck data the discrete model is *misspecified by
design*, and a misspecified MLE always loses. So:

> **For human mtDNA deep-sequencing studies, `--model continuous` is
> the primary estimator and discrete MLE is rarely the right answer.**
> mitoquest still ships discrete MLE because it has four legitimate,
> non-overlapping niches in which it is *the* correct choice.

#### Niche 1 — Hard biological bottlenecks (its native habitat)

There are real systems where the Multinomial sampling kernel is the
*biology*, not a modelling artefact:

| System | Why discrete is the correct model |
| :--- | :--- |
| Virus passage experiments (influenza, HIV, SARS-CoV-2 in cell culture) | Each passage really does pick `Ne` virions; founder count is integer. |
| Animal / plant breeding bottlenecks | Founder population is an integer number of individuals. |
| Bacterial colony picks, microfluidic dilution, FACS sub-sampling | Counts are integer cells. |
| Y-chromosome / single-locus founder events | Per-generation transmission is a single chromosome out of a discrete pool. |

For these systems, MLE_cont systematically *under*-estimates Ne
(mirror image of S1) because Beta cannot reproduce the `{k/Ne}` grid
clumping. Discrete MLE wins on the same generative process that
defines the system.

#### Niche 2 — Very small Ne (Ne ≲ 5)

The Beta-diffusion law is the `Ne → ∞` limit of Wright-Fisher. At
small Ne the limit is poor:

* The discrete sampler can put non-zero mass on `p_child = 0` and
  `p_child = 1` (true loss / fixation events).
* `Beta(α, β)` with `α + β = Ne − 1` has zero density at `0` and `1`
  → the continuous likelihood has to compensate by squeezing `Ne`
  even smaller, which is exactly the S1 failure (true Ne=5,
  MLE_cont=2.55).

For mtDNA-like data with deep coverage and `Ne` in the 5–30 range,
this hardly matters; for cell-line bottleneck experiments where `Ne`
might literally be 2 or 3, it matters a lot.

#### Niche 3 — Multi-generation per-generation Ne recovery (S11)

This is the most striking row in the validation table. Under iterated
**discrete** transmission with `g = 3, true_ne = 10`:

| Estimator | Reports | Recovers |
| :--- | :---: | :--- |
| MLE_cont    | 2.25 | neither apparent Ne (3.69) nor per-gen Ne (10) |
| Ne_Kimura   | 3.81 | apparent Ne (variance-implied) ≈ `1/F_g` |
| **MLE_disc**| **10** | **true per-generation Ne, despite 3 generations of drift** |

Reason (see §6.4): each iterated discrete step keeps the
post-bottleneck VAFs on the integer grid `{0, 1/Ne, …, 1}`, *regardless
of `g`*. The BetaBinomial-Binomial likelihood reads that grid
signature as decisive evidence for `Ne = 10`. So when you have
genuine discrete-bottleneck data and a deep pedigree, **discrete MLE
is the only one of the three estimators that can disentangle `Ne`
from `g`** — Kimura and continuous MLE both report the deflated
apparent Ne `1/F_g`. This is a niche but real biological use case for
passage / breeding studies with multi-generation pedigrees.

#### Niche 4 — Diagnostic for model misspecification

Even when you intend to report MLE_cont, the **gap** between MLE_cont
and MLE_disc is itself informative:

| Pattern | Diagnosis |
| :--- | :--- |
| `MLE_disc ≈ MLE_cont`, both small | Tight bottleneck, both models converge — high confidence (S7). |
| `MLE_disc ≫ MLE_cont` (e.g. 30 vs 5) | Continuous diffusion is the right physical model (mtDNA, S2–S10). The discrete grid had to inflate to fit. |
| `MLE_disc ≪ MLE_cont` | Hard discrete bottleneck (S1-flavoured). Switch primary estimate to discrete. |
| `MLE_disc` lands close to true Ne while MLE_cont and Kimura land at `1/F_g` | Iterated discrete sampling — flag pedigree depth (S11). |

You cannot read these patterns off a single estimate; you need both
numbers. Removing discrete MLE from the toolbox would lose the
ability to detect (a) model misspecification, (b) whether the
bottleneck is really discrete vs diffusion, and (c) the niche-3
pedigree-depth recovery.

#### Bottom line

* **Default for human mtDNA = `--model continuous --cross-check
  kimura`.** Discrete MLE is *not* a competitor here; it is
  intentionally misspecified for biologically continuous data.
* **Default for hard-bottleneck systems = `--model discrete`.**
  mitoquest exposes both because the same binary needs to handle both
  regimes.
* **Always look at all three numbers.** They each estimate a
  *different* quantity (likelihood under Beta, likelihood under
  Multinomial, variance-implied Ne); their pattern of agreement and
  disagreement is itself a primary diagnostic for misspecification,
  outlier contamination, and pedigree-depth bias.

---

## 7. `synthesize.py` reference

Standalone deterministic simulator (Python stdlib only). Supports
both bottleneck models, optional outlier injection, and
multi-generation transmission so you can reproduce `--model
continuous` / `--model discrete` ground truth at any pedigree depth.

```bash
# Default discrete bottleneck (backward compatible):
python3 synthesize.py --true-ne 5 --n-pairs 20

# Continuous Beta-diffusion bottleneck, allows fractional true_ne:
python3 synthesize.py --bottleneck-model continuous --true-ne 7.5

# Inject 5% high-drift outliers to test estimator robustness:
python3 synthesize.py --bottleneck-model continuous --true-ne 5 \
                     --outlier-frac 0.05

# Multi-generation pedigree (drift accumulates as F_g = 1−(1−1/Ne)^g):
python3 synthesize.py --bottleneck-model continuous --true-ne 10 \
                     --n-generations 5

# Tighter cohort, lower depth, custom seed:
python3 synthesize.py --true-ne 3 --n-pairs 100 --m-dp 200 --c-dp 200 \
                     --seed 42

# Pure biallelic dataset (no multi-allelic sites):
python3 synthesize.py --n-multiallelic 0

# Stress-test multi-allelic handling (every site tri-allelic):
python3 synthesize.py --n-multiallelic 12
```

Key flags:

| Flag                       | Default      | Notes                                                       |
| :---:                      | :---:        | :---                                                        |
| `--bottleneck-model`       | `discrete`   | `discrete` (Multinomial, integer Ne) or `continuous` (Beta-diffusion / Dirichlet, fractional Ne allowed). |
| `--true-ne`                | `5`          | Float; rounded to int under the discrete model.             |
| `--n-generations`          | `1`          | Number of bottleneck steps between mother and child. `g=1` = direct, `g≥2` = grandmother→grandchild and deeper. Drift accumulates as `F_g = 1−(1−1/Ne)^g`. |
| `--n-pairs`                | `20`         | Number of mother-child trios.                               |
| `--n-sites`                | `12`         | Number of mtDNA SNV sites.                                  |
| `--m-dp` / `--c-dp`        | `2000`       | Per-sample read depth.                                      |
| `--vaf-low` / `--vaf-high` | `0.15`/`0.85`| Maternal VAF range (Uniform).                               |
| `--n-multiallelic`         | `2`          | How many of `n_sites` should be tri-allelic.                |
| `--missing-gt-rate`        | `0.05`       | Probability of `GT='.'` truncation per (sample, site).      |
| `--outlier-frac`           | `0.0`        | Fraction of (site, pair) cells forced to `p_child ∈ {0, 1}`. |
| `--seed`                   | `20260528`   | RNG seed (deterministic output).                            |

Multi-allelic decomposition: `mitoquest trans-prep` splits each
multi-allelic VCF record into one row per ALT × per trio. For the
default 20-pair cohort this gives `10 biallelic × 20 × 1 + 2
tri-allelic × 20 × 2 = 280` rows. `ne-estimate` then treats each row
as an independent biallelic `(REF reads, ALT reads)` observation.

---

## 8. Reproducibility & version

* `synthesize.py` and `verify_popgen.py` are fully deterministic
  (`--seed`, default `20260528`); no third-party dependencies.
* `validate_estimators.py` calls `synthesize.py` and `mitoquest
  ne-estimate` only — no network, no external data, no Python packages
  beyond the stdlib.
* The harness was developed and validated against `mitoquest v1.8.3`.
* For an in-tree numerical correctness test, see
  [`tests/test_ne_estimate.cpp`](../../test_ne_estimate.cpp), which
  contains `NeEstContinuous.RealValuedOptimum` and
  `NeEstContinuous.RealValuedAgreesWithKimura` (Beta-diffusion
  ground truth, 800–1000 simulated pairs).

## 9. Notes & caveats

* All cohorts here are **synthetic** — positions are spread along
  `chrM` for plausibility but are not real mtDNA variants. Do not use
  any of these VCFs for biological interpretation.
* The discrete bottleneck simulator (`--bottleneck-model discrete`)
  matches the assumption used by `--model discrete` MLE exactly; the
  continuous simulator (`--bottleneck-model continuous`) matches the
  assumption used by `--model continuous`. Comparing each estimator
  against its matching ground truth is the cleanest way to read §6.
* For real data, mtDNA biology is much closer to the continuous
  Beta-diffusion model: the `--model continuous` MLE is the
  recommended primary estimate, and the Kimura cross-check is the
  recommended diagnostic.
* All three estimators currently fit a **single-generation** model.
  When pedigrees are deeper than mother→child, they recover the
  *apparent* Ne `1/F_g`, not the per-generation Ne (§3.3, §6.4).
  Multi-generation explicit modelling is on the roadmap.
