
## mitoquest v1.10.0 — Three-Generation G-M-C Trio Marginal Likelihood for Ne Estimation

### Highlights

This release extends the mtDNA transmission bottleneck estimator to support **three-generation pedigrees** (grandmother → mother → child) alongside the existing two-generation (mother → child) design, enabling mixed-cohort Ne estimation from a single VCF + two FAM files.

---

### `mitoquest trans-prep` — New `--gm-fam` flag

A second PLINK FAM file can now be supplied via `-g/--gm-fam` describing grandmother–mother relationships. For each mother–child pair, the MC mother is looked up as a child in the GM FAM; if matched, the GM FAM's `MOTHER_ID` is the pedigree grandmother.

The output TSV switches to a **wide 22-column format** that includes:

| New Column | Description |
|---|---|
| `GRANDMOTHER_ID` | Grandmother sample ID (`NA` when `HAS_G=0`) |
| `GRANDMOTHER_DP` | Grandmother total depth |
| `GRANDMOTHER_AD_REF` / `GRANDMOTHER_AD_ALT` | Grandmother allele depths |
| `GRANDMOTHER_VAF` | Grandmother ALT VAF |
| `HAS_G` | `1` = G-M-C trio (3-gen); `0` = MC pair (2-gen) |

```bash
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f mc_pairs.fam \
    -g gm_pairs.fam \        # grandmother-mother FAM
    -o trio_pairs.tsv
```

---

### `mitoquest ne-estimate` — Trio marginal likelihood

When the input TSV contains `HAS_G`, `GRANDMOTHER_DP`, and `GRANDMOTHER_AD_ALT` columns, `ne-estimate` auto-detects them and dispatches each trio row to the **closed-form G-M-C trio marginal likelihood**:

$$I_{\text{trio}}(N_e) = \binom{d_M}{k_M} \binom{d_C}{k_C} \cdot \frac{B(\alpha_G + k_M + k_C,\; \beta_G + (d_M - k_M) + (d_C - k_C))}{B(\alpha_G,\; \beta_G)}$$

where $\hat{p}_G = g_{\text{ad\_alt}}/g_{\text{dp}}$, $\alpha_G = \hat{p}_G(N_e - 1)$, $\beta_G = (1 - \hat{p}_G)(N_e - 1)$. The mother's latent heteroplasmy $p_M$ is analytically marginalised out. Homoplasmic grandmothers ($\hat{p}_G \in \{0, 1\}$) contribute a Ne-independent constant (0.0). Mixed cohorts (both `HAS_G=0` and `HAS_G=1` rows) are handled seamlessly via the composite marginal log-likelihood.

```bash
mitoquest ne-estimate \
    -i trio_pairs.tsv \
    --min-ne 2 --max-ne 100 \
    --model continuous \
    -o ne_result.json
```

---

### Validation

- **13 new unit tests** covering trio detection (`TransPrepThreeGen`) and trio likelihood correctness (`NeEstTrio`): closed-form vs. Gauss-Legendre quadrature argmax agreement, `HAS_G=0` fallback, homoplasmic-grandmother edge case, `Ne=1` boundary, manual formula verification, global LL dispatch, and `load_pairs` backward compatibility with legacy 16-column TSVs. All 21 relevant tests PASS.
- **End-to-end smoke test** on a synthetic 3-gen cohort (true $N_e = 20$, 289 informative rows): estimated $N_e = 17.72$ (95% CI 15.14–20.64); the true value falls inside the CI as expected.

---

### Bug Fixes

- Fixed a numerical bug in the Gauss-Legendre quadrature node/weight computation: the derivative $P_n'(z_i)$ was being taken from the Newton iteration residual ($\approx 0$ at convergence), producing infinite weights. Now recomputed explicitly via the recurrence $P_n'(z) = n(zP_n(z) - P_{n-1}(z))/(z^2 - 1)$ after convergence.

---

### Backward Compatibility

- Legacy 16-column TSVs (without `--gm-fam`) are read identically to previous versions; `HAS_G` defaults to `0` when absent.
- All existing CLI flags and JSON output fields are preserved.

---

### Files changed (10 files, +1288/−38 lines)

| File | Change |
|---|---|
| `src/trans_prep.{h,cpp}` | `--gm-fam` flag, `parse_gm_fam()`, `resolve_gm_for_trios()`, wide TSV output |
| `src/ne_estimate.{h,cpp}` | `PairData.g_*` fields, `compute_ll_trio_continuous()`, `compute_ll_trio_quadrature()`, trio dispatch in `compute_global_ll_continuous()`, backward-compatible `load_pairs()` |
| `tests/test_trans_prep.cpp` | 3 new `TransPrepThreeGen` tests |
| `tests/test_ne_estimate.cpp` | 10 new `NeEstTrio` tests |
| `tests/data/ne_pipeline/smoke_trio.tsv` | Synthetic 3-gen smoke test fixture |
| `README.md` | Trio model derivation, `--gm-fam` usage, wide TSV column table |
| `CMakeLists.txt`, `src/version.h` | Version bump 1.9.2 → 1.10.0 |
