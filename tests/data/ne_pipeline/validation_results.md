# Three-estimator validation results

Binary: `/Users/huangshujia/Projects/mitoquest/bin/mitoquest`

Each row reports the same `cohort.transmission_pairs.tsv` fed to
`mitoquest ne-estimate` three times: once with `--model continuous`
(MLE_cont), once with `--model discrete` (MLE_disc), with
`--cross-check kimura` enabled in both runs (Ne_Kimura).

Column `g` is the simulated pedigree depth (number of generations
between mother and child). Column `Ne_app_thy` is the Wright-Fisher
prediction of the *apparent* per-generation Ne that any single-
generation estimator should report: `1 / (1 - (1 - 1/Ne)^g)`.
For `g=1` this equals the true Ne; for `g>1` it is strictly smaller.

| Scenario | True_Ne | g | Ne_app_thy | N_pairs | m_dp/c_dp | Sim_model | Outliers | MLE_cont | CI_cont | MLE_disc | CI_disc | Ne_Kimura | CI_Kimura | Pairs_Used | Expected |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| S1_baseline_discrete | 5.00 | 1 | 5.00 | 200 | 2000/2000 | discrete | 0.00 | 2.55 | 2.49-2.61 | 5 | 5.00-5.00 | 4.99 | 4.76-5.23 | 4000 | MLE_discrete |
| S2_baseline_continuous | 5.00 | 1 | 5.00 | 200 | 2000/2000 | continuous | 0.00 | 5.12 | 4.97-5.28 | 31 | 31.00-31.00 | 5.06 | 4.86-5.25 | 4000 | MLE_continuous |
| S3_fractional_ne | 7.50 | 1 | 7.50 | 200 | 2000/2000 | continuous | 0.00 | 7.48 | 7.23-7.74 | 31 | 31.00-31.00 | 7.45 | 7.16-7.72 | 4000 | MLE_continuous |
| S4_small_cohort | 5.00 | 1 | 5.00 | 20 | 2000/2000 | continuous | 0.00 | 5.24 | 4.75-5.78 | 36 | 36.00-36.00 | 5.40 | 4.76-6.16 | 400 | MLE_continuous |
| S5_low_depth | 5.00 | 1 | 5.00 | 200 | 100/100 | continuous | 0.00 | 4.72 | 4.57-4.87 | 11 | 11.00-11.00 | 4.73 | 4.54-4.93 | 3984 | MLE_continuous |
| S6_outliers_5pct | 5.00 | 1 | 5.00 | 200 | 2000/2000 | continuous | 0.05 | 3.19 | 3.11-3.27 | 29 | 29.00-29.00 | 3.91 | 3.75-4.11 | 4000 | MLE_continuous |
| S7_large_ne | 30.00 | 1 | 30.00 | 200 | 2000/2000 | continuous | 0.00 | 29.27 | 28.06-30.51 | 37 | 37.00-37.00 | 29.71 | 28.46-31.06 | 4000 | all converge |
| S8_g2_ne10 | 10.00 | 2 | 5.26 | 200 | 2000/2000 | continuous | 0.00 | 4.91 | 4.76-5.06 | 32 | 32.00-32.00 | 5.30 | 5.07-5.52 | 4000 | MLE_cont matches 1/F_g |
| S9_g5_ne10 | 10.00 | 5 | 2.44 | 200 | 2000/2000 | continuous | 0.00 | 1.90 | 1.87-1.93 | 27 | 27.00-27.00 | 2.36 | 2.27-2.43 | 4000 | MLE_cont matches 1/F_g |
| S10_g10_ne30 | 30.00 | 10 | 3.48 | 200 | 2000/2000 | continuous | 0.00 | 2.66 | 2.60-2.72 | 27 | 27.00-27.00 | 3.50 | 3.38-3.62 | 4000 | MLE_cont matches 1/F_g |
| S11_g3_ne10_discrete | 10.00 | 3 | 3.69 | 200 | 2000/2000 | discrete | 0.00 | 2.25 | 2.20-2.29 | 10 | 10.00-10.00 | 3.81 | 3.66-3.97 | 4000 | all match 1/F_g |

Per-scenario raw JSON files: `scenarios/<name>/ne.continuous.json` and `scenarios/<name>/ne.discrete.json`.

