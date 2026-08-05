# MitoQuest v1.9.2 — `ne-estimate` Kimura cross-check stderr fixes

This is a patch release that fixes three bugs in the `mitoquest ne-estimate`
subcommand's stderr reporting of the Wonnapinij/Kimura cross-check when
`--kimura-trim` and/or `--top-drift-k` are used.

---

## Bug fixes

### 1. WARNING always compared MMLE against *untrimmed* Kimura

The `>3x disagreement` WARNING at the end of `run()` always used
`r.kimura.ne_kimura` (the untrimmed value), even when `--kimura-trim`
was set. This made the WARNING message identical with or without
`--kimura-trim`, giving the false impression that trimming had no effect.

**Fix:** when a trimmed Kimura estimate is available, the WARNING now
compares `Ne_MMLE` against `Ne_Kimura_trimmed` instead of the untrimmed
value. If the trimmed value still disagrees, the message explicitly names
`Ne_Kimura_trimmed` and suggests increasing the trim fraction.

### 2. Recommendation always said "re-run with --kimura-trim"

The WARNING text unconditionally recommended re-running with
`--kimura-trim 0.10`, even when the user had already applied trimming.

**Fix:** when trim was already applied and reconciled the two estimators
(trimmed Kimura now within 3x of MMLE), a `NOTE` is printed instead:

```
[ne-estimate] NOTE: trimming 20% of high-drift pairs reconciled
    Ne_Kimura_trimmed (...) with Ne_MMLE (...).
```

### 3. Trimmed Kimura and top-K outliers were invisible in stderr

The trimmed Kimura line was printed with a misplaced `\n` prefix, and the
top-K drift outliers (`--top-drift-k`) were written only to the JSON output
— never to stderr.

**Fix:**
- Trimmed Kimura now prints cleanly on its own line with `b_trimmed` and
  `Ne_Kimura_trimmed`.
- Top-K drift outliers are printed to stderr (pair index, M_VAF, C_VAF, F_i)
  for quick inspection without JSON parsing.

### 4. Silent no-op when Kimura options used without `--cross-check kimura`

`--kimura-trim`, `--top-drift-k`, and `--kimura-bootstrap` were silently
ignored when `--cross-check kimura` was not specified.

**Fix:** explicit WARNING is now emitted when any of these options are
set without `--cross-check kimura`.

---

## Example output (with fix)

```
[ne-estimate] Fitting Ne on 1234 pairs (maternal VAF in [0.01, 0.999], model=continuous).
[ne-estimate] Optimal Ne = 54.8158 (95% CI: 48.12 - 62.34), max marginal logL = -12345.6
[ne-estimate] Kimura cross-check: b = 0.80215432, Ne_kimura = 5.05325 on 987 informative pairs
[ne-estimate] Trimmed Kimura (trim 20%): b_trimmed = 0.97234, Ne_kimura_trimmed = 36.12 on 789 pairs
[ne-estimate] Top-20 drift outlier pairs (by F_i descending):
  #1  pair_index=42  M_VAF=0.123  C_VAF=0.891  F_i=0.98234
  #2  pair_index=891 M_VAF=0.456  C_VAF=0.012  F_i=0.95123
  ...
[ne-estimate] NOTE: trimming 20% of high-drift pairs reconciled
    Ne_Kimura_trimmed (36.12) with Ne_MMLE (54.8158).
```

---

## Full changelog

- **fix**: WARNING now compares MMLE against trimmed Kimura when available
- **fix**: reconciliation NOTE printed when trimming resolves the MMLE/Kimura gap
- **fix**: trimmed Kimura result printed cleanly to stderr with `b_trimmed`
- **fix**: top-K drift outliers printed to stderr (were JSON-only)
- **fix**: WARNING emitted when Kimura-specific options used without
  `--cross-check kimura`
- **build**: bump version 1.9.1 → 1.9.2

---

## Binaries

| Platform | Download |
|----------|----------|
| Linux (x86_64) | [mitoquest-linux-static](https://github.com/ShujiaHuang/mitoquest/releases/download/v1.9.2/mitoquest-linux-static) — requires glibc ≥ 2.35 |
| macOS (arm64 / Intel) | [mitoquest-macos-static](https://github.com/ShujiaHuang/mitoquest/releases/download/v1.9.2/mitoquest-macos-static) — requires macOS 12+ |

---

## Upgrading from v1.9.1

No breaking changes. The primary MMLE estimate (`Ne`, `CI_95_Low`,
`CI_95_High`) is unchanged. Only the stderr summary of the Kimura
cross-check is improved; the JSON output schema is unchanged (trimmed
fields were already present in v1.9.1).
