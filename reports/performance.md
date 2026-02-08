# Phase 3 Performance (Python vs R)

- Generated: 2026-02-06T07:46:16Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python benchmark logs: `benchmarks/phase3_baseline_raw.csv`, `benchmarks/phase3_baseline_summary.csv`, `benchmarks/phase3_baseline_meta.json`
- R benchmark logs: `benchmarks/phase3_baseline_r_raw.csv`, `benchmarks/phase3_baseline_r_summary.csv`, `benchmarks/phase3_baseline_r_meta.json`
- Cross-language comparison: `benchmarks/phase3_cross_language_comparison.csv`

## Environment

- Python version: `3.13.5`
- R version: `4.5.0`
- STAAR package version: `0.9.8`

## Cross-Language Results

| Scenario | Python median (s) | R median (s) | Python vs R speedup |
|---|---:|---:|---:|
| staar_unrelated_glm | 0.235836 | 0.437000 | 1.85x |
| staar_related_sparse_glmmkin_pure | 1.004456 | 1.625000 | 1.62x |
| staar_unrelated_binary_spa | 0.288611 | 0.085000 | 0.29x |
| staar_related_sparse_binary_spa_pure | 0.715168 | 4.725000 | 6.61x |
| staar_unrelated_glm_cond | 0.302497 | 2.136000 | 7.06x |
| indiv_score_unrelated_glm | 0.120387 | 0.382000 | 3.17x |
| ai_staar_unrelated_glm | 0.928488 | 0.384000 | 0.41x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 1.664256 | 1.624000 | 0.98x |

## Summary

- Geometric mean Python speedup vs R across scenarios: `1.64x`.
- Speedup is computed as `R median / Python median` (values > 1 mean Python is faster).
- This report is the Phase 3 baseline before Python optimizations.

## Optimization Update (2026-02-07)

Targeted Phase 3 optimization work for `STAAR-43` and `STAAR-44` is captured in:

- `benchmarks/phase3_opt_43_44_raw.csv`
- `benchmarks/phase3_opt_43_44_summary.csv`
- `benchmarks/phase3_opt_43_44_meta.json`
- `benchmarks/phase3_opt_43_44_comparison.csv`
- `reports/performance_opt_43_44.md`

Observed median improvements versus the Phase 3 baseline:

- `staar_unrelated_binary_spa`: `0.288611s` -> `0.161430s` (`~1.79x`)
- `ai_staar_unrelated_glm`: `0.928488s` -> `0.856711s` (`~1.08x`)

## Optimization Update (2026-02-07, STAAR-46)

Targeted Phase 3 optimization work for `STAAR-46` is captured in:

- `benchmarks/phase3_post42_probe_raw.csv`
- `benchmarks/phase3_post42_probe_summary.csv`
- `benchmarks/phase3_post42_probe_meta.json`
- `benchmarks/phase3_opt_46_raw.csv`
- `benchmarks/phase3_opt_46_summary.csv`
- `benchmarks/phase3_opt_46_meta.json`
- `benchmarks/phase3_opt_46_comparison.csv`
- `reports/performance_opt_46.md`

Observed median improvement versus the post-`STAAR-42` probe baseline:

- `staar_related_sparse_binary_spa_pure`: `2.343567s` -> `0.909528s` (`~2.58x`)

## Optimization Update (2026-02-07, STAAR-55)

Targeted Phase 3 optimization work for `STAAR-55` is captured in:

- `scripts/run_phase3_opt_55_benchmarks.py`
- `benchmarks/phase3_opt_55_raw.csv`
- `benchmarks/phase3_opt_55_summary.csv`
- `benchmarks/phase3_opt_55_meta.json`
- `benchmarks/phase3_opt_55_comparison.csv`
- `reports/performance_opt_55.md`

Observed median improvements in targeted no-cache vs cache benchmark mode:

- `staar_related_sparse_glmmkin_pure`: `0.968130s` -> `0.192804s` (`~5.02x`)
- `ai_staar_related_sparse_glmmkin_find_weight_pure`: `1.588007s` -> `0.817542s` (`~1.94x`)

## Optimization Update (2026-02-08, STAAR-56)

Targeted Phase 3 optimization work for `STAAR-56` is captured in:

- `scripts/run_phase3_opt_56_benchmarks.py`
- `benchmarks/phase3_opt_56_raw.csv`
- `benchmarks/phase3_opt_56_summary.csv`
- `benchmarks/phase3_opt_56_meta.json`
- `benchmarks/phase3_opt_56_comparison.csv`
- `reports/performance_opt_56.md`

Observed median improvements in targeted baseline-vs-cache repeated-call benchmark mode:

- `ai_staar_related_sparse_glmmkin`: `0.656068s` -> `0.000054s` (`~12047x`)
- `ai_staar_related_sparse_glmmkin_find_weight`: `0.630982s` -> `0.000019s` (`~33065x`)

## Release Hardening Update (2026-02-08, Cold vs Warm)

Cold-process versus warm-inprocess latency measurements are captured in:

- `scripts/run_phase3_cold_warm_benchmarks.py`
- `benchmarks/phase3_cold_warm_raw.csv`
- `benchmarks/phase3_cold_warm_summary.csv`
- `benchmarks/phase3_cold_warm_meta.json`
- `benchmarks/phase3_cold_warm_comparison.csv`
- `reports/performance_cold_warm.md`

Observed warm-process median speedups versus cold-process median:

- `staar_unrelated_glm`: `7.44x`
- `staar_related_sparse_glmmkin_pure`: `8.79x`
- `staar_related_sparse_binary_spa_pure`: `2.02x`
- `ai_staar_related_sparse_glmmkin_find_weight_pure`: `48245.53x` (reflects repeat-call cache hit path)

## Non-Example Coverage Update (2026-02-07, STAAR-45)

Non-`example` runtime-directory performance coverage is captured in:

- `benchmarks/phase3_nonexample_raw.csv`
- `benchmarks/phase3_nonexample_summary.csv`
- `benchmarks/phase3_nonexample_meta.json`
- `reports/performance_nonexample.md`

Observed median runtime for `staar_unrelated_glm_nonexample_dir`: `0.165067s`.
