# Phase 3 Performance Snapshot (STAAR-56)

- Generated: 2026-02-08T02:12:18Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario per mode
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Results by Mode

| Scenario | Baseline no-cache median (s) | Optimized cache median (s) | Speedup |
|---|---:|---:|---:|
| ai_staar_related_sparse_glmmkin | 0.656068 | 0.000054 | 12047.25x |
| ai_staar_related_sparse_glmmkin_find_weight | 0.630982 | 0.000019 | 33065.10x |

## Comparison Artifact

- `benchmarks/phase3_opt_56_comparison.csv`

## Notes

- Baseline mode bypasses AI result caching by wrapping `workflows.ai_staar` with a passthrough function.
- Optimized mode restores the original `ai_staar` callable and enables repeated-call AI result caching.
- Sentinels are expected to match between modes; values are recorded in `benchmarks/phase3_opt_56_meta.json`.
