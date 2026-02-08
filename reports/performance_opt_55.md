# Phase 3 Performance Snapshot (STAAR-55)

- Generated: 2026-02-07T16:53:07Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario per mode
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Results by Mode

| Scenario | Baseline no-cache median (s) | Optimized cache median (s) | Speedup |
|---|---:|---:|---:|
| staar_related_sparse_glmmkin_pure | 0.968130 | 0.192804 | 5.02x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 1.588007 | 0.817542 | 1.94x |

## Comparison Artifact

- `benchmarks/phase3_opt_55_comparison.csv`

## Notes

- Baseline mode disables related-null-model caching by replacing `workflows.fit_null_glmmkin` with a passthrough wrapper.
- Optimized mode restores the original `fit_null_glmmkin` and uses the new cached related-null-model path.
- Sentinels are expected to match between modes; values are recorded in `benchmarks/phase3_opt_55_meta.json`.
