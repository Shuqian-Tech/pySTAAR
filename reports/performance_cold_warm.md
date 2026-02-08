# Phase 3 Cold vs Warm Latency

- Generated: 2026-02-08T03:19:17Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy (warm mode): 1 run(s) discarded
- Measured runs: 3 per scenario per mode
- Cold mode includes interpreter startup + imports + single workflow call.
- Warm mode measures repeated calls in the same Python process.

## Results

| Scenario | Cold median (s) | Warm median (s) | Warm speedup vs cold |
|---|---:|---:|---:|
| staar_unrelated_glm | 1.601065 | 0.215315 | 7.44x |
| staar_related_sparse_glmmkin_pure | 2.433040 | 0.276721 | 8.79x |
| staar_related_sparse_binary_spa_pure | 1.966867 | 0.972656 | 2.02x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 3.244512 | 0.000067 | 48245.53x |

## Artifacts

- `benchmarks/phase3_cold_warm_raw.csv`
- `benchmarks/phase3_cold_warm_summary.csv`
- `benchmarks/phase3_cold_warm_comparison.csv`
- `benchmarks/phase3_cold_warm_meta.json`
