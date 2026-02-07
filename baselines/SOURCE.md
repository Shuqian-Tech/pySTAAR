# Baseline Source

- **R repo**: https://github.com/xihaoli/STAAR
- **Commit**: 9db9dd504905b9f469146f670e5f6dbe3e08d01a
- **Extraction date**: 2026-02-06
- **Extracted by**: codex
- **Hardware capture**: `baselines/hardware.txt`

## Commands used

```bash
R CMD check . --no-manual 2>&1 | tee baselines/r_cmd_check.log
Rscript baselines/scripts/environment_capture.R
bash baselines/scripts/capture_hardware.sh
Rscript baselines/scripts/fingerprint_example_data.R
Rscript baselines/scripts/prepare_example_data.R
Rscript baselines/scripts/extract_unrelated_glm.R
Rscript baselines/scripts/extract_unrelated_glm_cond.R
Rscript baselines/scripts/extract_unrelated_binary_spa.R
Rscript baselines/scripts/extract_unrelated_binary_spa_filter.R
Rscript baselines/scripts/extract_ai_staar_unrelated.R
Rscript baselines/scripts/extract_ai_staar_unrelated_find_weight.R
Rscript baselines/scripts/extract_ai_staar_related_sparse.R
Rscript baselines/scripts/extract_ai_staar_related_dense.R
Rscript baselines/scripts/extract_ai_staar_related_sparse_find_weight.R
Rscript baselines/scripts/extract_ai_staar_related_dense_find_weight.R
Rscript baselines/scripts/extract_indiv_score_unrelated.R
Rscript baselines/scripts/extract_related_sparse_glmmkin.R
Rscript baselines/scripts/extract_related_dense_glmmkin.R
Rscript baselines/scripts/extract_related_sparse_glmmkin_cond.R
Rscript baselines/scripts/extract_related_dense_glmmkin_cond.R
Rscript baselines/scripts/extract_related_sparse_binary_spa.R
Rscript baselines/scripts/extract_related_dense_binary_spa.R
Rscript baselines/scripts/extract_related_sparse_binary_spa_filter.R
Rscript baselines/scripts/extract_related_dense_binary_spa_filter.R
Rscript baselines/scripts/extract_indiv_score_related_sparse.R
Rscript baselines/scripts/extract_indiv_score_related_dense.R
```

## Notes

- R version: 4.5.0 (2025-04-11)
- Seed used: 600 for all stochastic operations
- `R CMD check --no-manual` completed with 3 WARNINGs and 4 NOTEs (see `baselines/r_cmd_check.log`).

## Additional Extraction (Non-Example Baseline)

- **R repo**: `../STAAR` (local clone of upstream project)
- **Commit**: `b06e478b03f00b7810ea8c75247a1b143f37bda4`
- **Extraction date**: 2026-02-07
- **Purpose**: add distinct non-clone `nonexample601` baseline expansions for `STAAR-47` through `STAAR-53`

### Commands used

Executed in `../STAAR`:

```bash
Rscript /tmp/extract_nonexample601.R
Rscript -e "<inline extraction for related_sparse_binary_spa_nonexample601 and related_sparse_glmmkin_cond_nonexample601 sentinels>"
Rscript -e "<inline extraction for related_dense_glmmkin_cond_nonexample601 and indiv_score_related_sparse_glmmkin_cond_nonexample601 sentinels>"
Rscript -e "<inline extraction for indiv_score_related_dense_glmmkin_cond_nonexample601 sentinels>"
```

Artifacts copied into this repo:

- `baselines/nonexample601_sim_data.rds`
- `baselines/nonexample601_sim_metadata.json`
- `baselines/nonexample601_fingerprint.json`
- `baselines/unrelated_glm_nonexample601_sentinels.json`
- `baselines/related_sparse_glmmkin_nonexample601_sentinels.json`
- `baselines/related_sparse_binary_spa_nonexample601_sentinels.json`
- `baselines/related_sparse_glmmkin_cond_nonexample601_sentinels.json`
- `baselines/related_dense_glmmkin_cond_nonexample601_sentinels.json`
- `baselines/indiv_score_related_sparse_glmmkin_cond_nonexample601_sentinels.json`
- `baselines/indiv_score_related_dense_glmmkin_cond_nonexample601_sentinels.json`
