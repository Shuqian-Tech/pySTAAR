## Workflows to extract

1. **staar_unrelated_glm** — `docs/STAAR_vignette.html` ("Analyzing the simulated data using STAAR" → "Unrelated samples")
- Inputs: STAAR package dataset `example` (genotype/maf/snploc), simulated covariates/annotations
- Key sentinels: `results_STAAR_O`, conventional test p-values, STAAR-S/B/A tables (first row), `num_variant`, null model coefficients

2. **staar_related_sparse_glmmkin** — `docs/STAAR_vignette.html` ("Related samples" → "Analyzing data using sparse GRM")
- Inputs: STAAR `example` dataset, simulated covariates/annotations, sparse kinship matrix
- Key sentinels: `results_STAAR_O`, conventional test p-values, STAAR-S/B/A tables (first row), `num_variant`, null model coefficients/variance

3. **staar_related_dense_glmmkin** — `docs/STAAR_vignette.html` ("Related samples" → "Analyzing data using dense GRM")
- Inputs: STAAR `example` dataset, simulated covariates/annotations, dense kinship matrix
- Key sentinels: `results_STAAR_O`, conventional test p-values, STAAR-S/B/A tables (first row), `num_variant`, null model coefficients/variance

4. **staar_unrelated_glm_cond** — conditional omnibus test using `STAAR_cond` with unrelated null model
- Inputs: STAAR `example` dataset, simulated covariates/annotations, conditioning variants from example genotype
- Key sentinels: `results_STAAR_O_cond`, `results_ACAT_O_cond`, STAAR-S/B/A conditional tables (first row), `num_variant`, `cMAC`

5. **staar_related_sparse_glmmkin_cond** — conditional omnibus test using `STAAR_cond` with sparse related null model
- Inputs: STAAR `example` dataset, simulated covariates/annotations, sparse kinship matrix, conditioning variants from example genotype
- Key sentinels: `results_STAAR_O_cond`, `results_ACAT_O_cond`, STAAR-S/B/A conditional tables (first row), `num_variant`, `cMAC`

6. **staar_related_dense_glmmkin_cond** — conditional omnibus test using `STAAR_cond` with dense related null model
- Inputs: STAAR `example` dataset, simulated covariates/annotations, dense kinship matrix, conditioning variants from example genotype
- Key sentinels: `results_STAAR_O_cond`, `results_ACAT_O_cond`, STAAR-S/B/A conditional tables (first row), `num_variant`, `cMAC`

7. **staar_unrelated_binary_spa** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA`
- Inputs: STAAR `example` dataset, binary phenotype derived from unrelated example covariates (`Y > quantile(Y, 0.95)`), functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

8. **staar_related_sparse_binary_spa** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA` with sparse related null model
- Inputs: STAAR `example` dataset, binary phenotype derived from related covariates (`Y > quantile(Y, 0.95)`), sparse kinship matrix, functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

9. **staar_related_dense_binary_spa** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA` with dense related null model
- Inputs: STAAR `example` dataset, binary phenotype derived from related covariates (`Y > quantile(Y, 0.95)`), dense kinship matrix, functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

10. **staar_unrelated_binary_spa_filter** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA` with `SPA_p_filter=TRUE` for unrelated samples
- Inputs: STAAR `example` dataset, binary phenotype derived from unrelated covariates (`Y > quantile(Y, 0.95)`), functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

11. **staar_related_sparse_binary_spa_filter** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA` with `SPA_p_filter=TRUE` and sparse related null model
- Inputs: STAAR `example` dataset, binary phenotype derived from related covariates (`Y > quantile(Y, 0.95)`), sparse kinship matrix, functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

12. **staar_related_dense_binary_spa_filter** — imbalanced case-control omnibus burden test using `STAAR_Binary_SPA` with `SPA_p_filter=TRUE` and dense related null model
- Inputs: STAAR `example` dataset, binary phenotype derived from related covariates (`Y > quantile(Y, 0.95)`), dense kinship matrix, functional annotations
- Key sentinels: `results_STAAR_B`, STAAR-B(1,25)/(1,1) tables (first row), `num_variant`, `cMAC`, `case_count`

13. **indiv_score_unrelated_glm** — individual-variant score test using `Indiv_Score_Test_Region` with unrelated null model
- Inputs: STAAR `example` dataset, unrelated phenotype/covariates
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled Score/SE/p-values

14. **indiv_score_related_sparse_glmmkin** — individual-variant score test using `Indiv_Score_Test_Region` with sparse related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, sparse kinship matrix
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled Score/SE/p-values

15. **indiv_score_related_dense_glmmkin** — individual-variant score test using `Indiv_Score_Test_Region` with dense related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, dense kinship matrix
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled Score/SE/p-values

16. **indiv_score_unrelated_glm_cond** — conditional individual-variant score test using `Indiv_Score_Test_Region_cond` with unrelated null model
- Inputs: STAAR `example` dataset, unrelated phenotype/covariates, conditioning variants from example genotype
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled conditional Score/SE/p-values

17. **indiv_score_related_sparse_glmmkin_cond** — conditional individual-variant score test using `Indiv_Score_Test_Region_cond` with sparse related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, sparse kinship matrix, conditioning variants from example genotype
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled conditional Score/SE/p-values

18. **indiv_score_related_dense_glmmkin_cond** — conditional individual-variant score test using `Indiv_Score_Test_Region_cond` with dense related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, dense kinship matrix, conditioning variants from example genotype
- Key sentinels: `num_variant`, `num_tested`, min/median p-value, top-hit index/p-value, sampled conditional Score/SE/p-values

19. **ai_staar_unrelated_glm** — ancestry-informed STAAR omnibus test using `AI_STAAR` with unrelated null model
- Inputs: STAAR `example` dataset, unrelated phenotype/covariates, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, STAAR-S/B/A tables (first row)

20. **ai_staar_related_sparse_glmmkin** — ancestry-informed STAAR omnibus test using `AI_STAAR` with sparse related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, sparse kinship matrix, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, STAAR-S/B/A tables (first row)

21. **ai_staar_related_dense_glmmkin** — ancestry-informed STAAR omnibus test using `AI_STAAR` with dense related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, dense kinship matrix, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, STAAR-S/B/A tables (first row)

22. **ai_staar_unrelated_glm_find_weight** — ancestry-informed STAAR omnibus test using `AI_STAAR(..., find_weight=TRUE)` with unrelated null model
- Inputs: STAAR `example` dataset, unrelated phenotype/covariates, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, ancestry weight matrices, per-base weighted STAAR metrics

23. **ai_staar_related_sparse_glmmkin_find_weight** — ancestry-informed STAAR omnibus test using `AI_STAAR(..., find_weight=TRUE)` with sparse related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, sparse kinship matrix, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, ancestry weight matrices, per-base weighted STAAR metrics

24. **ai_staar_related_dense_glmmkin_find_weight** — ancestry-informed STAAR omnibus test using `AI_STAAR(..., find_weight=TRUE)` with dense related null model
- Inputs: STAAR `example` dataset, related phenotype/covariates, dense kinship matrix, deterministic ancestry groups and base-test weight matrices
- Key sentinels: `num_variant`, `cMAC`, `results_STAAR_O`, `results_ACAT_O`, ancestry weight matrices, per-base weighted STAAR metrics
