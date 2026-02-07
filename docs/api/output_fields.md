# API: 输出字段参考

本页汇总各 workflow 返回 dict 的主要键名与含义，便于下游管道稳定对接。

## 1. STAAR 主流程

适用函数：

- `staar_unrelated_glm`
- `staar_related_sparse_glmmkin`
- `staar_related_dense_glmmkin`

| 键名 | 含义 |
|---|---|
| `num_variant` | 当前集合中通过稀有位点过滤后的位点数 |
| `results_STAAR_O` | STAAR-O 聚合 p-value |
| `results_STAAR_S_1_25` | STAAR-S(1,25) 映射结果（包含多个测试项） |
| `results_STAAR_S_1_1` | STAAR-S(1,1) 映射结果 |
| `results_STAAR_B_1_25` | STAAR-B(1,25) 映射结果 |
| `results_STAAR_B_1_1` | STAAR-B(1,1) 映射结果 |
| `results_STAAR_A_1_25` | STAAR-A(1,25) 映射结果 |
| `results_STAAR_A_1_1` | STAAR-A(1,1) 映射结果 |

额外键：

- `staar_unrelated_glm` 还返回：`nullmodel_beta`, `nullmodel_dispersion`
- related 版本还返回：`nullmodel_theta`（`dispersion`, `kins1`）

## 2. 条件分析流程

适用函数：

- `staar_unrelated_glm_cond`
- `staar_related_sparse_glmmkin_cond`
- `staar_related_dense_glmmkin_cond`

| 键名 | 含义 |
|---|---|
| `num_variant` | 稀有位点数量 |
| `cMAC` | 累积 minor allele count |
| `results_STAAR_O_cond` | 条件 STAAR-O |
| `results_ACAT_O_cond` | 条件 ACAT-O |
| `results_STAAR_S_1_25_cond` | 条件 STAAR-S(1,25) 映射结果 |
| `results_STAAR_S_1_1_cond` | 条件 STAAR-S(1,1) 映射结果 |
| `results_STAAR_B_1_25_cond` | 条件 STAAR-B(1,25) 映射结果 |
| `results_STAAR_B_1_1_cond` | 条件 STAAR-B(1,1) 映射结果 |
| `results_STAAR_A_1_25_cond` | 条件 STAAR-A(1,25) 映射结果 |
| `results_STAAR_A_1_1_cond` | 条件 STAAR-A(1,1) 映射结果 |

## 3. Binary SPA 流程

适用函数：

- `staar_unrelated_binary_spa`
- `staar_related_sparse_binary_spa`
- `staar_related_dense_binary_spa`

| 键名 | 含义 |
|---|---|
| `num_variant` | 稀有位点数量 |
| `cMAC` | 累积 minor allele count |
| `case_count` | 二值化后病例数量 |
| `results_STAAR_B` | STAAR-B 主统计量 p-value |
| `results_STAAR_B_1_25` | STAAR-B(1,25) 映射结果 |
| `results_STAAR_B_1_1` | STAAR-B(1,1) 映射结果 |

## 4. 单变异位点得分检验流程

适用函数：

- `indiv_score_unrelated_glm`
- `indiv_score_related_sparse_glmmkin`
- `indiv_score_related_dense_glmmkin`
- `indiv_score_unrelated_glm_cond`
- `indiv_score_related_sparse_glmmkin_cond`
- `indiv_score_related_dense_glmmkin_cond`

公共键：

| 键名 | 含义 |
|---|---|
| `num_variant` | 稀有位点数量 |
| `num_tested` | 实际参与检验的位点数量 |
| `pvalue_min` | 最小 p-value |
| `pvalue_median` | p-value 中位数 |
| `top_variant_index` | 最显著位点索引（1-based） |
| `top_pvalue` | 最显著位点 p-value |

非条件版额外映射：

- `score_samples`
- `se_samples`
- `pvalue_samples`

条件版额外映射：

- `score_cond_samples`
- `se_cond_samples`
- `pvalue_cond_samples`

说明：上述 `*_samples` 为抽样位点映射（例如 `v2`, `v3`），用于快速检查而非完整向量输出。

## 5. AI-STAAR 流程

适用函数：

- `ai_staar_unrelated_glm`
- `ai_staar_related_sparse_glmmkin`
- `ai_staar_related_dense_glmmkin`

| 键名 | 含义 |
|---|---|
| `num_variant` | 稀有位点数量 |
| `cMAC` | 累积 minor allele count |
| `results_STAAR_O` | AI-STAAR 的 STAAR-O |
| `results_ACAT_O` | AI-STAAR 的 ACAT-O |
| `results_STAAR_S_1_25` | AI-STAAR S(1,25) 映射结果 |
| `results_STAAR_S_1_1` | AI-STAAR S(1,1) 映射结果 |
| `results_STAAR_B_1_25` | AI-STAAR B(1,25) 映射结果 |
| `results_STAAR_B_1_1` | AI-STAAR B(1,1) 映射结果 |
| `results_STAAR_A_1_25` | AI-STAAR A(1,25) 映射结果 |
| `results_STAAR_A_1_1` | AI-STAAR A(1,1) 映射结果 |

## 6. AI-STAAR `find_weight=True`

适用函数：

- `ai_staar_unrelated_glm_find_weight`
- `ai_staar_related_sparse_glmmkin_find_weight`
- `ai_staar_related_dense_glmmkin_find_weight`

除主流程字段外，还会返回：

| 键名 | 含义 |
|---|---|
| `weight_all_1` | 各人群在第一组 base test 的权重映射 |
| `weight_all_2` | 各人群在第二组 base test 的权重映射 |
| `results_weight_staar_o` | 按权重组合后的 STAAR-O 映射 |
| `results_weight1_staar_o` | 第一组权重下 STAAR-O 映射 |
| `results_weight2_staar_o` | 第二组权重下 STAAR-O 映射 |
| `results_weight_staar_s_1_25` | 权重组合后的 STAAR-S(1,25) 映射 |
| `results_weight1_staar_s_1_25` | 第一组权重下 STAAR-S(1,25) 映射 |
| `results_weight2_staar_s_1_25` | 第二组权重下 STAAR-S(1,25) 映射 |

## 7. 稳定性建议

- 下游系统请优先依赖公共核心键，不要假设不同 workflow 的返回键完全相同。
- 若做长期存储，建议固定版本并保存 `seed`、`rare_maf_cutoff`、workflow 名称和关键参数。
