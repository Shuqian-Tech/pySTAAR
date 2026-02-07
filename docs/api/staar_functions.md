# API: STAAR / 条件分析 / AI-STAAR

本页整理核心分析函数与推荐入口。

## 1. 推荐入口（workflow API）

这些函数最适合端到端分析：

### 基础 STAAR

- `staar_unrelated_glm`
- `staar_related_sparse_glmmkin`
- `staar_related_dense_glmmkin`

### 条件分析

- `staar_unrelated_glm_cond`
- `staar_related_sparse_glmmkin_cond`
- `staar_related_dense_glmmkin_cond`

### Binary SPA

- `staar_unrelated_binary_spa`
- `staar_related_sparse_binary_spa`
- `staar_related_dense_binary_spa`

### 单变异得分检验

- `indiv_score_unrelated_glm`
- `indiv_score_related_sparse_glmmkin`
- `indiv_score_related_dense_glmmkin`
- `indiv_score_unrelated_glm_cond`
- `indiv_score_related_sparse_glmmkin_cond`
- `indiv_score_related_dense_glmmkin_cond`

### AI-STAAR

- `ai_staar_unrelated_glm`
- `ai_staar_related_sparse_glmmkin`
- `ai_staar_related_dense_glmmkin`
- `ai_staar_unrelated_glm_find_weight`
- `ai_staar_related_sparse_glmmkin_find_weight`
- `ai_staar_related_dense_glmmkin_find_weight`

## 2. 底层函数（高级用户）

- `staar`
- `staar_cond`
- `staar_binary_spa`
- `ai_staar`
- `indiv_score_test_region`
- `indiv_score_test_region_cond`

这些函数要求你手动构造 `obj_nullmodel`，适合做细粒度方法实验。

## 3. 典型参数

大多数 workflow 共享参数：

- `dataset`: `"example"` / 自定义名称 / 数据目录 / `STAARDataset`
- `seed`: 随机种子
- `rare_maf_cutoff`: 稀有变异阈值

条件分析额外参数：

- `method_cond`
- `adj_variant_indices`

Binary SPA 额外参数：

- `case_quantile`
- `SPA_p_filter`
- `p_filter_cutoff`

AI-STAAR 额外参数：

- `pop_groups`
- `pop_weights_1_1`
- `pop_weights_1_25`
- `pop_levels`

## 4. R 兼容函数名

`pystaar` 提供 R 风格别名：

- `STAAR` -> `staar`
- `STAAR_cond` -> `staar_cond`
- `STAAR_Binary_SPA` -> `staar_binary_spa`
- `AI_STAAR` -> `ai_staar`
- `Indiv_Score_Test_Region` -> `indiv_score_test_region`
- `Indiv_Score_Test_Region_cond` -> `indiv_score_test_region_cond`

## 5. 输出结构建议

不同 workflow 返回 dict，键名与 parity spec 保持一致。实践中建议只依赖你真正使用的字段，不要假设所有 workflow 返回同一组键。

## 6. 输出字段参考表

完整键名与含义对照见：[`output_fields.md`](output_fields.md)。
