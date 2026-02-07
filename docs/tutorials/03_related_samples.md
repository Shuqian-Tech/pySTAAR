# 教程 03: 亲属样本分析

亲属样本分析通常使用 GLMM+kinship 工作流。

## 1. sparse 与 dense 的选择

- `staar_related_sparse_glmmkin`: 大规模数据优先。
- `staar_related_dense_glmmkin`: 数据量较小时更直接。

## 2. 基础调用

```python
from pystaar import staar_related_sparse_glmmkin

res = staar_related_sparse_glmmkin(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["nullmodel_theta"])
```

## 3. 单变异位点检验

```python
from pystaar import indiv_score_related_sparse_glmmkin

res = indiv_score_related_sparse_glmmkin(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["pvalue_min"], res["top_variant_index"])
```

## 4. `use_precomputed_artifacts` 参数说明

- 当前推荐默认值：`False`。
- baseline parity 也已切换为 pure-path（`False`）。
- 该参数保留为兼容模式，不建议普通用户在新分析中启用。

## 5. 注意事项

- 亲属样本结果与 R 结果通常高度一致，但会存在小幅数值漂移。
- 当前项目已批准受影响哨兵指标容差上限 `rtol=5e-4`（见 `reports/deviations.md`）。
