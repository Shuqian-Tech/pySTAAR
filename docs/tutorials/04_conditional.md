# 教程 04: 条件分析 (STAAR_cond)

条件分析用于在给定已知位点后评估目标区域信号。

## 1. 非亲属条件分析

```python
from pystaar import staar_unrelated_glm_cond

res = staar_unrelated_glm_cond(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    method_cond="optimal",
    adj_variant_indices=(0, 3),
)

print(res["results_STAAR_O_cond"])
```

## 2. 相关样本条件分析

```python
from pystaar import staar_related_sparse_glmmkin_cond

res = staar_related_sparse_glmmkin_cond(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    method_cond="optimal",
    adj_variant_indices=(0, 3),
)

print(res["results_STAAR_O_cond"])
```

## 3. 参数说明

- `method_cond`: 条件设计矩阵构造策略，默认 `"optimal"`。
- `adj_variant_indices`: 作为条件协变量的变异位点索引（0-based）。

## 4. 错误排查

若出现范围错误：

- 检查 `adj_variant_indices` 是否为空。
- 检查索引是否超出 `geno` 的列范围。

## 5. 单变异条件检验

对应入口：

- `indiv_score_unrelated_glm_cond`
- `indiv_score_related_sparse_glmmkin_cond`
- `indiv_score_related_dense_glmmkin_cond`
