# 教程 02: 二分类表型与 Binary SPA

当表型是病例-对照（或需要阈值二值化）时，使用 Binary SPA 工作流。

## 1. 非亲属 Binary SPA

```python
from pystaar import staar_unrelated_binary_spa

res = staar_unrelated_binary_spa(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    case_quantile=0.95,
    SPA_p_filter=False,
)

print(res["case_count"], res["results_STAAR_B"])
```

## 2. 相关样本 Binary SPA

```python
from pystaar import staar_related_sparse_binary_spa

res = staar_related_sparse_binary_spa(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    case_quantile=0.95,
    SPA_p_filter=True,
    p_filter_cutoff=0.05,
)

print(res["results_STAAR_B"])
```

## 3. 关键参数

- `case_quantile`: 将连续 `Y` 转成二分类时的阈值分位数。
- `SPA_p_filter`: 是否启用 SPA 预筛选。
- `p_filter_cutoff`: 预筛选阈值（仅当 `SPA_p_filter=True` 生效）。

## 4. 结果字段

- `num_variant`
- `cMAC`
- `case_count`
- `results_STAAR_B`
- `results_STAAR_B_1_25`
- `results_STAAR_B_1_1`

## 5. 实践建议

- 二分类极度不平衡时，优先尝试更稳健的 `case_quantile` 设定并固定 `seed`。
- 若结果不稳定，先在 `example` 或小样本子集上验证流程。
