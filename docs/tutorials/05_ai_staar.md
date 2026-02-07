# 教程 05: AI-STAAR 多祖源分析

AI-STAAR 需要人群分组和分组权重信息。

## 1. 非亲属 AI-STAAR

```python
from pystaar import ai_staar_unrelated_glm

res = ai_staar_unrelated_glm(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["results_STAAR_O"], res["results_ACAT_O"])
```

## 2. 相关样本 AI-STAAR

```python
from pystaar import ai_staar_related_sparse_glmmkin

res = ai_staar_related_sparse_glmmkin(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["results_STAAR_O"])
```

## 3. `find_weight=True` 模式

```python
from pystaar import ai_staar_related_sparse_glmmkin_find_weight

res = ai_staar_related_sparse_glmmkin_find_weight(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["results_weight_staar_o"])
print(res["weight_all_1"])
```

## 4. 自定义人群元数据

当 `dataset != "example"` 时，建议显式传入：

- `pop_groups`
- `pop_weights_1_1`
- `pop_weights_1_25`
- 可选：`pop_levels`

## 5. 输出解释

AI-STAAR 常见输出：

- `results_STAAR_O`
- `results_ACAT_O`
- `results_STAAR_S_1_25`, `results_STAAR_S_1_1`
- `results_STAAR_B_1_25`, `results_STAAR_B_1_1`
- `results_STAAR_A_1_25`, `results_STAAR_A_1_1`

`find_weight=True` 额外输出：

- `weight_all_1`, `weight_all_2`
- `results_weight_staar_o`
- `results_weight1_staar_o`
- `results_weight2_staar_o`
