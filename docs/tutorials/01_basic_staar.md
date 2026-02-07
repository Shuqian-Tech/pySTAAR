# 教程 01: 基础 STAAR 分析流程

本教程演示最常见的连续性状 STAAR 流程。

## 1. 非亲属样本（GLM）

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print(res["num_variant"])
print(res["results_STAAR_O"])
print(res["results_STAAR_S_1_25"])
```

返回字段主要包括：

- `num_variant`
- `results_STAAR_O`
- `results_STAAR_S_1_25`, `results_STAAR_S_1_1`
- `results_STAAR_B_1_25`, `results_STAAR_B_1_1`
- `results_STAAR_A_1_25`, `results_STAAR_A_1_1`

## 2. 相关样本（GLMM kinship）

```python
from pystaar import staar_related_sparse_glmmkin, staar_related_dense_glmmkin

res_sparse = staar_related_sparse_glmmkin("example", seed=600, rare_maf_cutoff=0.05)
res_dense = staar_related_dense_glmmkin("example", seed=600, rare_maf_cutoff=0.05)

print(res_sparse["results_STAAR_O"], res_dense["results_STAAR_O"])
```

选择建议：

- `sparse`: 更节省内存，适合大规模 kinship。
- `dense`: 小中等规模下实现直接。

## 3. 参数解释

- `seed`: 随机种子，建议在对比实验中固定。
- `rare_maf_cutoff`: 稀有变异筛选阈值，常见值 `0.05` 或 `0.01`。

## 4. 用自定义数据目录运行

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="/path/to/my_dataset_dir",
    seed=600,
    rare_maf_cutoff=0.05,
)
```

目录结构见安装文档：[`../installation.md`](../installation.md)
