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
- 当前 release 口径下，受影响 related parity 哨兵指标使用收紧后的容差上限 `rtol=3.5e-4`（见 `reports/summary.md` 与 `reports/deviations.md`）。

## 6. 性能建议（sparse vs dense）

### 6.1 选择策略

- 样本规模较大、kinship 稀疏：优先 `sparse`。
- 样本规模中等且内存充足：可用 `dense` 对照验证。

### 6.2 内存与速度经验

- `dense` 通常内存占用更高，矩阵规模增长时压力明显。
- `sparse` 在大规模 kinship 下更可扩展，但不同后端下耗时表现可能波动。

### 6.3 实操建议

- 先在小子集上做功能验证，再放大全量数据。
- 固定 `seed` 后再比较速度，避免随机差异干扰。
- 同机对比时保持相同线程和依赖环境。

更多性能基线见：`reports/performance.md`。
