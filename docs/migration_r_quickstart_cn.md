# R 用户 15 分钟迁移清单

本页给已经熟悉 R STAAR 的用户一个最短迁移路径，目标是快速把现有脚本切到 `pySTAAR` workflow API。

## 1. 先确认三件事

1. 先固定参数：`seed`、`rare_maf_cutoff`、`adj_variant_indices`。
2. 先在 `dataset="example"` 跑通，再切换到自有数据目录。
3. 迁移时优先使用 workflow API，不要一开始就下沉到底层 null model。

## 2. 常用 R -> Python 对照

| R 常见调用 | Python 对应入口 |
|---|---|
| `STAAR(...)` | `staar_unrelated_glm` / `staar_related_sparse_glmmkin` / `staar_related_dense_glmmkin` |
| `STAAR_cond(...)` | `staar_unrelated_glm_cond` / `staar_related_sparse_glmmkin_cond` / `staar_related_dense_glmmkin_cond` |
| `STAAR_Binary_SPA(...)` | `staar_unrelated_binary_spa` / `staar_related_sparse_binary_spa` / `staar_related_dense_binary_spa` |
| `AI_STAAR(...)` | `ai_staar_unrelated_glm` / `ai_staar_related_sparse_glmmkin` / `ai_staar_related_dense_glmmkin` |
| `AI_STAAR(..., find_weight=TRUE)` | 对应 `*_find_weight` workflow |

## 3. 最小迁移示例

### 3.1 基础 STAAR（非亲属）

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)
print(res["results_STAAR_O"])
```

### 3.2 相关样本 + 条件分析

```python
from pystaar import staar_related_sparse_glmmkin_cond

res = staar_related_sparse_glmmkin_cond(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    method_cond="optimal",
    adj_variant_indices=(0, 3),  # 注意：0-based
)
print(res["results_STAAR_O_cond"])
```

### 3.3 Binary SPA

```python
from pystaar import staar_unrelated_binary_spa

res = staar_unrelated_binary_spa(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
    case_quantile=0.95,
    SPA_p_filter=True,
    p_filter_cutoff=0.05,
)
print(res["results_STAAR_B"])
```

## 4. R 用户最容易踩的 4 个坑

1. `adj_variant_indices` 是 0-based，不是 1-based。
2. `dataset` 传目录时，必须包含 6 个固定文件名。
3. Binary SPA workflow 会先按 `case_quantile` 对 `Y` 分位数二值化。
4. `use_precomputed_artifacts` 仅是兼容模式；baseline parity 默认走 pure-path（`False`）。

## 5. 迁移完成判据

1. 同一 `seed`、同一参数下，关键哨兵字段（如 `results_STAAR_O` / `results_STAAR_B`）与预期一致。
2. 跑通本仓库 parity 套件：`pytest tests/parity -q`。
3. 结果解释与下游统计决策保持一致。

## 6. 进一步阅读

- 详细 R 对照：[`migration_from_r.md`](migration_from_r.md)
- 安装与环境：[`installation.md`](installation.md)
- 数据目录模板：[`data_directory_template_cn.md`](data_directory_template_cn.md)
- 输出字段说明：[`api/output_fields.md`](api/output_fields.md)
