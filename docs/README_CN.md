# pySTAAR 中文快速入门

本项目是 STAAR R 包的 Python 迁移实现，适合希望从 R 工作流迁移到 Python 的用户。

## 1. 你能做什么

`pySTAAR` 当前已覆盖以下核心分析能力：

- STAAR（非亲属 / 亲属 sparse / 亲属 dense）
- 条件分析 `STAAR_cond`
- 二分类表型 `STAAR_Binary_SPA`
- 单变异位点得分检验 `Indiv_Score_Test_Region`
- 多祖源 `AI_STAAR`（含 `find_weight=True`）

## 2. 5 分钟上手

### 2.1 安装

```bash
pip install -e .
```

### 2.2 最小示例（非亲属 STAAR）

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print("num_variant:", res["num_variant"])
print("STAAR-O:", res["results_STAAR_O"])
```

### 2.3 亲属样本示例（稀疏 kinship）

```python
from pystaar import staar_related_sparse_glmmkin

res = staar_related_sparse_glmmkin(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)

print("num_variant:", res["num_variant"])
print("STAAR-O:", res["results_STAAR_O"])
```

## 3. 数据输入方式

`dataset` 支持三种形式：

- 数据集名称（例如 `"example"`）
- 目录路径（目录内必须包含约定文件，见安装文档）
- `STAARDataset` 对象（高级用法）

## 4. 常用工作流入口

- `staar_unrelated_glm`
- `staar_related_sparse_glmmkin`
- `staar_related_dense_glmmkin`
- `staar_unrelated_glm_cond`
- `staar_related_sparse_glmmkin_cond`
- `staar_related_dense_glmmkin_cond`
- `staar_unrelated_binary_spa`
- `staar_related_sparse_binary_spa`
- `staar_related_dense_binary_spa`
- `ai_staar_unrelated_glm`
- `ai_staar_related_sparse_glmmkin`
- `ai_staar_related_dense_glmmkin`

## 5. R 用户迁移建议

- 建议优先从 `workflows.py` 风格的高层 API 开始，而不是直接调用底层统计函数。
- 迁移时先固定 `seed`、`rare_maf_cutoff`、`adj_variant_indices` 等参数，再逐步扩展。
- 先用 `example` 验证流程跑通，再切换到自有数据目录。

详细对照见：[`docs/migration_from_r.md`](migration_from_r.md)

## 6. 与 R 结果的一致性说明

- 当前功能迁移已完成，可用于实际分析。
- 相关样本部分采用已批准的一致性策略：对受影响哨兵指标允许放宽到 `rtol<=5e-4`。
- 该策略记录在 `reports/deviations.md`（`DEV-001`）。

## 7. 文档导航

- 安装与环境：[`installation.md`](installation.md)
- 基础教程：[`tutorials/01_basic_staar.md`](tutorials/01_basic_staar.md)
- Binary SPA：[`tutorials/02_binary_spa.md`](tutorials/02_binary_spa.md)
- 亲属样本：[`tutorials/03_related_samples.md`](tutorials/03_related_samples.md)
- 条件分析：[`tutorials/04_conditional.md`](tutorials/04_conditional.md)
- AI-STAAR：[`tutorials/05_ai_staar.md`](tutorials/05_ai_staar.md)
- Null model API：[`api/null_models.md`](api/null_models.md)
- STAAR API：[`api/staar_functions.md`](api/staar_functions.md)
- 工具函数：[`api/utilities.md`](api/utilities.md)
