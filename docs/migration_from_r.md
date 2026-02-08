# 从 R STAAR 迁移到 pySTAAR

本指南面向已经使用 R STAAR 的用户，帮助你快速切换到 Python。

如果你只想先把流程跑通，先看 15 分钟版本：[`migration_r_quickstart_cn.md`](migration_r_quickstart_cn.md)。

## 0. 版本对应关系

- 当前 Python 迁移基线对齐的 R 包版本：`STAAR 0.9.8`
- R 基线来源仓库：`https://github.com/xihaoli/STAAR`
- 基线提取提交：`9db9dd504905b9f469146f670e5f6dbe3e08d01a`

参考：

- `baselines/SOURCE.md`
- `reports/performance.md`

## 1. 核心函数对照

| R 函数 | Python 对应 | 推荐入口 |
|---|---|---|
| `fit_null_glm` | `fit_null_glm` | 底层 API |
| `fit_null_glmmkin` | `fit_null_glmmkin` | 底层 API |
| `STAAR` | `staar` | `staar_*` workflows |
| `STAAR_cond` | `staar_cond` | `staar_*_cond` workflows |
| `STAAR_Binary_SPA` | `staar_binary_spa` | `staar_*_binary_spa` workflows |
| `Indiv_Score_Test_Region` | `indiv_score_test_region` | `indiv_score_*` workflows |
| `Indiv_Score_Test_Region_cond` | `indiv_score_test_region_cond` | `indiv_score_*_cond` workflows |
| `AI_STAAR` | `ai_staar` | `ai_staar_*` workflows |
| `CCT` | `CCT` | 通用聚合工具 |

## 2. 参数对照表（R -> Python）

| R 常见参数 | Python 参数 | 说明 |
|---|---|---|
| `geno` | `dataset`（workflow）或 `genotype`（底层） | workflow 会自动加载 `geno/phred/pheno/kins` |
| `obj_nullmodel` | `obj_nullmodel` | 底层函数一致 |
| `rare_maf_cutoff` | `rare_maf_cutoff` | 一致 |
| `annotation_phred` | `annotation_phred` | 一致 |
| `method_cond` | `method_cond` | 一致 |
| `geno_adj` / 条件位点 | `adj_variant_indices`（workflow）或 `genotype_adj`（底层） | workflow 用索引构造 |
| `find_weight` | `find_weight` 或 `*_find_weight` workflow | 一致 |
| 二分类阈值设定 | `case_quantile` | workflow 先对连续 `Y` 分位数二值化 |
| SPA 预筛选 | `SPA_p_filter`, `p_filter_cutoff` | 一致 |

## 3. 推荐迁移路径

1. 先用 workflow API 替代 R 高层调用。
2. 用 `example` 跑通，再切换自己的数据目录。
3. 需要逐步 debug 时再下沉到底层 API。

## 4. 端到端迁移示例

假设你在 R 中做的是“相关样本 + 条件分析”：

### R 思路

- `fit_null_glmmkin(...)`
- `STAAR_cond(...)`

### Python 对应

```python
from pystaar import staar_related_sparse_glmmkin_cond

res = staar_related_sparse_glmmkin_cond(
    dataset="/path/to/my_dataset_dir",
    seed=600,
    rare_maf_cutoff=0.05,
    method_cond="optimal",
    adj_variant_indices=(0, 3),
)

print(res["results_STAAR_O_cond"])
```

## 5. 与 R 结果差异如何理解

- 当前功能迁移已完成，核心工作流都可用。
- 相关样本路径在纯 Python 计算下存在可重复的小数值漂移。
- 当前 release 口径下，受影响 related parity 哨兵指标已收紧到 `rtol<=3.5e-4`（保留既有 per-sentinel `atol`）。
- `DEV-001` 已关闭，仅作为历史记录保留；当前状态以 `reports/summary.md` 和 `reports/deviations.md` 为准。

## 6. FAQ

### Q1: 为什么我的结果与 R 不完全一致？

A: 常见原因是底层数值后端差异（线性代数实现、优化路径、浮点细节）。如果差异在项目批准容差范围内，属于预期行为。

### Q2: 我必须使用 `use_precomputed_artifacts=True` 吗？

A: 不需要。当前 baseline parity 已在 pure-path（`False`）上运行。该参数仅保留兼容模式，不建议新分析默认开启。

### Q3: 如何判断迁移是否成功？

A: 建议先在同一 `seed`、相同参数下比较关键输出（例如 `results_STAAR_O`、`results_STAAR_B`），再看下游解释是否一致。

### Q4: 为什么条件分析结果和 R 对不上？

A: 最常见原因是条件位点索引。`pySTAAR` 的 `adj_variant_indices` 是 **0-based**，而不少 R 用户习惯按 1-based 思维记录位点顺序。迁移时请先确认索引基准。
