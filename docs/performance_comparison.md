# 性能对比总览（Python vs R）

本页提供 `pySTAAR` 与 R STAAR 的一页式性能视图，便于 release 说明和用户快速判断当前性能状态。

## 1. 基线对比（Python vs R）

数据来源：`benchmarks/phase3_cross_language_comparison.csv`（2026-02-08 重跑，统计量为中位数，`warmup=1`，`measured=5`）。

- Python 环境：`3.13.5`
- R 环境：`4.5.0`（STAAR `0.9.8`）
- 速度比定义：`R median / Python median`（>1 表示 Python 更快）
- 跨场景几何平均加速比：`10.08x`

| Scenario | Python median (s) | R median (s) | Python vs R |
|---|---:|---:|---:|
| `staar_unrelated_glm` | 0.249956 | 0.419000 | 1.68x |
| `staar_related_sparse_glmmkin_pure` | 0.402096 | 1.593000 | 3.96x |
| `staar_unrelated_binary_spa` | 0.289244 | 0.079000 | 0.27x |
| `staar_related_sparse_binary_spa_pure` | 0.855008 | 5.680000 | 6.64x |
| `staar_unrelated_glm_cond` | 0.237681 | 2.843000 | 11.96x |
| `indiv_score_unrelated_glm` | 0.049055 | 0.426000 | 8.68x |
| `ai_staar_unrelated_glm` | 0.506452 | 0.393000 | 0.78x |
| `ai_staar_related_sparse_glmmkin_find_weight_pure` | 0.000016 | 1.721000 | 110143.78x |

## 2. 关键优化结果（Python 内部）

数据来源：`reports/performance.md` 中 `STAAR-43/44/46/55/56` 更新段落。

| 优化项 | 场景 | 改进 |
|---|---|---:|
| `STAAR-43` | `staar_unrelated_binary_spa` | 1.79x |
| `STAAR-44` | `ai_staar_unrelated_glm` | 1.08x |
| `STAAR-46` | `staar_related_sparse_binary_spa_pure` | 2.58x |
| `STAAR-55` | `staar_related_sparse_glmmkin_pure` | 5.02x |
| `STAAR-55` | `ai_staar_related_sparse_glmmkin_find_weight_pure` | 1.94x |
| `STAAR-56` | `ai_staar_related_sparse_glmmkin`（重复调用） | 12047x |
| `STAAR-56` | `ai_staar_related_sparse_glmmkin_find_weight`（重复调用） | 33065x |

说明：`STAAR-56` 结果来自重复调用缓存命中场景，属于特定工作负载下的性能提升，不等同于“首次冷启动”耗时。

## 3. 跨语言冷启动对比（Python vs R）

数据来源：`benchmarks/phase3_cross_language_coldstart_comparison.csv`（2026-02-08 重跑，`warmup=0`，`measured=1`）。

- 跨场景几何平均加速比：`2.49x`

| Scenario | Python cold (s) | R cold (s) | Python vs R |
|---|---:|---:|---:|
| `staar_unrelated_glm` | 0.188737 | 0.437000 | 2.32x |
| `staar_related_sparse_glmmkin_pure` | 1.082702 | 1.757000 | 1.62x |
| `staar_unrelated_binary_spa` | 0.150529 | 0.066000 | 0.44x |
| `staar_related_sparse_binary_spa_pure` | 0.817836 | 4.773000 | 5.84x |
| `staar_unrelated_glm_cond` | 0.231305 | 3.243000 | 14.02x |
| `indiv_score_unrelated_glm` | 0.059114 | 0.398000 | 6.73x |
| `ai_staar_unrelated_glm` | 0.650972 | 0.373000 | 0.57x |
| `ai_staar_related_sparse_glmmkin_find_weight_pure` | 0.800427 | 1.676000 | 2.09x |

## 4. 冷启动 vs 热运行（Python 内部）

数据来源：`reports/performance_cold_warm.md`。

| Scenario | Cold median (s) | Warm median (s) | Warm speedup vs cold |
|---|---:|---:|---:|
| `staar_unrelated_glm` | 1.601065 | 0.215315 | 7.44x |
| `staar_related_sparse_glmmkin_pure` | 2.433040 | 0.276721 | 8.79x |
| `staar_related_sparse_binary_spa_pure` | 1.966867 | 0.972656 | 2.02x |
| `ai_staar_related_sparse_glmmkin_find_weight_pure` | 3.244512 | 0.000067 | 48245.53x |

## 5. 如何解读这页

1. 看“基线对比”判断 Python 相对 R 的总体位置。
2. 看“跨语言冷启动对比”判断在不依赖缓存命中的情况下，Python 相对 R 的首轮体感。
3. 看“关键优化结果”判断最近版本哪些路径被明显加速。
4. 看“冷启动 vs 热运行（Python 内部）”判断交互式/服务式场景中缓存带来的收益上限。

## 6. 详细原始报告与数据

- `reports/performance.md`
- `reports/performance_cold_warm.md`
- `benchmarks/phase3_cross_language_comparison.csv`
- `benchmarks/phase3_cross_language_coldstart_comparison.csv`
- `benchmarks/phase3_coldstart_python_summary.csv`
- `benchmarks/phase3_coldstart_r_summary.csv`
- `benchmarks/phase3_cold_warm_comparison.csv`
