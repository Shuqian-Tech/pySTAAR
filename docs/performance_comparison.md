# 性能对比总览（Python vs R）

本页提供 `pySTAAR` 与 R STAAR 的一页式性能视图，便于 release 说明和用户快速判断当前性能状态。

## 1. 基线对比（Python vs R）

数据来源：`reports/performance.md`（Phase 3 baseline，统计量为中位数，5 次测量，1 次 warm-up 丢弃）。

- Python 环境：`3.13.5`
- R 环境：`4.5.0`（STAAR `0.9.8`）
- 速度比定义：`R median / Python median`（>1 表示 Python 更快）
- 跨场景几何平均加速比：`1.64x`

| Scenario | Python median (s) | R median (s) | Python vs R |
|---|---:|---:|---:|
| `staar_unrelated_glm` | 0.235836 | 0.437000 | 1.85x |
| `staar_related_sparse_glmmkin_pure` | 1.004456 | 1.625000 | 1.62x |
| `staar_unrelated_binary_spa` | 0.288611 | 0.085000 | 0.29x |
| `staar_related_sparse_binary_spa_pure` | 0.715168 | 4.725000 | 6.61x |
| `staar_unrelated_glm_cond` | 0.302497 | 2.136000 | 7.06x |
| `indiv_score_unrelated_glm` | 0.120387 | 0.382000 | 3.17x |
| `ai_staar_unrelated_glm` | 0.928488 | 0.384000 | 0.41x |
| `ai_staar_related_sparse_glmmkin_find_weight_pure` | 1.664256 | 1.624000 | 0.98x |

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

## 3. 冷启动 vs 热运行（用户体验视角）

数据来源：`reports/performance_cold_warm.md`。

| Scenario | Cold median (s) | Warm median (s) | Warm speedup vs cold |
|---|---:|---:|---:|
| `staar_unrelated_glm` | 1.601065 | 0.215315 | 7.44x |
| `staar_related_sparse_glmmkin_pure` | 2.433040 | 0.276721 | 8.79x |
| `staar_related_sparse_binary_spa_pure` | 1.966867 | 0.972656 | 2.02x |
| `ai_staar_related_sparse_glmmkin_find_weight_pure` | 3.244512 | 0.000067 | 48245.53x |

## 4. 如何解读这页

1. 看“基线对比”判断 Python 相对 R 的总体位置。
2. 看“关键优化结果”判断最近版本哪些路径被明显加速。
3. 看“冷启动 vs 热运行”判断交互式/服务式场景的真实体感。

## 5. 详细原始报告与数据

- `reports/performance.md`
- `reports/performance_cold_warm.md`
- `benchmarks/phase3_cross_language_comparison.csv`
- `benchmarks/phase3_cold_warm_comparison.csv`
