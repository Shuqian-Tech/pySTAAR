# 1KG Local Parity Examples

本目录提供两类本地示例：

1. 数据级 parity（`genotype/snploc/maf`）  
2. 基于 1KG 基因型的可复现“模拟完整 workflow”对比（R vs Python）

## Backend Note

- `pySTAAR` 当前包含 backend-aware eigensolver 自动选择逻辑：
  - OpenBLAS 环境优先 SciPy 特征值路径（`scipy.linalg.eigh`）。
  - Accelerate 环境优先 NumPy 特征值路径（`np.linalg.eigvalsh`）。
- 示例报告 JSON 的 `meta` 会记录 `python_executable`、`numpy/scipy` 版本和 `blas_backend_hint`，便于审计。
- 官方跨平台性能口径仍以 OpenBLAS 为准；本目录结果用于本地 workflow 验证与性能观察。

## Prerequisites

- Python 环境可运行本仓库（推荐已执行 `pip install -e .[dev]`）。
- R 环境已安装 `STAAR`、`Matrix`、`jsonlite`（以及 `STAAR` 依赖）。
- 本地有输入文件：`~/Downloads/kgp3202_chr1_example_like.rda`

## 1) 数据级 parity + 性能

运行：

```bash
bash examples/1kg_parity/run_local.sh
```

也可显式指定解释器和报告标签（防止覆盖）：

```bash
PYTHON_BIN=/Library/Developer/CommandLineTools/usr/bin/python3 \
REPORT_TAG=mac_accelerate \
bash examples/1kg_parity/run_local.sh
```

输出：

- `examples/1kg_parity/reports/parity_and_perf.md`
- `examples/1kg_parity/reports/parity_and_perf.json`

说明：该模式只比较数据层，不运行完整 STAAR workflow。

## 2) 模拟完整 workflow（推荐发布展示）

运行：

```bash
bash examples/1kg_parity/run_full_workflow_local.sh
```

推荐显式指定解释器和报告标签：

```bash
PYTHON_BIN=/Library/Developer/CommandLineTools/usr/bin/python3 \
REPORT_TAG=mac_accelerate \
bash examples/1kg_parity/run_full_workflow_local.sh
```

该命令会：

1. 从 `.rda` 读取 `genotype/snploc/maf`。
2. 可复现地模拟 `phred/pheno/kins` 并生成本地数据目录。
3. 在同一模拟数据上分别运行 R 与 Python 的 workflow 套件并对比。
4. 输出数值一致性与性能对比报告。

输出：

- `examples/1kg_parity/reports/full_workflow_parity_and_perf.md`
- `examples/1kg_parity/reports/full_workflow_parity_and_perf.json`

## Full Workflow Parameters (env vars)

- `TARGET_NSAMPLE`：模拟数据样本数（默认 `1200`）
- `TARGET_NVAR`：模拟数据位点数（默认 `1500`）
- `SAMPLE_MAF_MAX`：抽样位点的 MAF 上限（默认 `0.2`）
- `SEED`：全流程随机种子（默认 `600`）
- `RARE_MAF_CUTOFF`：workflow 稀有变异阈值（默认 `0.05`）
- `ADJ_VARIANTS`：条件分析位点索引（1-based，默认 `1`）
- `RUNS` / `WARMUP`：性能评测次数（默认 `3` / `1`）
- `ATOL` / `RTOL`：R vs Python 数值比较阈值（默认 `1e-6` / `1e-3`）
- `PYTHON_BIN`：Python 解释器（默认 `python3`）
- `REPORT_TAG`：输出报告后缀标签（默认空）

示例：

```bash
TARGET_NSAMPLE=1600 TARGET_NVAR=2500 RUNS=3 WARMUP=1 \
bash examples/1kg_parity/run_full_workflow_local.sh
```

## Backend Comparison (sample on same machine)

建议固定参数：

```bash
TARGET_NSAMPLE=1200 TARGET_NVAR=1500 RUNS=3 WARMUP=1 SEED=600
```

本机示例（macOS，backend-aware 自动选择）：

| Backend | Python | Python total median (s) | R/Python speedup | Note |
|---|---|---:|---:|---|
| Accelerate (sample) | `3.9.6` | 10.6177 | 2.282x | `PYTHON_BIN=/Library/Developer/CommandLineTools/usr/bin/python3` |

说明：

- 上表是同机样例，数字会随硬件、线程设置、NumPy/SciPy 版本变化。
- 若要复现官方跨平台口径，请参考 OpenBLAS 基线（`docs/performance_comparison.md`）。
- 样例来源：`reports/full_workflow_parity_and_perf_mac_accelerate.json`（`RUNS=3`, `WARMUP=1`）。

## Notes

- `data/` 与 `reports/` 在本目录下默认本地忽略，不会提交到 Git。
- “模拟完整 workflow”中的 `phred/pheno/kins` 是可复现模拟数据，用于对比流程与性能，不代表真实临床/生产数据分布。
