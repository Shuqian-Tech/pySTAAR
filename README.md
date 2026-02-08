# pySTAAR

Python 版 STAAR（R 包）迁移项目，面向中文统计遗传/基因组分析用户。

For English docs, see [`docs/README.md`](docs/README.md).

## 项目定位

- 已完成计划内功能迁移（到 `STAAR-56`）。
- 默认 workflow 入口覆盖：STAAR、条件分析、Binary SPA、单变异得分检验、AI-STAAR。
- 当前 parity 基线为 pure-Python 路径（related workflows 不依赖预计算 R 协方差文件）。

## 快速安装

普通用户（发布版）：

```bash
pip install pystaar
```

本仓库开发模式：

```bash
pip install -e '.[dev]'
```

## 快速运行

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="example",
    seed=600,
    rare_maf_cutoff=0.05,
)
print("STAAR-O:", res["results_STAAR_O"])
```

## R 用户迁移入口

- 完整迁移说明：[`docs/migration_from_r.md`](docs/migration_from_r.md)
- 15 分钟迁移清单：[`docs/migration_r_quickstart_cn.md`](docs/migration_r_quickstart_cn.md)
- 数据目录模板：[`docs/data_directory_template_cn.md`](docs/data_directory_template_cn.md)

## 文档导航

- 中文快速入门：[`docs/README_CN.md`](docs/README_CN.md)
- 英文快速入门：[`docs/README.md`](docs/README.md)
- 安装与环境：[`docs/installation.md`](docs/installation.md)
- 性能对比总览（Python vs R）：[`docs/performance_comparison.md`](docs/performance_comparison.md)
- 性能口径说明：官方跨平台结论以 OpenBLAS backend 为准；macOS Accelerate 本地参考见 `examples/1kg_parity/README.md`。
- 本地 1KG 对比示例（数据级 + 模拟完整 workflow）：[`examples/1kg_parity/README.md`](examples/1kg_parity/README.md)
- 教程：
  - [`docs/tutorials/01_basic_staar.md`](docs/tutorials/01_basic_staar.md)
  - [`docs/tutorials/02_binary_spa.md`](docs/tutorials/02_binary_spa.md)
  - [`docs/tutorials/03_related_samples.md`](docs/tutorials/03_related_samples.md)
  - [`docs/tutorials/04_conditional.md`](docs/tutorials/04_conditional.md)
  - [`docs/tutorials/05_ai_staar.md`](docs/tutorials/05_ai_staar.md)
- API 文档：
  - [`docs/api/null_models.md`](docs/api/null_models.md)
  - [`docs/api/staar_functions.md`](docs/api/staar_functions.md)
  - [`docs/api/output_fields.md`](docs/api/output_fields.md)
  - [`docs/api/utilities.md`](docs/api/utilities.md)
  - [`docs/api/stability.md`](docs/api/stability.md)
- 变更记录：[`CHANGELOG.md`](CHANGELOG.md)

## 一致性说明

- 历史偏差记录 `DEV-001` 已关闭，仅保留历史背景。
- 当前状态与 release 口径以 [`reports/summary.md`](reports/summary.md) 和 [`reports/deviations.md`](reports/deviations.md) 为准。
