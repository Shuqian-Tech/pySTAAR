# 安装与环境配置

## 1. 环境要求

- Python `>=3.9`
- 推荐使用虚拟环境（`venv`/`conda`）

## 2. 本地开发安装

在仓库根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

安装开发测试依赖：

```bash
pip install -e '.[dev]'
```

## 3. 发布版用户安装（非开发模式）

如果你是普通使用者（不是本仓库开发者），推荐：

```bash
pip install pystaar
```

或安装指定版本：

```bash
pip install pystaar==1.0.0
```

## 4. 依赖版本

项目依赖在 `pyproject.toml` 中定义，核心依赖包括：

- `numpy>=1.24`
- `scipy>=1.10`
- `pandas>=2.0`
- `pyyaml>=6.0`

## 5. 性能复现环境（重要）

- `pip install pystaar` 只保证功能可用，不保证不同机器的性能数字一致。
- 若要复现官方跨平台性能口径，建议使用 OpenBLAS backend（见 `docs/performance_comparison.md`）。
- macOS 用户可使用 Apple Accelerate backend 作为本地性能参考（在部分场景会更快）。

后端自检命令：

```bash
python - <<'PY'
import numpy as np
np.__config__.show()
PY
```

建议在性能测试前固定线程环境变量（示例）：

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

可选：手动覆盖 eigensolver 选择策略（默认自动）：

```bash
export PYSTAAR_EIGENSOLVER=auto   # auto / numpy / scipy
export PYSTAAR_EIGENSOLVER_SIZE_THRESHOLD=256
```

macOS 本地示例（Accelerate/OpenBLAS 对比入口）见：

- `examples/1kg_parity/README.md`

## 6. 数据组织规范

`dataset` 支持目录路径。目录内需要以下文件：

- `geno.mtx`
- `phred.csv`
- `pheno_unrelated.csv`
- `pheno_related.csv`
- `kins_sparse.mtx`
- `kins_dense.mtx`

其中 phenotype 文件至少包含列：`Y`, `X1`, `X2`。

## 7. 验证安装

```bash
pytest -q
```

或只跑 parity：

```bash
pytest tests/parity -q
```

发布版（wheel/pip）用户可先做最小 smoke check：

```bash
python - <<'PY'
import pystaar
res = pystaar.staar_unrelated_glm(dataset="example", seed=600, rare_maf_cutoff=0.05)
print("SMOKE_OK", float(res["results_STAAR_O"]))
PY
```

如果你在仓库内做发布前验收，还可以运行：

```bash
python scripts/run_release_smoke_checks.py
```

## 8. 常见安装问题

### 8.1 `ModuleNotFoundError: pystaar`

通常是未执行 `pip install -e .`，或当前环境不是安装时的虚拟环境。

### 8.2 稀疏矩阵读取失败

请确认 `.mtx` 文件格式正确，且文件路径没有拼写错误。

### 8.3 版本冲突

优先新建干净虚拟环境，再重新安装。

## 9. 仓库维护者发布（PyPI/TestPyPI）

仓库已提供发布 workflow：`.github/workflows/publish.yml`。

- 自动正式发布：
  - 在 GitHub 创建并发布一个 Release（`published`）后，workflow 会构建并发布到 PyPI。
- 手动发布预演：
  - 在 Actions 中手动触发 `Publish` workflow，`target=testpypi` 发布到 TestPyPI。
- 手动正式发布：
  - 在 Actions 中手动触发 `Publish` workflow，`target=pypi` 发布到 PyPI。

发布前请先在仓库 Settings 中配置 GitHub Environments（`testpypi` / `pypi`）并按 Trusted Publishing 绑定对应 Python package。
