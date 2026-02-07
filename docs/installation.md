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
pip install -e .[dev]
```

## 3. 依赖版本

项目依赖在 `pyproject.toml` 中定义，核心依赖包括：

- `numpy>=1.24`
- `scipy>=1.10`
- `pandas>=2.0`
- `pyyaml>=6.0`

## 4. 数据组织规范

`dataset` 支持目录路径。目录内需要以下文件：

- `geno.mtx`
- `phred.csv`
- `pheno_unrelated.csv`
- `pheno_related.csv`
- `kins_sparse.mtx`
- `kins_dense.mtx`

其中 phenotype 文件至少包含列：`Y`, `X1`, `X2`。

## 5. 验证安装

```bash
pytest -q
```

或只跑 parity：

```bash
pytest tests/parity -q
```

## 6. 常见安装问题

### 6.1 `ModuleNotFoundError: pystaar`

通常是未执行 `pip install -e .`，或当前环境不是安装时的虚拟环境。

### 6.2 稀疏矩阵读取失败

请确认 `.mtx` 文件格式正确，且文件路径没有拼写错误。

### 6.3 版本冲突

优先新建干净虚拟环境，再重新安装。
