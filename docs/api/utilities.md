# API: 工具函数

## 1. `matrix_flip`

```python
from pystaar import matrix_flip

G_flip, flip_indicator, maf = matrix_flip(G)
```

用途：将基因型矩阵翻转到 MAF <= 0.5 的等价表示，并返回：

- `G_flip`: 翻转后矩阵
- `flip_indicator`: 每个位点是否翻转
- `maf`: 每个位点 MAF

常用于调试位点筛选和权重构建。

## 2. `CCT` / `cct`

```python
from pystaar import CCT

p = CCT([0.01, 0.2, 0.5])
```

用途：Cauchy combination test 聚合 p-value。

注意：

- 输入 p-value 需在 `(0, 1]` 范围。
- 可选 `weights`，长度需与 p-value 一致。

## 3. 数据工具（`pystaar.data`）

常用函数：

- `load_dataset`
- `load_example_dataset`
- `load_named_dataset`
- `load_dataset_from_directory`

### 数据目录格式

若通过目录加载，目录中应包含：

- `geno.mtx`
- `phred.csv`
- `pheno_unrelated.csv`
- `pheno_related.csv`
- `kins_sparse.mtx`
- `kins_dense.mtx`

## 4. 调试建议

- 先用 `load_example_dataset()` 检查流程。
- 如遇 shape 错误，优先检查 `geno/phred/pheno` 行列是否匹配。
