# 数据目录模板（中文）

当你用 `dataset="/path/to/dir"` 运行 workflow 时，目录必须包含以下 6 个文件：

- `geno.mtx`
- `phred.csv`
- `pheno_unrelated.csv`
- `pheno_related.csv`
- `kins_sparse.mtx`
- `kins_dense.mtx`

## 1. 最小目录结构

```text
my_dataset/
  geno.mtx
  phred.csv
  pheno_unrelated.csv
  pheno_related.csv
  kins_sparse.mtx
  kins_dense.mtx
```

## 2. 形状约束（必须一致）

- `geno.mtx`: `n_samples x n_variants`
- `phred.csv`: `n_variants x n_annotations`
- `pheno_unrelated.csv`: `n_samples x ...`（至少包含 `Y`, `X1`, `X2`）
- `pheno_related.csv`: `n_samples x ...`（至少包含 `Y`, `X1`, `X2`）
- `kins_sparse.mtx`: `n_samples x n_samples`
- `kins_dense.mtx`: `n_samples x n_samples`

如果这些维度不一致，后续通常会报 shape 相关错误。

## 3. phenotype 最小列模板

`pheno_unrelated.csv` 与 `pheno_related.csv` 至少包含：

```csv
Y,X1,X2
0.34,1.2,0
1.05,0.3,1
...
```

说明：

- 列名区分大小写。
- workflow 默认用 `Y` 作为表型，`X1/X2` 作为协变量（可在底层 API 自定义）。

## 4. 常见报错与对应排查

| 报错片段 | 常见原因 | 排查动作 |
|---|---|---|
| `missing required files` | 目录文件名不完整或拼写错误 | 检查 6 个固定文件名是否都存在 |
| `missing required columns: Y, X1, X2` | phenotype 缺列或列名大小写不一致 | 补齐列并统一列名 |
| `adj_variant_indices contains out-of-range indices` | 条件位点索引越界 | 检查 `geno` 变异位点数，使用 0-based 索引 |
| `Number of rare variant in the set is less than rv_num_cutoff!` | 当前阈值下稀有位点太少 | 调整 `rare_maf_cutoff`，检查输入数据质量 |

## 5. 一条快速自检命令

完成目录准备后，先跑一个最小调用：

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(
    dataset="/path/to/my_dataset",
    seed=600,
    rare_maf_cutoff=0.05,
)
print(res["results_STAAR_O"])
```
