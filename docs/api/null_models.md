# API: Null Models

本页描述 null model 拟合函数，供高级用户或需要与 R 流程逐步对齐的用户使用。

## 1. `fit_null_glm`

```python
fit_null_glm(
    df,
    outcome_col="Y",
    covariate_cols=("X1", "X2"),
    add_intercept=True,
)
```

用途：非亲属连续性状 null model。

返回：`NullModelGLM`

关键字段：

- `X`, `y`
- `fitted`, `residuals`
- `beta`, `dispersion`
- `weights`, `family`

## 2. `fit_null_glm_binary_spa`

```python
fit_null_glm_binary_spa(
    df,
    max_iter=100,
    tol=1e-8,
    outcome_col="Y",
    covariate_cols=("X1", "X2"),
    add_intercept=True,
)
```

用途：非亲属二分类 SPA null model。

返回：`NullModelGLM`

注意：`Y` 需为 0/1。

## 3. `fit_null_glmmkin`

```python
fit_null_glmmkin(
    df,
    kins,
    sparse_kins,
    precomputed_cov=None,
    precomputed_scaled_residuals=None,
    precomputed_theta=None,
    outcome_col="Y",
    covariate_cols=("X1", "X2"),
    add_intercept=True,
)
```

用途：相关样本连续性状 GLMM kinship null model。

返回：`NullModelGLMMKin`

关键字段：

- `theta`（`dispersion` 与 `kins` 组件）
- `scaled_residuals`
- `sigma_solver`, `Sigma_iX`, `cov`

## 4. `fit_null_glmmkin_binary_spa`

```python
fit_null_glmmkin_binary_spa(
    df,
    kins,
    sparse_kins,
    max_iter=80,
    tol=1e-8,
    outcome_col="Y",
    covariate_cols=("X1", "X2"),
    add_intercept=True,
)
```

用途：相关样本二分类 SPA null model。

返回：`NullModelGLMMKinBinarySPA`

关键字段：

- `weights`
- `XW`, `XXWX_inv`
- `scaled_residuals`

## 5. R 名称兼容别名

- `fit_null_glm_Binary_SPA` -> `fit_null_glm_binary_spa`
- `fit_null_glmmkin_Binary_SPA` -> `fit_null_glmmkin_binary_spa`

## 6. 推荐实践

多数用户建议直接调用 `workflows.py` 高层入口（如 `staar_unrelated_glm`），除非你明确需要控制 null model 构造细节。
