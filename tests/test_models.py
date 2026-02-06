import numpy as np
import pytest

from pystaar.data import load_example_dataset
from pystaar.models import (
    fit_null_glm,
    fit_null_glm_binary_spa,
    fit_null_glmmkin,
    fit_null_glmmkin_binary_spa,
)


def _renamed_pheno(df):
    return df.rename(columns={"Y": "trait", "X1": "cov1", "X2": "cov2"})


def test_fit_null_glm_accepts_custom_columns():
    data = load_example_dataset()
    df = _renamed_pheno(data.pheno_unrelated.copy())

    obj = fit_null_glm(
        df,
        outcome_col="trait",
        covariate_cols=("cov1", "cov2"),
    )

    assert obj.X.shape == (df.shape[0], 3)
    np.testing.assert_allclose(obj.y, df["trait"].to_numpy(dtype=float))


def test_fit_null_glmkin_accepts_custom_columns():
    data = load_example_dataset()
    df = _renamed_pheno(data.pheno_related.copy())

    obj = fit_null_glmmkin(
        df,
        kins=data.kins_sparse,
        sparse_kins=True,
        outcome_col="trait",
        covariate_cols=("cov1", "cov2"),
    )

    assert obj.X.shape == (df.shape[0], 3)
    np.testing.assert_allclose(obj.y, df["trait"].to_numpy(dtype=float))


def test_binary_null_models_accept_custom_columns():
    data = load_example_dataset()

    df_unrelated = _renamed_pheno(data.pheno_unrelated.copy())
    threshold_unrelated = float(np.quantile(df_unrelated["trait"].to_numpy(dtype=float), 0.95))
    df_unrelated["trait"] = (df_unrelated["trait"].to_numpy(dtype=float) > threshold_unrelated).astype(
        float
    )
    obj_glm = fit_null_glm_binary_spa(
        df_unrelated,
        outcome_col="trait",
        covariate_cols=("cov1", "cov2"),
    )
    assert obj_glm.X.shape == (df_unrelated.shape[0], 3)

    df_related = _renamed_pheno(data.pheno_related.copy())
    threshold_related = float(np.quantile(df_related["trait"].to_numpy(dtype=float), 0.95))
    df_related["trait"] = (df_related["trait"].to_numpy(dtype=float) > threshold_related).astype(
        float
    )
    obj_glmmkin = fit_null_glmmkin_binary_spa(
        df_related,
        kins=data.kins_sparse,
        sparse_kins=True,
        outcome_col="trait",
        covariate_cols=("cov1", "cov2"),
    )
    assert obj_glmmkin.X.shape == (df_related.shape[0], 3)


def test_design_matrix_validation_errors():
    data = load_example_dataset()
    df = data.pheno_unrelated.copy()

    with pytest.raises(ValueError, match="Outcome column"):
        fit_null_glm(df, outcome_col="not_present")

    with pytest.raises(ValueError, match="Covariate columns not found"):
        fit_null_glm(df, covariate_cols=("X1", "missing_col"))

    with pytest.raises(ValueError, match="Design matrix must include at least one column"):
        fit_null_glm(df, covariate_cols=(), add_intercept=False)
