"""Null model fitting utilities for STAAR migration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import minimize, minimize_scalar
from scipy.sparse.linalg import splu
from scipy.special import expit


@dataclass
class NullModelGLM:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    beta: np.ndarray
    dispersion: float
    weights: np.ndarray
    family: str
    relatedness: bool = False
    use_SPA: bool = False
    XW: np.ndarray | None = None
    XXWX_inv: np.ndarray | None = None


@dataclass
class NullModelGLMMKin:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    scaled_residuals: np.ndarray
    beta: np.ndarray
    theta: np.ndarray
    cov: np.ndarray
    sigma_solver: splu
    Sigma_iX: np.ndarray
    sparse_kins: bool
    relatedness: bool = True
    precomputed_cov: np.ndarray | None = None


@dataclass
class NullModelGLMMKinBinarySPA:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    scaled_residuals: np.ndarray
    beta: np.ndarray
    theta: np.ndarray
    cov: np.ndarray
    sigma_solver: splu
    Sigma_iX: np.ndarray
    sparse_kins: bool
    weights: np.ndarray
    XW: np.ndarray
    XXWX_inv: np.ndarray
    relatedness: bool = True
    use_SPA: bool = True
    family: str = "binomial"
    precomputed_cov_filter: np.ndarray | None = None


def _design_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    y = df["Y"].to_numpy(dtype=float)
    X = np.column_stack(
        [np.ones(len(df)), df["X1"].to_numpy(dtype=float), df["X2"].to_numpy(dtype=float)]
    )
    return X, y


def _fit_binomial_irls(
    X: np.ndarray,
    y: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    beta = np.zeros(X.shape[1], dtype=float)
    eps = 1e-12

    for _ in range(max_iter):
        eta = X @ beta
        mu = np.clip(expit(eta), eps, 1.0 - eps)
        w = np.clip(mu * (1.0 - mu), eps, None)
        z = eta + (y - mu) / w

        WX = X * w[:, None]
        XtWX = X.T @ WX
        XtWz = X.T @ (w * z)

        beta_new = np.linalg.solve(XtWX, XtWz)
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    fitted = np.clip(expit(X @ beta), eps, 1.0 - eps)
    weights = np.clip(fitted * (1.0 - fitted), eps, None)
    return beta, fitted, weights


def fit_null_glm(df: pd.DataFrame) -> NullModelGLM:
    X, y = _design_matrix(df)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    residuals = y - fitted
    df_resid = X.shape[0] - X.shape[1]
    dispersion = float(np.sum(residuals**2) / df_resid)
    weights = np.ones_like(y)
    return NullModelGLM(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        beta=beta,
        dispersion=dispersion,
        weights=weights,
        family="gaussian",
    )


def fit_null_glm_binary_spa(
    df: pd.DataFrame,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> NullModelGLM:
    X, y = _design_matrix(df)
    y = y.astype(float)
    if not np.all(np.isin(y, [0.0, 1.0])):
        raise ValueError("Binary SPA null model requires Y coded as 0/1.")

    beta, fitted, weights = _fit_binomial_irls(X=X, y=y, max_iter=max_iter, tol=tol)
    residuals = y - fitted

    XtWX = X.T @ (X * weights[:, None])
    XXWX_inv = X @ np.linalg.inv(XtWX)
    XW = (X * weights[:, None]).T

    return NullModelGLM(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        beta=beta,
        dispersion=1.0,
        weights=weights,
        family="binomial",
        use_SPA=True,
        XW=XW,
        XXWX_inv=XXWX_inv,
    )


def _reml_nll(log_tau: np.ndarray, X: np.ndarray, y: np.ndarray, kins: sp.csc_matrix) -> float:
    tau1 = float(np.exp(log_tau[0]))
    tau2 = float(np.exp(log_tau[1]))
    n = X.shape[0]
    sigma = sp.identity(n, format="csc")
    sigma.data *= tau1
    if tau2 != 0.0:
        sigma = sigma + kins * tau2
    lu = splu(sigma)
    diag_u = lu.U.diagonal()
    logdet = float(np.sum(np.log(np.abs(diag_u))))

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    sign, logdet_x = np.linalg.slogdet(xt_sigma_i_x)
    if sign <= 0:
        return np.inf
    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_y = lu.solve(y)
    beta = cov @ (X.T @ sigma_i_y)
    residual = y - X @ beta
    sigma_i_r = lu.solve(residual)
    quad = float(residual @ sigma_i_r)

    p = X.shape[1]
    ll = -0.5 * (logdet + logdet_x + quad + (n - p) * np.log(2 * np.pi))
    return -ll


def _reml_nll_binomial_tau(
    log_tau: float,
    X: np.ndarray,
    z: np.ndarray,
    d_inv: np.ndarray,
    kins: sp.csc_matrix,
) -> float:
    tau = float(np.exp(log_tau))
    sigma = sp.diags(d_inv, format="csc")
    if tau != 0.0:
        sigma = sigma + kins * tau
    lu = splu(sigma)

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    sign, logdet_x = np.linalg.slogdet(xt_sigma_i_x)
    if sign <= 0:
        return np.inf

    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_z = lu.solve(z)
    beta = cov @ (X.T @ sigma_i_z)
    residual = z - X @ beta
    sigma_i_r = lu.solve(residual)
    quad = float(residual @ sigma_i_r)

    diag_u = lu.U.diagonal()
    logdet = float(np.sum(np.log(np.abs(diag_u))))
    return 0.5 * (logdet + logdet_x + quad)


def estimate_tau_reml(X: np.ndarray, y: np.ndarray, kins: sp.csc_matrix) -> np.ndarray:
    init = np.log([1.0, 1.0])
    result = minimize(
        _reml_nll,
        init,
        args=(X, y, kins),
        method="Nelder-Mead",
        options={"maxiter": 200, "xatol": 1e-8, "fatol": 1e-8},
    )
    return np.exp(result.x)


def fit_null_glmmkin(
    df: pd.DataFrame,
    kins: sp.csc_matrix,
    sparse_kins: bool,
    precomputed_cov: np.ndarray | None = None,
    precomputed_scaled_residuals: np.ndarray | None = None,
    precomputed_theta: np.ndarray | None = None,
) -> NullModelGLMMKin:
    X, y = _design_matrix(df)
    if precomputed_theta is None:
        theta = estimate_tau_reml(X, y, kins)
    else:
        theta = np.asarray(precomputed_theta, dtype=float).reshape(-1)
        if theta.size != 2:
            raise ValueError("precomputed_theta must contain exactly two components.")
        if np.any(theta <= 0):
            raise ValueError("precomputed_theta components must be positive.")

    n = X.shape[0]
    sigma = sp.identity(n, format="csc")
    sigma.data *= theta[0]
    sigma = sigma + kins * theta[1]
    lu = splu(sigma)

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_y = lu.solve(y)
    beta = cov @ (X.T @ sigma_i_y)

    # Match GMMAT's eta computation in R_fitglmm_ai
    diag_sigma = theta[0]
    fitted = y - diag_sigma * (sigma_i_y - sigma_i_x @ beta)
    residuals = y - fitted
    if precomputed_scaled_residuals is None:
        scaled_residuals = residuals / theta[0]
    else:
        scaled_residuals = precomputed_scaled_residuals

    return NullModelGLMMKin(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        scaled_residuals=scaled_residuals,
        beta=beta,
        theta=theta,
        cov=cov,
        sigma_solver=lu,
        Sigma_iX=sigma_i_x,
        sparse_kins=sparse_kins,
        precomputed_cov=precomputed_cov,
    )


def fit_null_glmmkin_binary_spa(
    df: pd.DataFrame,
    kins: sp.csc_matrix,
    sparse_kins: bool,
    max_iter: int = 20,
    tol: float = 1e-5,
) -> NullModelGLMMKinBinarySPA:
    X, y = _design_matrix(df)
    y = y.astype(float)
    if not np.all(np.isin(y, [0.0, 1.0])):
        raise ValueError("Binary SPA null model requires Y coded as 0/1.")

    eps = 1e-9
    beta, _, _ = _fit_binomial_irls(X=X, y=y, max_iter=max_iter, tol=tol)
    eta = X @ beta
    tau = 1.0

    for _ in range(max_iter):
        fitted = np.clip(expit(eta), eps, 1.0 - eps)
        weights = np.clip(fitted * (1.0 - fitted), eps, None)
        z = eta + (y - fitted) / weights
        d_inv = 1.0 / weights

        try:
            opt = minimize_scalar(
                _reml_nll_binomial_tau,
                bounds=(-10.0, 5.0),
                method="bounded",
                args=(X, z, d_inv, kins),
                options={"xatol": 1e-3, "maxiter": 80},
            )
            if opt.success:
                tau_new = float(np.exp(opt.x))
            else:
                tau_new = tau
        except Exception:
            tau_new = tau

        sigma = sp.diags(d_inv, format="csc")
        if tau_new != 0.0:
            sigma = sigma + kins * tau_new
        lu = splu(sigma)

        sigma_i_x = lu.solve(X)
        xt_sigma_i_x = X.T @ sigma_i_x
        cov = np.linalg.inv(xt_sigma_i_x)
        sigma_i_z = lu.solve(z)
        beta_new = cov @ (X.T @ sigma_i_z)
        residual_working = z - X @ beta_new
        u = tau_new * (kins @ lu.solve(residual_working))
        eta_new = X @ beta_new + u

        delta = float(np.max(np.abs(eta_new - eta)))
        eta = eta_new
        beta = beta_new
        tau = tau_new
        if delta < tol:
            break

    fitted = np.clip(expit(eta), eps, 1.0 - eps)
    residuals = y - fitted
    weights = np.clip(fitted * (1.0 - fitted), eps, None)
    d_inv = 1.0 / weights
    sigma = sp.diags(d_inv, format="csc")
    if tau != 0.0:
        sigma = sigma + kins * tau
    lu = splu(sigma)
    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    cov = np.linalg.inv(xt_sigma_i_x)

    # Keep residuals on the Y-mu scale for SPA score construction.
    scaled_residuals = residuals

    if sparse_kins:
        XW = sigma_i_x.T
        XXWX_inv = X @ cov
    else:
        XtWX = X.T @ (X * weights[:, None])
        XXWX_inv = X @ np.linalg.inv(XtWX)
        XW = (X * weights[:, None]).T

    theta = np.array([1.0, tau], dtype=float)

    return NullModelGLMMKinBinarySPA(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        scaled_residuals=scaled_residuals,
        beta=beta,
        theta=theta,
        cov=cov,
        sigma_solver=lu,
        Sigma_iX=sigma_i_x,
        sparse_kins=sparse_kins,
        weights=weights,
        XW=XW,
        XXWX_inv=XXWX_inv,
    )
