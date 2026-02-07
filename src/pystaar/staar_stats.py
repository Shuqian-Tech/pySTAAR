"""Statistical helpers for STAAR computations."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy import special, stats


def _normalize_weights(weights: np.ndarray, expected_size: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=float).reshape(-1)
    if weights.size != expected_size:
        raise ValueError("The length of weights should be the same as that of the p-values!")
    if np.any(np.isnan(weights)):
        raise ValueError("Cannot have NAs in the weights!")
    if np.any(weights < 0):
        raise ValueError("All the weights must be non-negative!")

    weight_sum = float(np.sum(weights))
    if (not np.isfinite(weight_sum)) or weight_sum <= 0.0:
        raise ValueError("Sum of weights must be positive!")
    return weights / weight_sum


def cct_pval(pvals: np.ndarray, weights: np.ndarray) -> float:
    """Cauchy combination test with weights (matches CCT_pval.cpp)."""
    pvals = np.asarray(pvals, dtype=float).reshape(-1)
    if pvals.size == 0:
        raise ValueError("At least one p-value is required.")
    if np.any(np.isnan(pvals)):
        raise ValueError("Cannot have NAs in the p-values!")
    if np.any(pvals < 0) or np.any(pvals > 1):
        raise ValueError("All p-values must be between 0 and 1!")

    weights = _normalize_weights(weights, expected_size=pvals.size)

    cct_stat = 0.0
    for p, w in zip(pvals, weights):
        if p < 1e-16:
            cct_stat += w / p / math.pi
        else:
            cct_stat += w * math.tan((0.5 - p) * math.pi)

    if cct_stat > 1e15:
        return (1.0 / cct_stat) / math.pi
    return float(stats.cauchy.sf(cct_stat, loc=0.0, scale=1.0))


def cct(pvals: Iterable[float], weights: Iterable[float] | None = None) -> float:
    """Cauchy combination test (matches CCT.R)."""
    pvals = np.asarray(list(pvals), dtype=float).reshape(-1)
    if pvals.size == 0:
        raise ValueError("At least one p-value is required.")

    if np.any(np.isnan(pvals)):
        raise ValueError("Cannot have NAs in the p-values!")
    if np.any(pvals < 0) or np.any(pvals > 1):
        raise ValueError("All p-values must be between 0 and 1!")

    is_zero = np.any(pvals == 0)
    is_one = np.any(pvals == 1)
    if is_zero and is_one:
        raise ValueError("Cannot have both 0 and 1 p-values!")
    if is_zero:
        return 0.0
    if is_one:
        return 1.0

    if weights is None:
        weights = np.full_like(pvals, 1.0 / pvals.size)
    else:
        weights = _normalize_weights(list(weights), expected_size=pvals.size)

    is_small = pvals < 1e-16
    if not np.any(is_small):
        cct_stat = np.sum(weights * np.tan((0.5 - pvals) * math.pi))
    else:
        cct_stat = np.sum(weights[is_small] / pvals[is_small] / math.pi)
        cct_stat += np.sum(weights[~is_small] * np.tan((0.5 - pvals[~is_small]) * math.pi))

    if cct_stat > 1e15:
        return (1.0 / cct_stat) / math.pi
    return float(stats.cauchy.sf(cct_stat, loc=0.0, scale=1.0))


def _K(x: float, eigenvals: np.ndarray) -> float:
    return float(-0.5 * np.sum(np.log(1.0 - 2.0 * eigenvals * x)))


def _K1(x: float, eigenvals: np.ndarray, q: float) -> float:
    return float(np.sum(eigenvals / (1.0 - 2.0 * eigenvals * x)) - q)


def _K2(x: float, eigenvals: np.ndarray) -> float:
    return float(2.0 * np.sum((eigenvals**2) / (1.0 - 2.0 * eigenvals * x) ** 2))


def _bisection(eigenvals: np.ndarray, q: float, xmin: float, xmax: float) -> float:
    xupper = xmax
    xlower = xmin
    x0 = 0.0
    while abs(xupper - xlower) > 1e-8:
        x0 = (xupper + xlower) / 2.0
        k1 = _K1(x0, eigenvals, q)
        if k1 == 0:
            break
        if k1 > 0:
            xupper = x0
        else:
            xlower = x0
    return x0


def saddle(q: float, eigenvals: np.ndarray) -> float:
    eigenvals = np.asarray(eigenvals, dtype=float)
    lambdamax = float(np.max(eigenvals))
    q = q / lambdamax
    eigenvals = eigenvals / lambdamax
    lambdamax = 1.0

    if q > np.sum(eigenvals):
        xmin = -0.01
    else:
        xmin = -len(eigenvals) / (2.0 * q)

    xmax = 1.0 / (2.0 * lambdamax) * 0.99999

    xhat = _bisection(eigenvals, q, xmin, xmax)
    w = math.sqrt(2.0 * (xhat * q - _K(xhat, eigenvals)))
    if xhat < 0:
        w = -w
    v = xhat * math.sqrt(_K2(xhat, eigenvals))

    if abs(xhat) < 1e-4:
        return 2.0

    z = w + math.log(v / w) / w
    return float(special.ndtr(-z))
