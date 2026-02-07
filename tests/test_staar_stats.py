import numpy as np
import pytest

from pystaar.staar_stats import cct, cct_pval


def test_cct_rejects_empty_pvalues():
    with pytest.raises(ValueError, match="At least one p-value is required"):
        cct([])


def test_cct_rejects_zero_sum_weights():
    with pytest.raises(ValueError, match="Sum of weights must be positive"):
        cct([0.1, 0.2], weights=[0.0, 0.0])


def test_cct_pval_rejects_zero_sum_weights():
    with pytest.raises(ValueError, match="Sum of weights must be positive"):
        cct_pval(np.array([0.1, 0.2]), np.array([0.0, 0.0]))


def test_cct_pval_ignores_exact_one_pvalues():
    pvals = np.array([1.0, 1e-5, 5e-5], dtype=float)
    weights = np.array([0.2, 0.4, 0.4], dtype=float)

    observed = cct_pval(pvals, weights)
    expected = cct_pval(pvals[1:], weights[1:] / np.sum(weights[1:]))

    assert observed == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_cct_pval_all_one_pvalues_returns_one():
    assert cct_pval(np.array([1.0, 1.0]), np.array([0.4, 0.6])) == pytest.approx(1.0)
