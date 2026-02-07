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
