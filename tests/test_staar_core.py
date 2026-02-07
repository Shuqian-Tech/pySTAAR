import numpy as np

from pystaar import staar_core


def test_matrix_flip_all_missing_columns_match_legacy_behavior():
    genotype = np.full((4, 3), -1.0, dtype=float)

    flipped, af, maf = staar_core.matrix_flip(genotype)

    np.testing.assert_allclose(flipped, np.zeros_like(flipped))
    assert np.isnan(af).all()
    assert np.isnan(maf).all()
