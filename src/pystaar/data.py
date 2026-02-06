"""Data loading utilities for STAAR migration examples."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread


DATA_DIR = Path(__file__).resolve().parents[2] / "data"


@dataclass
class ExampleDataset:
    geno: np.ndarray
    phred: np.ndarray
    pheno_unrelated: pd.DataFrame
    pheno_related: pd.DataFrame
    kins_sparse: sp.csc_matrix
    kins_dense: sp.csc_matrix


def _load_sparse(path: Path) -> sp.csc_matrix:
    mat = mmread(path)
    if not sp.issparse(mat):
        mat = sp.csc_matrix(mat)
    else:
        mat = mat.tocsc()
    return mat


def load_example_dataset() -> ExampleDataset:
    geno = _load_sparse(DATA_DIR / "example_geno.mtx").toarray().astype(float)
    phred = pd.read_csv(DATA_DIR / "example_phred.csv").to_numpy(dtype=float)
    pheno_unrelated = pd.read_csv(DATA_DIR / "example_pheno_unrelated.csv")
    pheno_related = pd.read_csv(DATA_DIR / "example_pheno_related.csv")
    kins_sparse = _load_sparse(DATA_DIR / "example_kins_sparse.mtx")
    kins_dense = _load_sparse(DATA_DIR / "example_kins_dense.mtx")

    return ExampleDataset(
        geno=geno,
        phred=phred,
        pheno_unrelated=pheno_unrelated,
        pheno_related=pheno_related,
        kins_sparse=kins_sparse,
        kins_dense=kins_dense,
    )
