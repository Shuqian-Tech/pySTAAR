"""Data loading utilities for STAAR migration datasets."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread


DATA_DIR = Path(__file__).resolve().parents[2] / "data"


@dataclass
class STAARDataset:
    geno: np.ndarray
    phred: np.ndarray
    pheno_unrelated: pd.DataFrame
    pheno_related: pd.DataFrame
    kins_sparse: sp.csc_matrix
    kins_dense: sp.csc_matrix


ExampleDataset = STAARDataset


def _load_sparse(path: Path) -> sp.csc_matrix:
    mat = mmread(path)
    if not sp.issparse(mat):
        mat = sp.csc_matrix(mat)
    else:
        mat = mat.tocsc()
    return mat


def _validate_pheno_columns(df: pd.DataFrame, source: Path) -> pd.DataFrame:
    required = ("Y", "X1", "X2")
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"Phenotype file '{source}' is missing required columns: {', '.join(missing)}."
        )
    return df


def _load_dataset_from_paths(
    *,
    dataset_label: str,
    geno_path: Path,
    phred_path: Path,
    pheno_unrelated_path: Path,
    pheno_related_path: Path,
    kins_sparse_path: Path,
    kins_dense_path: Path,
) -> STAARDataset:
    required_paths = [
        geno_path,
        phred_path,
        pheno_unrelated_path,
        pheno_related_path,
        kins_sparse_path,
        kins_dense_path,
    ]
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise ValueError(
            f"Dataset '{dataset_label}' is missing required files: {', '.join(missing)}"
        )

    geno = _load_sparse(geno_path).toarray().astype(float)
    phred = pd.read_csv(phred_path).to_numpy(dtype=float)
    pheno_unrelated = _validate_pheno_columns(pd.read_csv(pheno_unrelated_path), pheno_unrelated_path)
    pheno_related = _validate_pheno_columns(pd.read_csv(pheno_related_path), pheno_related_path)
    kins_sparse = _load_sparse(kins_sparse_path)
    kins_dense = _load_sparse(kins_dense_path)

    return STAARDataset(
        geno=geno,
        phred=phred,
        pheno_unrelated=pheno_unrelated,
        pheno_related=pheno_related,
        kins_sparse=kins_sparse,
        kins_dense=kins_dense,
    )


def load_named_dataset(name: str) -> STAARDataset:
    return _load_named_dataset_cached(name)


@lru_cache(maxsize=16)
def _load_named_dataset_cached(name: str) -> STAARDataset:
    base = DATA_DIR
    return _load_dataset_from_paths(
        dataset_label=name,
        geno_path=base / f"{name}_geno.mtx",
        phred_path=base / f"{name}_phred.csv",
        pheno_unrelated_path=base / f"{name}_pheno_unrelated.csv",
        pheno_related_path=base / f"{name}_pheno_related.csv",
        kins_sparse_path=base / f"{name}_kins_sparse.mtx",
        kins_dense_path=base / f"{name}_kins_dense.mtx",
    )


def load_dataset_from_directory(path: str | Path) -> STAARDataset:
    directory = Path(path).expanduser().resolve()
    if not directory.exists() or not directory.is_dir():
        raise ValueError(f"Dataset directory does not exist: {directory}")
    return _load_dataset_from_directory_cached(str(directory))


@lru_cache(maxsize=16)
def _load_dataset_from_directory_cached(directory_str: str) -> STAARDataset:
    directory = Path(directory_str)
    return _load_dataset_from_paths(
        dataset_label=directory_str,
        geno_path=directory / "geno.mtx",
        phred_path=directory / "phred.csv",
        pheno_unrelated_path=directory / "pheno_unrelated.csv",
        pheno_related_path=directory / "pheno_related.csv",
        kins_sparse_path=directory / "kins_sparse.mtx",
        kins_dense_path=directory / "kins_dense.mtx",
    )


def load_example_dataset() -> ExampleDataset:
    return load_named_dataset("example")


def _cache_info_to_dict(cache_info) -> dict[str, int]:
    return {
        "hits": int(cache_info.hits),
        "misses": int(cache_info.misses),
        "maxsize": int(cache_info.maxsize),
        "currsize": int(cache_info.currsize),
    }


def get_dataset_cache_info() -> dict[str, dict[str, int]]:
    return {
        "named_dataset": _cache_info_to_dict(_load_named_dataset_cached.cache_info()),
        "directory_dataset": _cache_info_to_dict(_load_dataset_from_directory_cached.cache_info()),
    }


def clear_dataset_cache() -> dict[str, dict[str, dict[str, int]]]:
    before = get_dataset_cache_info()
    _load_named_dataset_cached.cache_clear()
    _load_dataset_from_directory_cached.cache_clear()
    after = get_dataset_cache_info()
    return {"before": before, "after": after}


def load_dataset(dataset: str | STAARDataset) -> STAARDataset:
    if isinstance(dataset, STAARDataset):
        return dataset
    if isinstance(dataset, str):
        if dataset == "example":
            return load_example_dataset()
        dataset_path = Path(dataset).expanduser()
        if dataset_path.exists() and dataset_path.is_dir():
            return load_dataset_from_directory(dataset_path)
        return load_named_dataset(dataset)
    raise TypeError(
        "dataset must be a dataset name, a dataset directory path, or a STAARDataset object."
    )
