from pathlib import Path

from pystaar import data


REQUIRED_EXAMPLE_FILES = (
    "example_geno.mtx",
    "example_phred.csv",
    "example_pheno_unrelated.csv",
    "example_pheno_related.csv",
    "example_kins_sparse.mtx",
    "example_kins_dense.mtx",
)


def test_data_dir_contains_required_example_files():
    missing = [name for name in REQUIRED_EXAMPLE_FILES if not (data.DATA_DIR / name).exists()]
    assert not missing, f"DATA_DIR is missing required files: {missing}"


def test_packaged_data_bundle_contains_required_example_files():
    bundled_dir = Path(data.__file__).resolve().parent / "_data"
    missing = [name for name in REQUIRED_EXAMPLE_FILES if not (bundled_dir / name).exists()]
    assert not missing, f"Bundled package data is missing files: {missing}"
