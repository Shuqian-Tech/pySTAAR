#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) {
    return(args[i + 1])
  }
  default
}

get_arg_int <- function(key, default = NULL) {
  value <- get_arg(key, NULL)
  if (is.null(value)) {
    return(default)
  }
  as.integer(value)
}

get_arg_num <- function(key, default = NULL) {
  value <- get_arg(key, NULL)
  if (is.null(value)) {
    return(default)
  }
  as.numeric(value)
}

input_rda <- get_arg("--input-rda")
outdir <- get_arg("--outdir")
seed <- get_arg_int("--seed", 600L)
target_nsample <- get_arg_int("--target-nsample", 1200L)
target_nvar <- get_arg_int("--target-nvar", 1500L)
sample_maf_max <- get_arg_num("--sample-maf-max", 0.2)
num_annotations <- get_arg_int("--num-annotations", 10L)
kin_sparse_threshold <- get_arg_num("--kin-sparse-threshold", 0.02)

if (is.null(input_rda)) {
  stop("Missing required argument: --input-rda")
}
if (is.null(outdir)) {
  stop("Missing required argument: --outdir")
}
if (target_nsample < 2) {
  stop("--target-nsample must be >= 2")
}
if (target_nvar < 10) {
  stop("--target-nvar must be >= 10")
}
if (num_annotations < 1) {
  stop("--num-annotations must be >= 1")
}
if (sample_maf_max <= 0 || sample_maf_max > 0.5) {
  stop("--sample-maf-max must be in (0, 0.5]")
}
if (kin_sparse_threshold < 0) {
  stop("--kin-sparse-threshold must be >= 0")
}

input_rda <- normalizePath(input_rda, mustWork = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)

e <- new.env(parent = emptyenv())
load(input_rda, envir = e)
obj_names <- ls(e)
if (length(obj_names) == 0) {
  stop("No objects found in RDA: ", input_rda)
}

source_obj <- NULL
source_name <- NULL
for (nm in obj_names) {
  cand <- get(nm, envir = e)
  if (
    is.list(cand) &&
      all(c("genotype", "snploc", "maf") %in% names(cand))
  ) {
    source_obj <- cand
    source_name <- nm
    break
  }
}
if (is.null(source_obj)) {
  stop("No list object with genotype/snploc/maf found in ", input_rda)
}

geno_all <- source_obj$genotype
if (!methods::is(geno_all, "dgCMatrix")) {
  geno_all <- as(geno_all, "dgCMatrix")
}
snploc_all <- as.data.frame(source_obj$snploc, stringsAsFactors = FALSE)
maf_all <- as.numeric(source_obj$maf)

if (length(maf_all) != ncol(geno_all)) {
  stop("MAF length does not match genotype columns.")
}
if (nrow(snploc_all) != ncol(geno_all)) {
  stop("snploc rows do not match genotype columns.")
}

variant_ids <- colnames(geno_all)
if (is.null(variant_ids) || length(variant_ids) != ncol(geno_all)) {
  variant_ids <- names(source_obj$maf)
}
if (is.null(variant_ids) || length(variant_ids) != ncol(geno_all)) {
  variant_ids <- as.character(seq_len(ncol(geno_all)))
}

set.seed(seed)
sample_idx <- seq_len(nrow(geno_all))
if (target_nsample < length(sample_idx)) {
  sample_idx <- sort(sample(sample_idx, size = target_nsample, replace = FALSE))
}
geno_sample <- geno_all[sample_idx, , drop = FALSE]

maf_candidate <- which(!is.na(maf_all) & maf_all > 0 & maf_all <= sample_maf_max)
if (length(maf_candidate) == 0) {
  stop("No variants available under sample_maf_max=", sample_maf_max)
}

variant_idx <- maf_candidate
if (target_nvar < length(variant_idx)) {
  variant_idx <- sort(sample(variant_idx, size = target_nvar, replace = FALSE))
}

geno <- geno_sample[, variant_idx, drop = FALSE]
snploc <- snploc_all[variant_idx, , drop = FALSE]
maf <- maf_all[variant_idx]
variant_ids_sel <- variant_ids[variant_idx]

if (!is.null(rownames(geno_all))) {
  rownames(geno) <- rownames(geno_all)[sample_idx]
}
colnames(geno) <- variant_ids_sel

n_samples <- nrow(geno)
n_variants <- ncol(geno)

set.seed(seed + 101L)
z <- matrix(
  rnorm(n_variants * num_annotations),
  nrow = n_variants,
  ncol = num_annotations
)
rank_mat <- apply(z, 2, rank, ties.method = "average")
phred_base <- 1 - rank_mat / (n_variants + 1)
phred_base[phred_base < 1e-15] <- 1e-15
phred <- -10 * log10(phred_base)
colnames(phred) <- paste0("Z", seq_len(num_annotations))

geno_dense <- as.matrix(geno)
maf_est <- colMeans(geno_dense, na.rm = TRUE) / 2

# Stable pedigree-style kinship blocks (family size = 4),
# aligned with the baseline simulation style used in this repo.
unitmat <- matrix(0.5, nrow = 4, ncol = 4)
diag(unitmat) <- 1
unitmat[1, 2] <- 0
unitmat[2, 1] <- 0

kins_dense <- diag(n_samples)
n_family <- n_samples %/% 4
if (n_family > 0) {
  for (fam in seq_len(n_family)) {
    idx_start <- (fam - 1) * 4 + 1
    idx_end <- idx_start + 3
    kins_dense[idx_start:idx_end, idx_start:idx_end] <- unitmat
  }
}
diag(kins_dense) <- diag(kins_dense) + 1e-6

kins_sparse_mat <- kins_dense
if (kin_sparse_threshold > 0) {
  kins_sparse_mat[abs(kins_sparse_mat) < kin_sparse_threshold] <- 0
}
diag(kins_sparse_mat) <- diag(kins_dense)
kins_sparse <- Matrix(kins_sparse_mat, sparse = TRUE)

sample_ids <- as.character(seq_len(n_samples))
rownames(geno) <- sample_ids
dimnames(kins_dense) <- list(sample_ids, sample_ids)
dimnames(kins_sparse) <- list(sample_ids, sample_ids)

set.seed(seed + 202L)
x1 <- rnorm(n_samples)
x2 <- rbinom(n_samples, 1, 0.5)
eps <- rnorm(n_samples)

set.seed(seed + 303L)
causal_count <- min(max(10L, as.integer(floor(0.02 * n_variants))), n_variants)
causal_idx <- sort(sample(seq_len(n_variants), size = causal_count, replace = FALSE))
maf_safe <- pmax(maf_est, 1e-4)
beta <- -0.1 * log10(maf_safe[causal_idx])
genetic_effect <- as.vector(geno_dense[, causal_idx, drop = FALSE] %*% beta)
genetic_effect <- as.numeric(scale(genetic_effect))

set.seed(seed + 404L)
u <- numeric(n_samples)
if (n_family > 0) {
  chol_unit <- chol(0.5 * unitmat)
  for (fam in seq_len(n_family)) {
    idx_start <- (fam - 1) * 4 + 1
    idx_end <- idx_start + 3
    z <- rnorm(4)
    u[idx_start:idx_end] <- as.vector(z %*% chol_unit)
  }
}
if (n_family * 4 < n_samples) {
  rem_idx <- (n_family * 4 + 1):n_samples
  u[rem_idx] <- rnorm(length(rem_idx), sd = sqrt(0.5))
}
u <- as.numeric(scale(u))

y_unrelated <- 0.2 + 0.5 * x1 + 0.4 * x2 + genetic_effect + eps
y_related <- 0.2 + 0.5 * x1 + 0.4 * x2 + genetic_effect + 0.5 * u + eps

pheno_unrelated <- data.frame(
  Y = as.numeric(y_unrelated),
  X1 = as.numeric(x1),
  X2 = as.numeric(x2),
  stringsAsFactors = FALSE
)
pheno_related <- data.frame(
  Y = as.numeric(y_related),
  X1 = as.numeric(x1),
  X2 = as.numeric(x2),
  id = sample_ids,
  stringsAsFactors = FALSE
)

invisible(writeMM(geno, file.path(outdir, "geno.mtx")))
write.csv(phred, file.path(outdir, "phred.csv"), row.names = FALSE, quote = FALSE)
write.csv(pheno_unrelated, file.path(outdir, "pheno_unrelated.csv"), row.names = FALSE, quote = FALSE)
write.csv(pheno_related, file.path(outdir, "pheno_related.csv"), row.names = FALSE, quote = FALSE)
invisible(writeMM(kins_sparse, file.path(outdir, "kins_sparse.mtx")))
invisible(writeMM(Matrix(kins_dense, sparse = TRUE), file.path(outdir, "kins_dense.mtx")))

sim_rds_path <- file.path(outdir, "sim_workflow_data.rds")
saveRDS(
  list(
    Geno = geno,
    PHRED = phred,
    pheno_unrelated = pheno_unrelated,
    pheno_related = pheno_related,
    kins_sparse = kins_sparse,
    kins_dense = kins_dense,
    maf = maf,
    snploc = snploc,
    variant_id = variant_ids_sel,
    sample_idx = sample_idx,
    variant_idx = variant_idx
  ),
  file = sim_rds_path
)

summary_payload <- list(
  input_rda = input_rda,
  input_object = source_name,
  seed = as.integer(seed),
  target_nsample = as.integer(target_nsample),
  target_nvar = as.integer(target_nvar),
  sample_maf_max = as.numeric(sample_maf_max),
  num_annotations = as.integer(num_annotations),
  kin_sparse_threshold = as.numeric(kin_sparse_threshold),
  n_samples = as.integer(n_samples),
  n_variants = as.integer(n_variants),
  nnz_geno = as.integer(length(geno@x)),
  nnz_kins_sparse = as.integer(length(kins_sparse@x)),
  sim_rds = sim_rds_path
)

write_json(
  summary_payload,
  path = file.path(outdir, "sim_workflow_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

cat("Simulated workflow dataset created.\n")
cat("  outdir:", outdir, "\n")
cat("  n_samples:", n_samples, " n_variants:", n_variants, "\n")
cat("  sim_rds:", sim_rds_path, "\n")
