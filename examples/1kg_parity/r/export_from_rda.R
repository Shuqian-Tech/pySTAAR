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

input_rda <- get_arg("--input-rda")
outdir <- get_arg("--outdir", ".")

if (is.null(input_rda)) {
  stop("Missing required argument: --input-rda")
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

obj <- NULL
obj_name <- NULL
for (nm in obj_names) {
  cand <- get(nm, envir = e)
  if (
    is.list(cand) &&
      all(c("genotype", "snploc", "maf") %in% names(cand))
  ) {
    obj <- cand
    obj_name <- nm
    break
  }
}
if (is.null(obj)) {
  stop("No list object with genotype/snploc/maf found in ", input_rda)
}

genotype <- obj$genotype
if (!methods::is(genotype, "dgCMatrix")) {
  genotype <- as(genotype, "dgCMatrix")
}
snploc <- as.data.frame(obj$snploc, stringsAsFactors = FALSE)
maf <- as.numeric(obj$maf)
maf_names <- names(obj$maf)

if (length(maf) != ncol(genotype)) {
  stop("MAF length does not match genotype columns.")
}
if (nrow(snploc) != ncol(genotype)) {
  stop("snploc rows do not match genotype columns.")
}

if (is.null(maf_names) || length(maf_names) != length(maf)) {
  maf_names <- colnames(genotype)
}
if (is.null(maf_names) || length(maf_names) != length(maf)) {
  maf_names <- as.character(seq_along(maf))
}

snploc_path <- file.path(outdir, "snploc.csv.gz")
maf_path <- file.path(outdir, "maf.csv.gz")
geno_path <- file.path(outdir, "genotype.mtx")
r_summary_path <- file.path(outdir, "r_summary.json")

write.csv(snploc, gzfile(snploc_path), row.names = FALSE, quote = TRUE)
maf_df <- data.frame(
  variant_id = maf_names,
  maf = maf,
  stringsAsFactors = FALSE
)
write.csv(maf_df, gzfile(maf_path), row.names = FALSE, quote = FALSE)
invisible(writeMM(genotype, geno_path))

matrix_market_path <- geno_path
gzip_bin <- Sys.which("gzip")
if (nzchar(gzip_bin)) {
  status <- system2(gzip_bin, c("-f", geno_path))
  if (identical(status, 0L)) {
    matrix_market_path <- paste0(geno_path, ".gz")
  }
}

maf_q <- stats::quantile(maf, probs = c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE)
maf_summary <- as.list(unname(maf_q))
names(maf_summary) <- c("min", "q01", "median", "q99", "max")

payload <- list(
  input_rda = input_rda,
  object_name = obj_name,
  n_samples = as.integer(nrow(genotype)),
  n_variants = as.integer(ncol(genotype)),
  nnz = as.integer(length(genotype@x)),
  genotype_sum = as.numeric(sum(genotype@x)),
  maf_mean = as.numeric(mean(maf)),
  maf_summary = maf_summary,
  snploc_columns = colnames(snploc),
  files = list(
    snploc = snploc_path,
    maf = maf_path,
    genotype_matrix_market = matrix_market_path
  )
)

write_json(
  payload,
  path = r_summary_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

cat("Export complete.\n")
cat("  input:", input_rda, "\n")
cat("  outdir:", outdir, "\n")
cat("  n_samples:", nrow(genotype), " n_variants:", ncol(genotype), " nnz:", length(genotype@x), "\n")
cat("  matrix:", matrix_market_path, "\n")
cat("  summary:", r_summary_path, "\n")
