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

matrix_path <- get_arg("--matrix")
maf_path <- get_arg("--maf")
snploc_path <- get_arg("--snploc")
out_path <- get_arg("--out")
runs <- get_arg_int("--runs", 5)
warmup <- get_arg_int("--warmup", 1)

if (is.null(matrix_path) || is.null(maf_path) || is.null(snploc_path) || is.null(out_path)) {
  stop("Missing required args. Need --matrix --maf --snploc --out")
}
if (runs < 1) {
  stop("--runs must be >= 1")
}
if (warmup < 0) {
  stop("--warmup must be >= 0")
}

matrix_path <- normalizePath(matrix_path, mustWork = TRUE)
maf_path <- normalizePath(maf_path, mustWork = TRUE)
snploc_path <- normalizePath(snploc_path, mustWork = TRUE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
out_path <- normalizePath(out_path, mustWork = FALSE)

read_mm <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(readMM(gzfile(path)))
  }
  readMM(path)
}

run_once <- function() {
  t0 <- proc.time()[3]
  geno <- read_mm(matrix_path)
  maf <- read.csv(maf_path, stringsAsFactors = FALSE)
  snploc <- read.csv(snploc_path, stringsAsFactors = FALSE)

  maf_values <- as.numeric(maf$maf)
  maf_q <- stats::quantile(maf_values, probs = c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE)

  nnz <- if (methods::is(geno, "sparseMatrix")) length(geno@x) else sum(geno != 0)
  variant_match <- FALSE
  if ("variant_id" %in% names(maf) && "variant_id" %in% names(snploc)) {
    variant_match <- identical(as.character(maf$variant_id), as.character(snploc$variant_id))
  }

  elapsed <- as.numeric(proc.time()[3] - t0)
  list(
    elapsed_seconds = elapsed,
    summary = list(
      n_samples = as.integer(nrow(geno)),
      n_variants = as.integer(ncol(geno)),
      nnz = as.integer(nnz),
      genotype_sum = as.numeric(sum(geno@x)),
      maf_mean = as.numeric(mean(maf_values)),
      maf_summary = list(
        min = as.numeric(maf_q[[1]]),
        q01 = as.numeric(maf_q[[2]]),
        median = as.numeric(maf_q[[3]]),
        q99 = as.numeric(maf_q[[4]]),
        max = as.numeric(maf_q[[5]])
      ),
      snploc_rows = as.integer(nrow(snploc)),
      variant_id_match = variant_match
    )
  )
}

if (warmup > 0) {
  for (i in seq_len(warmup)) {
    invisible(run_once())
  }
}

timings <- numeric(0)
last_summary <- NULL
for (i in seq_len(runs)) {
  res <- run_once()
  timings <- c(timings, res$elapsed_seconds)
  last_summary <- res$summary
}

payload <- list(
  benchmark = list(
    language = "R",
    target = "matrix_market_load_and_summary",
    runs = as.integer(runs),
    warmup = as.integer(warmup),
    timings_seconds = as.list(timings),
    median_seconds = as.numeric(stats::median(timings)),
    mean_seconds = as.numeric(mean(timings)),
    min_seconds = as.numeric(min(timings)),
    max_seconds = as.numeric(max(timings))
  ),
  summary = last_summary
)

write_json(payload, path = out_path, pretty = TRUE, auto_unbox = TRUE, digits = 12)
cat("R benchmark complete:", out_path, "\n")
