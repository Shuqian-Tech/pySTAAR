#!/usr/bin/env Rscript

options(mc.cores = 1)

suppressPackageStartupMessages({
  library(STAAR)
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

parse_int_list <- function(value, default) {
  if (is.null(value) || nchar(trimws(value)) == 0) {
    return(default)
  }
  pieces <- unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE)
  as.integer(trimws(pieces))
}

sim_rds <- get_arg("--sim-rds")
out_results <- get_arg("--out-results")
out_benchmark <- get_arg("--out-benchmark")
runs <- get_arg_int("--runs", 3L)
warmup <- get_arg_int("--warmup", 1L)
seed <- get_arg_int("--seed", 600L)
rare_maf_cutoff <- get_arg_num("--rare-maf-cutoff", 0.05)
adj_variants <- parse_int_list(get_arg("--adj-variants", "1"), c(1L))

if (is.null(sim_rds) || is.null(out_results) || is.null(out_benchmark)) {
  stop("Required args: --sim-rds --out-results --out-benchmark")
}
if (runs < 1) {
  stop("--runs must be >= 1")
}
if (warmup < 0) {
  stop("--warmup must be >= 0")
}
if (length(adj_variants) < 1) {
  stop("--adj-variants must contain at least one index")
}

sim_rds <- normalizePath(sim_rds, mustWork = TRUE)
dir.create(dirname(out_results), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_benchmark), recursive = TRUE, showWarnings = FALSE)
out_results <- normalizePath(out_results, mustWork = FALSE)
out_benchmark <- normalizePath(out_benchmark, mustWork = FALSE)

sim <- readRDS(sim_rds)
n_variants <- ncol(sim$Geno)
if (any(adj_variants < 1) || any(adj_variants > n_variants)) {
  stop("adj_variants must be within [1, n_variants].")
}
adj_variants_requested <- adj_variants

geno_dense_for_adj <- as.matrix(sim$Geno)
var_positive <- which(apply(geno_dense_for_adj, 2, function(v) stats::var(as.numeric(v)) > 0))
if (length(var_positive) < 1) {
  stop("Need at least one non-constant variant for conditional workflows.")
}
adj_variants <- unique(adj_variants[adj_variants %in% var_positive])
target_adj_len <- max(1L, length(adj_variants_requested))
if (length(adj_variants) < target_adj_len) {
  extras <- var_positive[!var_positive %in% adj_variants]
  need <- target_adj_len - length(adj_variants)
  adj_variants <- c(adj_variants, extras[seq_len(min(need, length(extras)))])
}
adj_variants <- as.integer(adj_variants[seq_len(min(length(adj_variants), target_adj_len))])

extract_staar <- function(pvalues) {
  list(
    num_variant = as.numeric(pvalues$num_variant[1]),
    results_STAAR_O = as.numeric(pvalues$results_STAAR_O[1])
  )
}

extract_cond <- function(pvalues) {
  list(
    num_variant = as.numeric(pvalues$num_variant[1]),
    cMAC = as.numeric(pvalues$cMAC[1]),
    results_STAAR_O_cond = as.numeric(pvalues$results_STAAR_O_cond[1]),
    results_ACAT_O_cond = as.numeric(pvalues$results_ACAT_O_cond[1])
  )
}

extract_binary <- function(pvalues, case_count) {
  list(
    num_variant = as.numeric(pvalues$num_variant[1]),
    cMAC = as.numeric(pvalues$cMAC[1]),
    case_count = as.numeric(case_count),
    results_STAAR_B = as.numeric(pvalues$results_STAAR_B[1])
  )
}

quiet_eval <- function(expr) {
  capture.output(result <- eval.parent(substitute(expr)))
  result
}

run_suite <- function() {
  set.seed(seed)
  timings <- list()
  results <- list()
  suite_t0 <- proc.time()[3]

  # STAAR: unrelated
  t0 <- proc.time()[3]
  obj_unrelated <- fit_null_glm(
    Y ~ X1 + X2,
    data = sim$pheno_unrelated,
    family = "gaussian"
  )
  p_unrelated <- STAAR(
    genotype = sim$Geno,
    obj_nullmodel = obj_unrelated,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff
  )
  timings$staar_unrelated_glm <- as.numeric(proc.time()[3] - t0)
  results$staar_unrelated_glm <- extract_staar(p_unrelated)

  # STAAR: related sparse
  t0 <- proc.time()[3]
  obj_related_sparse <- quiet_eval(fit_null_glmmkin(
    Y ~ X1 + X2,
    data = sim$pheno_related,
    family = gaussian(link = "identity"),
    id = "id",
    kins = sim$kins_sparse
  ))
  p_related_sparse <- STAAR(
    genotype = sim$Geno,
    obj_nullmodel = obj_related_sparse,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff
  )
  timings$staar_related_sparse_glmmkin <- as.numeric(proc.time()[3] - t0)
  results$staar_related_sparse_glmmkin <- extract_staar(p_related_sparse)

  # STAAR: related dense
  t0 <- proc.time()[3]
  obj_related_dense <- quiet_eval(fit_null_glmmkin(
    Y ~ X1 + X2,
    data = sim$pheno_related,
    family = gaussian(link = "identity"),
    id = "id",
    kins = sim$kins_dense
  ))
  p_related_dense <- STAAR(
    genotype = sim$Geno,
    obj_nullmodel = obj_related_dense,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff
  )
  timings$staar_related_dense_glmmkin <- as.numeric(proc.time()[3] - t0)
  results$staar_related_dense_glmmkin <- extract_staar(p_related_dense)

  genotype_adj <- as.matrix(sim$Geno[, adj_variants, drop = FALSE])

  # Conditional: unrelated
  t0 <- proc.time()[3]
  p_unrelated_cond <- STAAR_cond(
    genotype = sim$Geno,
    genotype_adj = genotype_adj,
    obj_nullmodel = obj_unrelated,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    method_cond = "optimal"
  )
  timings$staar_unrelated_glm_cond <- as.numeric(proc.time()[3] - t0)
  results$staar_unrelated_glm_cond <- extract_cond(p_unrelated_cond)

  # Conditional: related sparse
  t0 <- proc.time()[3]
  p_related_sparse_cond <- STAAR_cond(
    genotype = sim$Geno,
    genotype_adj = genotype_adj,
    obj_nullmodel = obj_related_sparse,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    method_cond = "optimal"
  )
  timings$staar_related_sparse_glmmkin_cond <- as.numeric(proc.time()[3] - t0)
  results$staar_related_sparse_glmmkin_cond <- extract_cond(p_related_sparse_cond)

  # Conditional: related dense
  t0 <- proc.time()[3]
  p_related_dense_cond <- STAAR_cond(
    genotype = sim$Geno,
    genotype_adj = genotype_adj,
    obj_nullmodel = obj_related_dense,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    method_cond = "optimal"
  )
  timings$staar_related_dense_glmmkin_cond <- as.numeric(proc.time()[3] - t0)
  results$staar_related_dense_glmmkin_cond <- extract_cond(p_related_dense_cond)

  # Binary SPA: unrelated
  t0 <- proc.time()[3]
  pheno_unrelated_bin <- sim$pheno_unrelated
  threshold_unrelated <- as.numeric(stats::quantile(pheno_unrelated_bin$Y, probs = 0.95, type = 7))
  pheno_unrelated_bin$Y <- as.integer(pheno_unrelated_bin$Y > threshold_unrelated)
  case_count_unrelated <- sum(pheno_unrelated_bin$Y)
  obj_unrelated_bin <- fit_null_glm_Binary_SPA(
    Y ~ X1 + X2,
    data = pheno_unrelated_bin,
    family = binomial(link = "logit")
  )
  p_unrelated_bin <- STAAR_Binary_SPA(
    genotype = sim$Geno,
    obj_nullmodel = obj_unrelated_bin,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    SPA_p_filter = FALSE
  )
  timings$staar_unrelated_binary_spa <- as.numeric(proc.time()[3] - t0)
  results$staar_unrelated_binary_spa <- extract_binary(p_unrelated_bin, case_count_unrelated)

  # Binary SPA: related sparse
  t0 <- proc.time()[3]
  pheno_related_bin <- sim$pheno_related
  threshold_related <- as.numeric(stats::quantile(pheno_related_bin$Y, probs = 0.95, type = 7))
  pheno_related_bin$Y <- as.integer(pheno_related_bin$Y > threshold_related)
  case_count_related <- sum(pheno_related_bin$Y)
  obj_related_sparse_bin <- quiet_eval(fit_null_glmmkin_Binary_SPA(
    Y ~ X1 + X2,
    data = pheno_related_bin,
    kins = sim$kins_sparse,
    id = "id",
    family = binomial(link = "logit")
  ))
  p_related_sparse_bin <- STAAR_Binary_SPA(
    genotype = sim$Geno,
    obj_nullmodel = obj_related_sparse_bin,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    SPA_p_filter = FALSE
  )
  timings$staar_related_sparse_binary_spa <- as.numeric(proc.time()[3] - t0)
  results$staar_related_sparse_binary_spa <- extract_binary(p_related_sparse_bin, case_count_related)

  # Binary SPA: related dense
  t0 <- proc.time()[3]
  obj_related_dense_bin <- quiet_eval(fit_null_glmmkin_Binary_SPA(
    Y ~ X1 + X2,
    data = pheno_related_bin,
    kins = sim$kins_dense,
    id = "id",
    family = binomial(link = "logit")
  ))
  p_related_dense_bin <- STAAR_Binary_SPA(
    genotype = sim$Geno,
    obj_nullmodel = obj_related_dense_bin,
    annotation_phred = sim$PHRED,
    rare_maf_cutoff = rare_maf_cutoff,
    SPA_p_filter = FALSE
  )
  timings$staar_related_dense_binary_spa <- as.numeric(proc.time()[3] - t0)
  results$staar_related_dense_binary_spa <- extract_binary(p_related_dense_bin, case_count_related)

  total_elapsed <- as.numeric(proc.time()[3] - suite_t0)
  list(results = results, timings = timings, total_elapsed_seconds = total_elapsed)
}

if (warmup > 0) {
  for (i in seq_len(warmup)) {
    invisible(run_suite())
  }
}

timing_rows <- list()
total_timings <- numeric(0)
final_results <- NULL
workflow_names <- NULL
for (i in seq_len(runs)) {
  run_out <- run_suite()
  final_results <- run_out$results
  timing_rows[[i]] <- run_out$timings
  total_timings <- c(total_timings, run_out$total_elapsed_seconds)
  if (is.null(workflow_names)) {
    workflow_names <- names(run_out$timings)
  }
}

benchmark_workflows <- list()
for (name in workflow_names) {
  values <- vapply(timing_rows, function(row) as.numeric(row[[name]]), numeric(1))
  benchmark_workflows[[name]] <- list(
    timings_seconds = as.list(values),
    median_seconds = as.numeric(stats::median(values)),
    mean_seconds = as.numeric(mean(values)),
    min_seconds = as.numeric(min(values)),
    max_seconds = as.numeric(max(values))
  )
}

benchmark_payload <- list(
  benchmark = list(
    language = "R",
    suite = "full_workflow_simulated",
    runs = as.integer(runs),
    warmup = as.integer(warmup),
    rare_maf_cutoff = as.numeric(rare_maf_cutoff),
    adj_variants_requested = as.list(as.integer(adj_variants_requested)),
    adj_variants_used = as.list(as.integer(adj_variants)),
    total_timings_seconds = as.list(total_timings),
    total_median_seconds = as.numeric(stats::median(total_timings)),
    total_mean_seconds = as.numeric(mean(total_timings)),
    workflows = benchmark_workflows
  )
)

write_json(final_results, path = out_results, pretty = TRUE, auto_unbox = TRUE, digits = 12)
write_json(benchmark_payload, path = out_benchmark, pretty = TRUE, auto_unbox = TRUE, digits = 12)

cat("R full-workflow benchmark complete.\n")
cat("  results:", out_results, "\n")
cat("  benchmark:", out_benchmark, "\n")
