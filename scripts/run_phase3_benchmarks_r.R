#!/usr/bin/env Rscript

options(mc.cores = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    return(default)
  }
  if (idx == length(args)) {
    stop(sprintf("Missing value for %s", flag))
  }
  args[[idx + 1]]
}

as_int <- function(x, name) {
  val <- as.integer(x)
  if (is.na(val)) {
    stop(sprintf("Invalid integer for %s: %s", name, x))
  }
  val
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Unable to determine script path from commandArgs.")
}
script_path <- sub("^--file=", "", script_arg[[1]])
script_dir <- dirname(normalizePath(script_path))
root <- normalizePath(file.path(script_dir, ".."))

warmup_runs <- as_int(arg_value("--warmup-runs", "1"), "warmup-runs")
measured_runs <- as_int(arg_value("--measured-runs", "5"), "measured-runs")
if (warmup_runs < 0) {
  stop("warmup-runs must be >= 0")
}
if (measured_runs < 1) {
  stop("measured-runs must be >= 1")
}

raw_output <- arg_value(
  "--raw-output",
  file.path(root, "benchmarks", "phase3_baseline_r_raw.csv")
)
summary_output <- arg_value(
  "--summary-output",
  file.path(root, "benchmarks", "phase3_baseline_r_summary.csv")
)
meta_output <- arg_value(
  "--meta-output",
  file.path(root, "benchmarks", "phase3_baseline_r_meta.json")
)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite package is required.")
}
if (!requireNamespace("STAAR", quietly = TRUE)) {
  stop("STAAR package is required.")
}

library(STAAR)

sim <- readRDS(file.path(root, "baselines", "example_sim_data.rds"))

load_ai_metadata <- function(root_dir) {
  pop_group_df <- utils::read.csv(
    file.path(root_dir, "data", "example_ai_pop_groups.csv"),
    stringsAsFactors = FALSE
  )
  pop_weights_1_1_df <- utils::read.csv(
    file.path(root_dir, "data", "example_ai_pop_weights_1_1.csv"),
    check.names = FALSE
  )
  pop_weights_1_25_df <- utils::read.csv(
    file.path(root_dir, "data", "example_ai_pop_weights_1_25.csv"),
    check.names = FALSE
  )

  if (!"pop_group" %in% colnames(pop_group_df)) {
    stop("example_ai_pop_groups.csv must include 'pop_group'.")
  }
  if (!"population" %in% colnames(pop_weights_1_1_df) ||
      !"population" %in% colnames(pop_weights_1_25_df)) {
    stop("AI weight files must include 'population' column.")
  }

  pop_groups <- as.character(pop_group_df$pop_group)
  pop_levels <- unique(pop_groups)
  w11 <- pop_weights_1_1_df[match(pop_levels, pop_weights_1_1_df$population), , drop = FALSE]
  w125 <- pop_weights_1_25_df[match(pop_levels, pop_weights_1_25_df$population), , drop = FALSE]

  list(
    pop_groups = pop_groups,
    pop_weights_1_1 = as.matrix(w11[, setdiff(colnames(w11), "population"), drop = FALSE]),
    pop_weights_1_25 = as.matrix(w125[, setdiff(colnames(w125), "population"), drop = FALSE])
  )
}

ai_meta_unrelated <- load_ai_metadata(root)
ai_meta_related <- load_ai_metadata(root)

scenario_fns <- list(
  staar_unrelated_glm = function() {
    obj_nullmodel <- fit_null_glm(
      Y ~ X1 + X2,
      data = sim$pheno_unrelated,
      family = "gaussian"
    )
    pvalues <- STAAR(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05
    )
    as.numeric(pvalues$results_STAAR_O[1])
  },
  staar_related_sparse_glmmkin_pure = function() {
    obj_nullmodel <- fit_null_glmmkin(
      Y ~ X1 + X2,
      data = sim$pheno_related,
      family = gaussian(link = "identity"),
      id = "id",
      kins = sim$kins_sparse
    )
    pvalues <- STAAR(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05
    )
    as.numeric(pvalues$results_STAAR_O[1])
  },
  staar_unrelated_binary_spa = function() {
    pheno_bin <- sim$pheno_unrelated
    threshold <- as.numeric(stats::quantile(pheno_bin$Y, probs = 0.95, type = 7))
    pheno_bin$Y <- as.integer(pheno_bin$Y > threshold)
    obj_nullmodel <- fit_null_glm_Binary_SPA(
      Y ~ X1 + X2,
      data = pheno_bin,
      family = binomial(link = "logit")
    )
    pvalues <- STAAR_Binary_SPA(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05,
      SPA_p_filter = FALSE
    )
    as.numeric(pvalues$results_STAAR_B[1])
  },
  staar_related_sparse_binary_spa_pure = function() {
    pheno_bin <- sim$pheno_related
    threshold <- as.numeric(stats::quantile(pheno_bin$Y, probs = 0.95, type = 7))
    pheno_bin$Y <- as.integer(pheno_bin$Y > threshold)
    obj_nullmodel <- fit_null_glmmkin_Binary_SPA(
      Y ~ X1 + X2,
      data = pheno_bin,
      kins = sim$kins_sparse,
      id = "id",
      family = binomial(link = "logit")
    )
    pvalues <- STAAR_Binary_SPA(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05,
      SPA_p_filter = FALSE
    )
    as.numeric(pvalues$results_STAAR_B[1])
  },
  staar_unrelated_glm_cond = function() {
    obj_nullmodel <- fit_null_glm(
      Y ~ X1 + X2,
      data = sim$pheno_unrelated,
      family = "gaussian"
    )
    genotype_adj <- as.matrix(sim$Geno[, c(1, 4), drop = FALSE])
    pvalues <- STAAR_cond(
      genotype = sim$Geno,
      genotype_adj = genotype_adj,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05,
      method_cond = "optimal"
    )
    as.numeric(pvalues$results_STAAR_O_cond[1])
  },
  indiv_score_unrelated_glm = function() {
    obj_nullmodel <- fit_null_glm(
      Y ~ X1 + X2,
      data = sim$pheno_unrelated,
      family = gaussian(link = "identity")
    )
    res <- Indiv_Score_Test_Region(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      rare_maf_cutoff = 0.05
    )
    as.numeric(min(res$pvalue, na.rm = TRUE))
  },
  ai_staar_unrelated_glm = function() {
    obj_nullmodel <- fit_null_glm(
      Y ~ X1 + X2,
      data = sim$pheno_unrelated,
      family = gaussian(link = "identity")
    )
    obj_nullmodel$pop.groups <- ai_meta_unrelated$pop_groups
    obj_nullmodel$pop_weights_1_1 <- ai_meta_unrelated$pop_weights_1_1
    obj_nullmodel$pop_weights_1_25 <- ai_meta_unrelated$pop_weights_1_25
    pvalues <- AI_STAAR(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05,
      find_weight = FALSE
    )
    as.numeric(pvalues$results_STAAR_O[1])
  },
  ai_staar_related_sparse_glmmkin_find_weight_pure = function() {
    obj_nullmodel <- fit_null_glmmkin(
      Y ~ X1 + X2,
      data = sim$pheno_related,
      family = gaussian(link = "identity"),
      id = "id",
      kins = sim$kins_sparse
    )
    obj_nullmodel$pop.groups <- ai_meta_related$pop_groups
    obj_nullmodel$pop_weights_1_1 <- ai_meta_related$pop_weights_1_1
    obj_nullmodel$pop_weights_1_25 <- ai_meta_related$pop_weights_1_25
    pvalues <- AI_STAAR(
      genotype = sim$Geno,
      obj_nullmodel = obj_nullmodel,
      annotation_phred = sim$PHRED,
      rare_maf_cutoff = 0.05,
      find_weight = TRUE
    )
    as.numeric(pvalues$results_STAAR_O[1])
  }
)

raw_rows <- data.frame(
  scenario_id = character(),
  run_type = character(),
  run_index = integer(),
  seconds = double(),
  sentinel_value = double(),
  stringsAsFactors = FALSE
)

summary_rows <- data.frame(
  scenario_id = character(),
  warmup_runs = integer(),
  measured_runs = integer(),
  median_seconds = double(),
  mean_seconds = double(),
  min_seconds = double(),
  max_seconds = double(),
  std_seconds = double(),
  relative_speedup_vs_baseline = double(),
  stringsAsFactors = FALSE
)

sentinel_checks <- list()

for (scenario_id in names(scenario_fns)) {
  fn <- scenario_fns[[scenario_id]]
  for (warm_idx in seq_len(warmup_runs)) {
    start <- proc.time()[["elapsed"]]
    sentinel_value <- fn()
    elapsed <- proc.time()[["elapsed"]] - start
    raw_rows <- rbind(
      raw_rows,
      data.frame(
        scenario_id = scenario_id,
        run_type = "warmup",
        run_index = warm_idx,
        seconds = elapsed,
        sentinel_value = sentinel_value,
        stringsAsFactors = FALSE
      )
    )
    sentinel_checks[[paste0(scenario_id, ".warmup")]] <- sentinel_value
  }

  measured <- numeric()
  for (run_idx in seq_len(measured_runs)) {
    start <- proc.time()[["elapsed"]]
    sentinel_value <- fn()
    elapsed <- proc.time()[["elapsed"]] - start
    measured <- c(measured, elapsed)
    raw_rows <- rbind(
      raw_rows,
      data.frame(
        scenario_id = scenario_id,
        run_type = "measured",
        run_index = run_idx,
        seconds = elapsed,
        sentinel_value = sentinel_value,
        stringsAsFactors = FALSE
      )
    )
    sentinel_checks[[paste0(scenario_id, ".run", run_idx)]] <- sentinel_value
  }

  summary_rows <- rbind(
    summary_rows,
    data.frame(
      scenario_id = scenario_id,
      warmup_runs = warmup_runs,
      measured_runs = measured_runs,
      median_seconds = as.numeric(stats::median(measured)),
      mean_seconds = as.numeric(mean(measured)),
      min_seconds = as.numeric(min(measured)),
      max_seconds = as.numeric(max(measured)),
      std_seconds = if (length(measured) > 1) as.numeric(stats::sd(measured)) else 0.0,
      relative_speedup_vs_baseline = 1.0,
      stringsAsFactors = FALSE
    )
  )
}

dir.create(dirname(raw_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(meta_output), recursive = TRUE, showWarnings = FALSE)

utils::write.csv(raw_rows, raw_output, row.names = FALSE)
utils::write.csv(summary_rows, summary_output, row.names = FALSE)

generated_utc <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
meta <- list(
  generated_utc = generated_utc,
  warmup_runs = warmup_runs,
  measured_runs = measured_runs,
  dataset_fingerprint_file = "baselines/example_fingerprint.json",
  reference_backend_file = "reports/reference_backend.md",
  r_environment_file = "baselines/r_environment.json",
  r_version = as.character(getRversion()),
  staar_version = as.character(utils::packageVersion("STAAR")),
  platform = as.list(Sys.info()),
  scenarios = names(scenario_fns),
  sentinel_checks = sentinel_checks
)
jsonlite::write_json(meta, meta_output, pretty = TRUE, auto_unbox = TRUE)
