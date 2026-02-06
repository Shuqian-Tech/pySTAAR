options(mc.cores = 1)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("STAAR", quietly = TRUE)) {
  stop("STAAR package is required.")
}

library(STAAR)

sim <- readRDS("baselines/example_sim_data.rds")
obj_nullmodel <- fit_null_glm(
  Y ~ X1 + X2,
  data = sim$pheno_unrelated,
  family = gaussian(link = "identity")
)

pop_group_df <- utils::read.csv("data/example_ai_pop_groups.csv", stringsAsFactors = FALSE)
pop_weights_1_1_df <- utils::read.csv("data/example_ai_pop_weights_1_1.csv", check.names = FALSE)
pop_weights_1_25_df <- utils::read.csv("data/example_ai_pop_weights_1_25.csv", check.names = FALSE)

if (!"pop_group" %in% colnames(pop_group_df)) {
  stop("AI pop group file must include 'pop_group' column.")
}
if (!"population" %in% colnames(pop_weights_1_1_df) || !"population" %in% colnames(pop_weights_1_25_df)) {
  stop("AI pop-weight files must include 'population' column.")
}

pop_groups <- as.character(pop_group_df$pop_group)
pop_levels <- unique(pop_groups)
pop_weights_1_1 <- as.matrix(pop_weights_1_1_df[, setdiff(colnames(pop_weights_1_1_df), "population"), drop = FALSE])
pop_weights_1_25 <- as.matrix(pop_weights_1_25_df[, setdiff(colnames(pop_weights_1_25_df), "population"), drop = FALSE])

obj_nullmodel$pop.groups <- pop_groups
obj_nullmodel$pop_weights_1_1 <- pop_weights_1_1
obj_nullmodel$pop_weights_1_25 <- pop_weights_1_25

pvalues <- AI_STAAR(
  genotype = sim$Geno,
  obj_nullmodel = obj_nullmodel,
  annotation_phred = sim$PHRED,
  rare_maf_cutoff = 0.05,
  find_weight = TRUE
)

flatten_weight_matrix <- function(weight_matrix, levels) {
  out <- list()
  for (i in seq_len(nrow(weight_matrix))) {
    for (b in seq_len(ncol(weight_matrix))) {
      out[[sprintf("%s.B%d", levels[i], b)]] <- as.numeric(weight_matrix[i, b])
    }
  }
  out
}

extract_metric_by_base <- function(weight_matrix, metric_name) {
  idx <- match(metric_name, rownames(weight_matrix))
  if (is.na(idx)) {
    stop(sprintf("Metric '%s' not found in results matrix.", metric_name))
  }
  out <- list()
  for (b in seq_len(ncol(weight_matrix))) {
    out[[sprintf("B%d", b)]] <- as.numeric(weight_matrix[idx, b])
  }
  out
}

sentinels <- list(
  num_variant = as.numeric(pvalues$num_variant[1]),
  cMAC = as.numeric(pvalues$cMAC[1]),
  results_STAAR_O = as.numeric(pvalues$results_STAAR_O[1]),
  results_ACAT_O = as.numeric(pvalues$results_ACAT_O[1]),
  weight_all_1 = flatten_weight_matrix(pvalues$weight_all_1, pop_levels),
  weight_all_2 = flatten_weight_matrix(pvalues$weight_all_2, pop_levels),
  results_weight_staar_o = extract_metric_by_base(pvalues$results_weight, "results_STAAR_O"),
  results_weight1_staar_o = extract_metric_by_base(pvalues$results_weight1, "results_STAAR_O"),
  results_weight2_staar_o = extract_metric_by_base(pvalues$results_weight2, "results_STAAR_O"),
  results_weight_staar_s_1_25 = extract_metric_by_base(
    pvalues$results_weight,
    "results_STAAR_S_1_25.STAAR-S(1,25)"
  ),
  results_weight1_staar_s_1_25 = extract_metric_by_base(
    pvalues$results_weight1,
    "results_STAAR_S_1_25.STAAR-S(1,25)"
  ),
  results_weight2_staar_s_1_25 = extract_metric_by_base(
    pvalues$results_weight2,
    "results_STAAR_S_1_25.STAAR-S(1,25)"
  )
)

jsonlite::write_json(
  sentinels,
  "baselines/ai_staar_unrelated_glm_find_weight_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
