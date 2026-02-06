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

set.seed(600)
pop_levels <- c("EUR", "AFR", "EAS")
n_pop <- length(pop_levels)
B <- 2

pop_groups <- sample(pop_levels, nrow(sim$pheno_unrelated), replace = TRUE, prob = c(0.5, 0.3, 0.2))
pop_groups[seq_len(n_pop)] <- pop_levels

normalize_cols <- function(mat) {
  sweep(mat, 2, colSums(mat), "/")
}
pop_weights_1_1 <- normalize_cols(matrix(runif(n_pop * B, min = 0.2, max = 1.0), nrow = n_pop, ncol = B))
pop_weights_1_25 <- normalize_cols(matrix(runif(n_pop * B, min = 0.2, max = 1.0), nrow = n_pop, ncol = B))

obj_nullmodel$pop.groups <- pop_groups
obj_nullmodel$pop_weights_1_1 <- pop_weights_1_1
obj_nullmodel$pop_weights_1_25 <- pop_weights_1_25

pvalues <- AI_STAAR(
  genotype = sim$Geno,
  obj_nullmodel = obj_nullmodel,
  annotation_phred = sim$PHRED,
  rare_maf_cutoff = 0.05,
  find_weight = FALSE
)

row_to_named_list <- function(df_row) {
  vec <- as.numeric(df_row)
  names(vec) <- colnames(df_row)
  as.list(vec)
}

sentinels <- list(
  num_variant = as.numeric(pvalues$num_variant[1]),
  cMAC = as.numeric(pvalues$cMAC[1]),
  results_STAAR_O = as.numeric(pvalues$results_STAAR_O[1]),
  results_ACAT_O = as.numeric(pvalues$results_ACAT_O[1]),
  results_STAAR_S_1_25 = row_to_named_list(pvalues$results_STAAR_S_1_25[1, , drop = FALSE]),
  results_STAAR_S_1_1 = row_to_named_list(pvalues$results_STAAR_S_1_1[1, , drop = FALSE]),
  results_STAAR_B_1_25 = row_to_named_list(pvalues$results_STAAR_B_1_25[1, , drop = FALSE]),
  results_STAAR_B_1_1 = row_to_named_list(pvalues$results_STAAR_B_1_1[1, , drop = FALSE]),
  results_STAAR_A_1_25 = row_to_named_list(pvalues$results_STAAR_A_1_25[1, , drop = FALSE]),
  results_STAAR_A_1_1 = row_to_named_list(pvalues$results_STAAR_A_1_1[1, , drop = FALSE])
)

jsonlite::write_json(
  sentinels,
  "baselines/ai_staar_unrelated_glm_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)

pop_group_df <- data.frame(pop_group = pop_groups)
utils::write.csv(pop_group_df, "data/example_ai_pop_groups.csv", row.names = FALSE)

pop_weights_1_1_df <- data.frame(
  population = pop_levels,
  pop_weights_1_1,
  check.names = FALSE
)
colnames(pop_weights_1_1_df) <- c("population", paste0("B", seq_len(B)))
utils::write.csv(pop_weights_1_1_df, "data/example_ai_pop_weights_1_1.csv", row.names = FALSE)

pop_weights_1_25_df <- data.frame(
  population = pop_levels,
  pop_weights_1_25,
  check.names = FALSE
)
colnames(pop_weights_1_25_df) <- c("population", paste0("B", seq_len(B)))
utils::write.csv(pop_weights_1_25_df, "data/example_ai_pop_weights_1_25.csv", row.names = FALSE)
