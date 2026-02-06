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
  family = "gaussian"
)

pvalues <- STAAR(
  genotype = sim$Geno,
  obj_nullmodel = obj_nullmodel,
  annotation_phred = sim$PHRED,
  rare_maf_cutoff = 0.05
)

row_to_named_list <- function(df_row) {
  vec <- as.numeric(df_row)
  names(vec) <- colnames(df_row)
  as.list(vec)
}

beta_vec <- if (!is.null(obj_nullmodel$beta)) as.numeric(obj_nullmodel$beta) else NULL
beta_list <- if (!is.null(beta_vec)) {
  names(beta_vec) <- names(obj_nullmodel$beta)
  as.list(beta_vec)
} else {
  NULL
}

sentinels <- list(
  num_variant = as.numeric(pvalues$num_variant[1]),
  results_STAAR_O = as.numeric(pvalues$results_STAAR_O[1]),
  results_STAAR_S_1_25 = row_to_named_list(pvalues$results_STAAR_S_1_25[1, , drop = FALSE]),
  results_STAAR_S_1_1 = row_to_named_list(pvalues$results_STAAR_S_1_1[1, , drop = FALSE]),
  results_STAAR_B_1_25 = row_to_named_list(pvalues$results_STAAR_B_1_25[1, , drop = FALSE]),
  results_STAAR_B_1_1 = row_to_named_list(pvalues$results_STAAR_B_1_1[1, , drop = FALSE]),
  results_STAAR_A_1_25 = row_to_named_list(pvalues$results_STAAR_A_1_25[1, , drop = FALSE]),
  results_STAAR_A_1_1 = row_to_named_list(pvalues$results_STAAR_A_1_1[1, , drop = FALSE]),
  nullmodel_beta = beta_list,
  nullmodel_dispersion = as.numeric(obj_nullmodel$dispersion)
)

jsonlite::write_json(
  sentinels,
  "baselines/unrelated_glm_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
