options(mc.cores = 1)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("STAAR", quietly = TRUE)) {
  stop("STAAR package is required.")
}

library(STAAR)

sim <- readRDS("baselines/example_sim_data.rds")

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

theta_vec <- if (!is.null(obj_nullmodel$theta)) as.numeric(obj_nullmodel$theta) else NULL
theta_list <- if (!is.null(theta_vec)) {
  names(theta_vec) <- names(obj_nullmodel$theta)
  as.list(theta_vec)
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
  nullmodel_theta = theta_list
)

jsonlite::write_json(
  sentinels,
  "baselines/related_sparse_glmmkin_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
