options(mc.cores = 1)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("STAAR", quietly = TRUE)) {
  stop("STAAR package is required.")
}

library(STAAR)

sim <- readRDS("baselines/example_sim_data.rds")
pheno_bin <- sim$pheno_related
threshold <- as.numeric(stats::quantile(pheno_bin$Y, probs = 0.95, type = 7))
pheno_bin$Y <- as.integer(pheno_bin$Y > threshold)

obj_nullmodel <- fit_null_glmmkin_Binary_SPA(
  Y ~ X1 + X2,
  data = pheno_bin,
  kins = sim$kins_dense,
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

row_to_named_list <- function(df_row) {
  vec <- as.numeric(df_row)
  names(vec) <- colnames(df_row)
  as.list(vec)
}

sentinels <- list(
  num_variant = as.numeric(pvalues$num_variant[1]),
  cMAC = as.numeric(pvalues$cMAC[1]),
  case_count = as.numeric(sum(pheno_bin$Y)),
  results_STAAR_B = as.numeric(pvalues$results_STAAR_B[1]),
  results_STAAR_B_1_25 = row_to_named_list(pvalues$results_STAAR_B_1_25[1, , drop = FALSE]),
  results_STAAR_B_1_1 = row_to_named_list(pvalues$results_STAAR_B_1_1[1, , drop = FALSE])
)

jsonlite::write_json(
  sentinels,
  "baselines/related_dense_binary_spa_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)

utils::write.csv(as.matrix(obj_nullmodel$fitted.values), "data/example_glmmkin_binary_spa_dense_fitted.csv", row.names = FALSE)
utils::write.csv(as.matrix(obj_nullmodel$scaled.residuals), "data/example_glmmkin_binary_spa_dense_scaled_residuals.csv", row.names = FALSE)
utils::write.csv(as.matrix(obj_nullmodel$XW), "data/example_glmmkin_binary_spa_dense_XW.csv", row.names = FALSE)
utils::write.csv(as.matrix(obj_nullmodel$XXWX_inv), "data/example_glmmkin_binary_spa_dense_XXWX_inv.csv", row.names = FALSE)
