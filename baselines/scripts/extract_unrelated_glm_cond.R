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

# Condition on two very-rare variants from the same example genotype matrix.
genotype_adj <- as.matrix(sim$Geno[, c(1, 4), drop = FALSE])

pvalues <- STAAR_cond(
  genotype = sim$Geno,
  genotype_adj = genotype_adj,
  obj_nullmodel = obj_nullmodel,
  annotation_phred = sim$PHRED,
  rare_maf_cutoff = 0.05,
  method_cond = "optimal"
)

row_to_named_list <- function(df_row) {
  vec <- as.numeric(df_row)
  names(vec) <- colnames(df_row)
  as.list(vec)
}

sentinels <- list(
  num_variant = as.numeric(pvalues$num_variant[1]),
  cMAC = as.numeric(pvalues$cMAC[1]),
  results_STAAR_O_cond = as.numeric(pvalues$results_STAAR_O_cond[1]),
  results_ACAT_O_cond = as.numeric(pvalues$results_ACAT_O_cond[1]),
  results_STAAR_S_1_25_cond = row_to_named_list(pvalues$results_STAAR_S_1_25_cond[1, , drop = FALSE]),
  results_STAAR_S_1_1_cond = row_to_named_list(pvalues$results_STAAR_S_1_1_cond[1, , drop = FALSE]),
  results_STAAR_B_1_25_cond = row_to_named_list(pvalues$results_STAAR_B_1_25_cond[1, , drop = FALSE]),
  results_STAAR_B_1_1_cond = row_to_named_list(pvalues$results_STAAR_B_1_1_cond[1, , drop = FALSE]),
  results_STAAR_A_1_25_cond = row_to_named_list(pvalues$results_STAAR_A_1_25_cond[1, , drop = FALSE]),
  results_STAAR_A_1_1_cond = row_to_named_list(pvalues$results_STAAR_A_1_1_cond[1, , drop = FALSE])
)

jsonlite::write_json(
  sentinels,
  "baselines/unrelated_glm_cond_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
