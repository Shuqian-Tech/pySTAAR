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

genotype_flip <- matrix_flip(as.matrix(sim$Geno))
maf <- genotype_flip$MAF
rv_label <- as.vector((maf < 0.05) & (maf > 0))
G <- as.matrix(genotype_flip$Geno[, rv_label, drop = FALSE])

residuals_phenotype <- obj_nullmodel$scaled.residuals
residuals_fit <- lm(residuals_phenotype ~ genotype_adj + obj_nullmodel$X - 1)
X_adj <- model.matrix(residuals_fit)

Sigma_i <- obj_nullmodel$Sigma_i
Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
cov_mat <- obj_nullmodel$cov
xtx_inv <- solve(t(X_adj) %*% X_adj)

Sigma_iG <- as.matrix(Sigma_i %*% G)
Sigma_iXG <- t(Sigma_iX) %*% G
PX_adj <- as.matrix(Sigma_i %*% X_adj - Sigma_iX %*% cov_mat %*% (t(Sigma_iX) %*% X_adj))

Cov_cond <- t(Sigma_iG) %*% G -
  t(Sigma_iXG) %*% cov_mat %*% Sigma_iXG -
  t(t(X_adj) %*% G) %*% xtx_inv %*% (t(PX_adj) %*% G) -
  t(t(PX_adj) %*% G) %*% xtx_inv %*% (t(X_adj) %*% G) +
  t(t(X_adj) %*% G) %*% xtx_inv %*% t(PX_adj) %*% X_adj %*% xtx_inv %*% (t(X_adj) %*% G)

utils::write.csv(
  as.matrix(Cov_cond),
  "data/example_glmmkin_cov_cond_sparse.csv",
  row.names = FALSE
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
  "baselines/related_sparse_glmmkin_cond_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
