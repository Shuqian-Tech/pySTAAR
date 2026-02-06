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
  kins = sim$kins_sparse,
  id = "id",
  family = binomial(link = "logit")
)

pvalues <- STAAR_Binary_SPA(
  genotype = sim$Geno,
  obj_nullmodel = obj_nullmodel,
  annotation_phred = sim$PHRED,
  rare_maf_cutoff = 0.05,
  SPA_p_filter = TRUE,
  p_filter_cutoff = 0.05
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
  "baselines/related_sparse_binary_spa_filter_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)

geno_flip <- matrix_flip(as.matrix(sim$Geno))
rv_label <- as.vector((geno_flip$MAF < 0.05) & (geno_flip$MAF > 0))
G <- geno_flip$Geno[, rv_label, drop = FALSE]
G_mat <- as.matrix(G)

Sigma_i <- obj_nullmodel$Sigma_i
Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
cov <- as.matrix(obj_nullmodel$cov)
Sigma_iG <- as.matrix(Sigma_i %*% G_mat)
tSigma_iX_G <- t(Sigma_iX) %*% G_mat
Cov_filter <- t(Sigma_iG) %*% G_mat - t(tSigma_iX_G) %*% cov %*% tSigma_iX_G

utils::write.csv(
  as.matrix(Cov_filter),
  "data/example_glmmkin_binary_spa_sparse_cov_filter.csv",
  row.names = FALSE
)
