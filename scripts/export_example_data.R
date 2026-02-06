library(Matrix)
library(STAAR)

sim <- readRDS("baselines/example_sim_data.rds")

if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

writeMM(sim$Geno, "data/example_geno.mtx")
writeMM(sim$kins_sparse, "data/example_kins_sparse.mtx")
writeMM(Matrix(sim$kins_dense, sparse = TRUE), "data/example_kins_dense.mtx")

write.csv(sim$PHRED, "data/example_phred.csv", row.names = FALSE)
write.csv(sim$pheno_unrelated, "data/example_pheno_unrelated.csv", row.names = FALSE)
write.csv(sim$pheno_related, "data/example_pheno_related.csv", row.names = FALSE)

# Export related null-model artifacts for parity (scaled residuals + Cov)
obj <- fit_null_glmmkin(
  Y ~ X1 + X2,
  data = sim$pheno_related,
  family = gaussian(link = "identity"),
  id = "id",
  kins = sim$kins_sparse
)

compute_cov_for_cutoff <- function(genotype, obj_nullmodel, rare_maf_cutoff) {
  geno_flip <- matrix_flip(as.matrix(genotype))
  maf <- geno_flip$MAF
  rv_label <- (maf < rare_maf_cutoff) & (maf > 0)
  G <- as(geno_flip$Geno[, rv_label], "dgCMatrix")

  Sigma_i <- obj_nullmodel$Sigma_i
  Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
  cov <- obj_nullmodel$cov
  tSigma_iX_G <- t(Sigma_iX) %*% G
  Sigma_iG <- Sigma_i %*% G
  t(Sigma_iG) %*% G - t(tSigma_iX_G) %*% cov %*% tSigma_iX_G
}

Cov <- compute_cov_for_cutoff(sim$Geno, obj, rare_maf_cutoff = 0.05)
Cov_rare_maf_0_01 <- compute_cov_for_cutoff(sim$Geno, obj, rare_maf_cutoff = 0.01)

write.csv(obj$scaled.residuals, "data/example_glmmkin_scaled_residuals.csv", row.names = FALSE)
write.csv(as.matrix(Cov), "data/example_glmmkin_cov.csv", row.names = FALSE)
write.csv(as.matrix(Cov_rare_maf_0_01), "data/example_glmmkin_cov_rare_maf_0_01.csv", row.names = FALSE)
