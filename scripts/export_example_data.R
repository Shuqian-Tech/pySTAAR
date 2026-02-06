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

G <- as.matrix(sim$Geno)
G <- matrix_flip(G)
MAF <- G$MAF
RV_label <- (MAF < 0.05) & (MAF > 0)
G <- G$Geno[, RV_label]
G <- as(G, "dgCMatrix")

Sigma_i <- obj$Sigma_i
Sigma_iX <- as.matrix(obj$Sigma_iX)
cov <- obj$cov
tSigma_iX_G <- t(Sigma_iX) %*% G
Sigma_iG <- Sigma_i %*% G
Cov <- t(Sigma_iG) %*% G - t(tSigma_iX_G) %*% cov %*% tSigma_iX_G

write.csv(obj$scaled.residuals, "data/example_glmmkin_scaled_residuals.csv", row.names = FALSE)
write.csv(as.matrix(Cov), "data/example_glmmkin_cov.csv", row.names = FALSE)
