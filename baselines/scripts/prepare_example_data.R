options(mc.cores = 1)
set.seed(600)

required_pkgs <- c("Matrix", "STAAR", "kinship2", "MASS", "matrixStats", "rje")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste("Missing required packages:", paste(missing_pkgs, collapse = ", ")))
}

library(Matrix)
library(STAAR)
library(kinship2)
library(MASS)
library(matrixStats)
library(rje)

# Load example dataset
suppressMessages(data("example", package = "STAAR"))

genotype <- example$genotype
maf <- example$maf
snploc <- example$snploc

# Kinship matrix setup
grid <- 1
Npercell <- 10000
ndiv <- 1
vfam <- 0.5
N <- round(grid * grid * Npercell / ndiv)

unitmat <- matrix(0.5, 4, 4)
diag(unitmat) <- 1
unitmat[1, 2] <- unitmat[2, 1] <- 0

ped <- data.frame(
  famid = rep(as.integer(1:2500), each = 4),
  id = as.integer(1:10000L),
  fa = rep(0, 10000),
  mo = rep(0, 10000)
)
for (i in 1:2500) {
  ped$fa[4 * i - (0:1)] <- ped$id[4 * i - 3]
  ped$mo[4 * i - (0:1)] <- ped$id[4 * i - 2]
}
kins <- makekinship(ped$famid, ped$id, ped$fa, ped$mo)

# Covariates
X1 <- rnorm(N)
X2 <- rbinom(N, 1, 0.5)
eps <- rnorm(N)

# Functional annotations
numVar <- dim(snploc)[1]
Z1 <- rnorm(numVar)
Z2 <- rnorm(numVar)
Z3 <- rnorm(numVar)
Z4 <- rnorm(numVar)
Z5 <- rnorm(numVar)
Z6 <- rnorm(numVar)
Z7 <- rnorm(numVar)
Z8 <- rnorm(numVar)
Z9 <- rnorm(numVar)
Z10 <- rnorm(numVar)
Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10)

rank <- colRanks(Z, preserveShape = TRUE)
PHRED <- -10 * log10(1 - rank / dim(rank)[1])

# Signal region selection
maxLength <- 1000000
sigLength <- 5000
startloc <- sample(1:(maxLength - sigLength + 1), 1)
endloc <- startloc + sigLength - 1
snplist <- which(snploc$CHROM_POS >= startloc & snploc$CHROM_POS <= endloc)
snpRegion <- snploc[snplist, ]
numSNPs <- dim(snpRegion)[1]
Geno <- genotype[, snplist]

# Simulate causal variants
b0 <- rje::logit(0.015)
b <- rep(log(5), 10)
causalprob <- apply(Z[snplist, ], 1, function(z) {
  ind <- sample(1:10, 5)
  rje::expit(b0 + b[ind] %*% z[ind])
})
isCausal <- rbinom(numSNPs, 1, causalprob)

# Effect sizes and phenotype generation
alpha0 <- 0
alpha1 <- 0.5
alpha2 <- 0.5
c0 <- 0.13
beta <- -c0 * log10(maf[snplist][which(isCausal == 1)])

Geno_mat <- as.matrix(Geno)
causal_idx <- which(isCausal == 1)

Y_unrelated <- alpha0 + alpha1 * X1 + alpha2 * X2 +
  Geno_mat[, causal_idx, drop = FALSE] %*% beta + eps

randomfam <- mvrnorm(N / 4, rep(0, 4), vfam * unitmat)
randomfam <- as.vector(t(randomfam))
id0 <- 1:N

Y_related <- alpha0 + alpha1 * X1 + alpha2 * X2 +
  Geno_mat[, causal_idx, drop = FALSE] %*% beta + randomfam + eps

pheno_unrelated <- data.frame(Y = as.vector(Y_unrelated), X1 = X1, X2 = X2)
pheno_related <- data.frame(Y = as.vector(Y_related), X1 = X1, X2 = X2, id = id0)

saveRDS(
  list(
    Geno = Geno,
    PHRED = PHRED[snplist, , drop = FALSE],
    pheno_unrelated = pheno_unrelated,
    pheno_related = pheno_related,
    kins_sparse = kins,
    kins_dense = as.matrix(kins),
    snplist = snplist,
    startloc = startloc,
    endloc = endloc,
    numSNPs = numSNPs,
    params = list(
      alpha0 = alpha0,
      alpha1 = alpha1,
      alpha2 = alpha2,
      c0 = c0,
      beta = beta,
      causal_variant_count = length(causal_idx),
      vfam = vfam,
      N = N,
      seed = 600
    )
  ),
  file = "baselines/example_sim_data.rds"
)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}

metadata <- list(
  seed = 600,
  N = N,
  startloc = startloc,
  endloc = endloc,
  numSNPs = numSNPs,
  causal_variant_count = length(causal_idx),
  snplist_length = length(snplist)
)
jsonlite::write_json(
  metadata,
  "baselines/example_sim_metadata.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
