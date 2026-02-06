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
  kins = sim$kins_dense
)

pop_group_df <- utils::read.csv("data/example_ai_pop_groups.csv", stringsAsFactors = FALSE)
pop_weights_1_1_df <- utils::read.csv("data/example_ai_pop_weights_1_1.csv", check.names = FALSE)
pop_weights_1_25_df <- utils::read.csv("data/example_ai_pop_weights_1_25.csv", check.names = FALSE)

if (!"pop_group" %in% colnames(pop_group_df)) {
  stop("AI pop group file must include 'pop_group' column.")
}
if (!"population" %in% colnames(pop_weights_1_1_df) || !"population" %in% colnames(pop_weights_1_25_df)) {
  stop("AI pop-weight files must include 'population' column.")
}

pop_groups <- as.character(pop_group_df$pop_group)
if (length(pop_groups) != nrow(sim$pheno_related)) {
  stop("AI pop groups length does not match related sample size.")
}

pop_weights_1_1 <- as.matrix(pop_weights_1_1_df[, setdiff(colnames(pop_weights_1_1_df), "population"), drop = FALSE])
pop_weights_1_25 <- as.matrix(pop_weights_1_25_df[, setdiff(colnames(pop_weights_1_25_df), "population"), drop = FALSE])

obj_nullmodel$pop.groups <- pop_groups
obj_nullmodel$pop_weights_1_1 <- pop_weights_1_1
obj_nullmodel$pop_weights_1_25 <- pop_weights_1_25

genotype_flip <- matrix_flip(as.matrix(sim$Geno))
MAF <- genotype_flip$MAF
RV_label <- as.vector((MAF < 0.05) & (MAF > 0))
Geno_rare <- genotype_flip$Geno[, RV_label, drop = FALSE]

pop_levels <- unique(pop_groups)
n_pop <- length(pop_levels)
indices <- vector("list", n_pop)
a_p <- matrix(0, nrow = 1, ncol = n_pop)
for (i in seq_len(n_pop)) {
  eth <- pop_levels[i]
  indices[[i]] <- which(pop_groups %in% eth)
  a_p[, i] <- mean(apply(as.matrix(Geno_rare[indices[[i]], , drop = FALSE]), 2, function(x) {
    min(mean(x) / 2, 1 - mean(x) / 2)
  }))
}
a_p <- ifelse(a_p > 0, dbeta(a_p, 1, 25), a_p)
B <- ncol(pop_weights_1_1)

for (b in seq_len(B)) {
  w_b_1 <- pop_weights_1_1[, b]
  w_b_2 <- as.vector(t(a_p %*% diag(pop_weights_1_25[, b])))

  Geno_rare_1 <- Geno_rare_2 <- Geno_rare
  for (i in seq_len(n_pop)) {
    Geno_rare_1[indices[[i]], ] <- w_b_1[i] * Geno_rare[indices[[i]], ]
    Geno_rare_2[indices[[i]], ] <- w_b_2[i] * Geno_rare[indices[[i]], ]
  }

  G1_sparse <- as(Geno_rare_1, "dgCMatrix")
  G2_sparse <- as(Geno_rare_2, "dgCMatrix")
  G1 <- as.matrix(G1_sparse)
  G2 <- as.matrix(G2_sparse)

  PG1 <- as.matrix(obj_nullmodel$P %*% G1_sparse)
  PG2 <- as.matrix(obj_nullmodel$P %*% G2_sparse)
  cov_1 <- t(PG1) %*% G1
  cov_2 <- t(PG2) %*% G2
  cov_1 <- (cov_1 + t(cov_1)) / 2
  cov_2 <- (cov_2 + t(cov_2)) / 2

  utils::write.csv(cov_1, sprintf("data/example_ai_cov_dense_s1_b%d.csv", b), row.names = FALSE)
  utils::write.csv(cov_2, sprintf("data/example_ai_cov_dense_s2_b%d.csv", b), row.names = FALSE)
}

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
  "baselines/ai_staar_related_dense_glmmkin_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
