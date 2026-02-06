options(mc.cores = 1)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("STAAR", quietly = TRUE)) {
  stop("STAAR package is required.")
}

library(STAAR)

sample_indices <- c(2, 3, 11, 51, 101, 161)
adj_variant_indices <- c(1, 4)

sample_map <- function(values, idx) {
  vec <- as.numeric(values[idx])
  names(vec) <- paste0("v", idx)
  as.list(vec)
}

build_sentinels <- function(results, score_col, se_col, pvalue_col, score_key, se_key, pvalue_key) {
  tested <- which(!is.na(results[[pvalue_col]]))
  pvals <- as.numeric(results[[pvalue_col]][tested])
  top_idx <- tested[which.min(pvals)]

  out <- list(
    num_variant = as.numeric(length(tested)),
    num_tested = as.numeric(length(tested)),
    pvalue_min = as.numeric(min(pvals)),
    pvalue_median = as.numeric(stats::median(pvals)),
    top_variant_index = as.numeric(top_idx),
    top_pvalue = as.numeric(results[[pvalue_col]][top_idx])
  )
  out[[score_key]] <- sample_map(results[[score_col]], sample_indices)
  out[[se_key]] <- sample_map(results[[se_col]], sample_indices)
  out[[pvalue_key]] <- sample_map(results[[pvalue_col]], sample_indices)
  out
}

sim <- readRDS("baselines/example_sim_data.rds")
obj_nullmodel <- fit_null_glmmkin(
  Y ~ X1 + X2,
  data = sim$pheno_related,
  kins = sim$kins_sparse,
  id = "id",
  family = gaussian(link = "identity")
)

res <- Indiv_Score_Test_Region(
  genotype = sim$Geno,
  obj_nullmodel = obj_nullmodel,
  rare_maf_cutoff = 0.05
)

res_cond <- Indiv_Score_Test_Region_cond(
  genotype = sim$Geno,
  genotype_adj = as.matrix(sim$Geno[, adj_variant_indices]),
  obj_nullmodel = obj_nullmodel,
  rare_maf_cutoff = 0.05,
  method_cond = "optimal"
)

jsonlite::write_json(
  build_sentinels(
    res,
    score_col = "Score",
    se_col = "SE",
    pvalue_col = "pvalue",
    score_key = "score_samples",
    se_key = "se_samples",
    pvalue_key = "pvalue_samples"
  ),
  "baselines/indiv_score_related_sparse_glmmkin_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)

jsonlite::write_json(
  build_sentinels(
    res_cond,
    score_col = "Score_cond",
    se_col = "SE_cond",
    pvalue_col = "pvalue_cond",
    score_key = "score_cond_samples",
    se_key = "se_cond_samples",
    pvalue_key = "pvalue_cond_samples"
  ),
  "baselines/indiv_score_related_sparse_glmmkin_cond_sentinels.json",
  pretty = TRUE,
  digits = 15,
  auto_unbox = TRUE
)
