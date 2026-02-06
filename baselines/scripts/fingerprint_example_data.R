if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest", repos = "https://cloud.r-project.org")
}

suppressMessages(data("example", package = "STAAR"))

fingerprint <- list(
  name = "example",
  source = "STAAR package dataset",
  components = list(
    genotype = list(
      class = class(example$genotype),
      dim = dim(example$genotype)
    ),
    maf = list(
      class = class(example$maf),
      length = length(example$maf)
    ),
    snploc = list(
      class = class(example$snploc),
      dim = dim(example$snploc),
      colnames = colnames(example$snploc)
    )
  ),
  checksum = digest::digest(example, algo = "sha256")
)

jsonlite::write_json(
  fingerprint,
  "baselines/example_fingerprint.json",
  pretty = TRUE,
  digits = 15
)
