if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}

env_info <- list(
  sys_info = as.list(Sys.info()),
  cores_logical = parallel::detectCores(logical = TRUE),
  cores_physical = parallel::detectCores(logical = FALSE),
  r_version = R.version.string,
  platform = R.version$platform,
  session_info = capture.output(sessionInfo()),
  blas = extSoftVersion()["BLAS"],
  lapack = extSoftVersion()["LAPACK"],
  locale = Sys.getlocale(),
  timestamp = as.character(Sys.time())
)

jsonlite::write_json(env_info, "baselines/r_environment.json", pretty = TRUE)
