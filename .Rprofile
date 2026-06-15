local_lib <- file.path(getwd(), "Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib, winslash = "/", mustWork = TRUE), .libPaths()))
}

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE,
  timeout = 1200
)
