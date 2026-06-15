proj <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_lib <- file.path(proj, "Rlibs")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(local_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Project:", proj, "\n")
cat("Local R library:", local_lib, "\n")
cat("Current .libPaths():\n")
print(.libPaths())

cran_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "tibble", "stringr",
  "purrr", "data.table", "pheatmap", "ggpubr", "patchwork",
  "RColorBrewer", "matrixStats", "scales", "viridisLite",
  "ggrepel", "cowplot", "digest", "cluster", "mclust", "NMF",
  "magick", "circlize"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("\nInstalling CRAN package:", p, "\n")
    install.packages(p, lib = local_lib, dependencies = TRUE)
  } else {
    cat("Already installed:", p, "\n")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\nInstalling BiocManager\n")
  install.packages("BiocManager", lib = local_lib, dependencies = TRUE)
}

.libPaths(c(local_lib, .libPaths()))

bioc_pkgs <- c(
  "GSVA", "GSEABase", "ComplexHeatmap",
  "ConsensusClusterPlus", "limma", "sva", "BiocParallel",
  "edgeR", "IntNMF"
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("\nInstalling Bioconductor package:", p, "\n")
    BiocManager::install(p, lib = local_lib, ask = FALSE, update = FALSE)
  } else {
    cat("Already installed:", p, "\n")
  }
}

cat("\n===== Final package check =====\n")
all_pkgs <- c(cran_pkgs, "BiocManager", bioc_pkgs)
print(sapply(all_pkgs, requireNamespace, quietly = TRUE))

cat("\nDone.\n")
