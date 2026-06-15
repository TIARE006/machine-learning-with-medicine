options(stringsAsFactors = FALSE)

.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(edgeR)
  library(limma)
  library(matrixStats)
  library(digest)
})

rna_path <- "data/raw/RNA_seq/GSE254877_raw_counts_expression.csv"
lncrna_path <- "data/raw/lncRNA_seq/GSE254877_lncRNA_raw_counts_expression.csv"
smallrna_path <- "data/raw/small_RNA_seq/GSE254878_smallRNAs_raw_counts_expression.csv"
pheno_path <- "data/raw/RNA_seq/GSE254877_pheno.csv"
feature_spec_path <- "scripts/fig1/01_expression_inputs/locked_intnmf_feature_spec.csv"
outdir <- "results/fig1_inputs/generated"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(rna_path), file.exists(lncrna_path), file.exists(smallrna_path), file.exists(pheno_path), file.exists(feature_spec_path))

pheno <- read_csv(pheno_path, show_col_types = FALSE, name_repair = "minimal")
sample_id <- str_extract(as.character(pheno$title), "AB_[0-9]+")
subgroup <- str_trim(str_replace(as.character(pheno$title), "^AB_[0-9]+,", ""))

meta <- tibble(
  sample_id = sample_id,
  geo_accession = as.character(pheno$geo_accession),
  title = as.character(pheno$title),
  subgroup = subgroup
) %>%
  filter(!is.na(sample_id)) %>%
  arrange(as.integer(str_remove(sample_id, "AB_")))

if (!identical(meta$sample_id, paste0("AB_", 1:84))) {
  stop("GSE254877 phenotype does not resolve uniquely to AB_1 through AB_84.")
}

meta_path <- file.path(outdir, "intNMF_sample_metadata.csv")
write_csv(meta, meta_path)

read_count_block <- function(path, label, allowed_types = NULL, uppercase_features = FALSE) {
  raw <- read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    name_repair = "unique"
  )

  first_row <- as.character(unlist(raw[1, , drop = FALSE], use.names = FALSE))
  sample_cols <- which(str_detect(first_row, "^AB_[0-9]+$"))

  if (length(sample_cols) != 84 || !all(diff(sample_cols) == 1)) {
    stop(label, ": expected 84 contiguous AB sample columns.")
  }

  dat <- raw[-1, , drop = FALSE]
  feature_col <- names(dat)[1]
  type_col <- if (min(sample_cols) >= 3) names(dat)[2] else NULL
  features <- as.character(dat[[feature_col]])
  if (uppercase_features) features <- toupper(features)

  keep <- !is.na(features) & features != ""
  if (!is.null(allowed_types)) {
    keep <- keep & as.character(dat[[type_col]]) %in% allowed_types
  }

  mat <- as.matrix(dat[, sample_cols, drop = FALSE])
  mat <- apply(mat, 2, function(x) suppressWarnings(as.numeric(x)))
  mat[is.na(mat)] <- 0
  rownames(mat) <- features
  colnames(mat) <- first_row[sample_cols]
  mat <- mat[keep, meta$sample_id, drop = FALSE]
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  storage.mode(mat) <- "numeric"
  mat
}

make_qc <- function(counts, prefix) {
  tibble(
    sample_id = colnames(counts),
    "{prefix}_library_size" := colSums(counts),
    "log10_{prefix}_library_size" := log10(colSums(counts) + 1),
    "{prefix}_detected_features" := colSums(counts > 0),
    "{prefix}_zero_fraction" := colMeans(counts == 0)
  )
}

prepare_block <- function(counts, prefix, min_count, min_samples, locked_features, locked_names) {
  counts <- counts[rowSums(counts >= min_count) >= min_samples, , drop = FALSE]
  dge <- DGEList(counts = round(counts))
  dge <- normLibSizes(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)

  covariates <- qc %>%
    arrange(match(sample_id, meta$sample_id)) %>%
    transmute(
      detected = as.numeric(scale(.data[[paste0(prefix, "_detected_features")]])),
      loglib = as.numeric(scale(.data[[paste0("log10_", prefix, "_library_size")]]))
    ) %>%
    as.matrix()

  design <- model.matrix(~ subgroup, data = meta)
  adjusted <- removeBatchEffect(logcpm, covariates = covariates, design = design)
  missing_locked <- setdiff(locked_features, rownames(adjusted))
  if (length(missing_locked) > 0) {
    stop(prefix, " locked features missing after filtering: ", paste(missing_locked, collapse = ", "))
  }
  selected <- adjusted[locked_features, , drop = FALSE]

  block <- t(selected)
  if (min(block) < 0) block <- block + abs(min(block))
  block <- pmax(block, .Machine$double.eps)
  block <- block / max(block)
  rownames(block) <- meta$sample_id
  colnames(block) <- locked_names
  block
}

mrna <- read_count_block(rna_path, "mRNA", allowed_types = "protein_coding", uppercase_features = TRUE)
lncrna <- read_count_block(lncrna_path, "lncRNA", uppercase_features = TRUE)
smallrna <- read_count_block(smallrna_path, "smallRNA")

qc <- make_qc(mrna, "RNA") %>%
  left_join(make_qc(lncrna, "lncRNA"), by = "sample_id") %>%
  left_join(make_qc(smallrna, "smallRNA"), by = "sample_id") %>%
  left_join(meta, by = "sample_id")

feature_spec <- read_csv(feature_spec_path, show_col_types = FALSE) %>%
  arrange(block, position)

mrna_spec <- feature_spec %>% filter(modality == "mRNA")
lncrna_spec <- feature_spec %>% filter(modality == "lncRNA")
smallrna_spec <- feature_spec %>% filter(modality == "smallRNA")

mrna_block <- prepare_block(
  mrna, "RNA", min_count = 10, min_samples = 5,
  locked_features = mrna_spec$feature,
  locked_names = mrna_spec$model_name
)

lncrna_block <- prepare_block(
  lncrna, "lncRNA", min_count = 5, min_samples = 5,
  locked_features = lncrna_spec$feature,
  locked_names = lncrna_spec$model_name
)

dat <- list(
  RNA_lncRNA = cbind(mrna_block, lncrna_block),
  smallRNA = prepare_block(
    smallrna, "smallRNA", min_count = 5, min_samples = 5,
    locked_features = smallrna_spec$feature,
    locked_names = smallrna_spec$model_name
  )
)

input_path <- file.path(outdir, "two_block_intNMF_inputs.rds")
qc_path <- file.path(outdir, "two_block_sample_QC.csv")
saveRDS(dat, input_path)
write_csv(qc, qc_path)

files <- c(rna_path, lncrna_path, smallrna_path, pheno_path, feature_spec_path, meta_path, input_path, qc_path)
write_csv(
  tibble(
    role = c("raw_RNA", "raw_lncRNA", "raw_smallRNA", "GEO_phenotype", "locked_feature_spec", "generated_metadata", "generated_intNMF_input", "generated_QC"),
    path = files,
    sha256 = vapply(files, digest, algo = "sha256", file = TRUE, FUN.VALUE = character(1))
  ),
  file.path(outdir, "intNMF_input_lineage_manifest.csv")
)

cat("Generated intNMF source inputs:", nrow(dat$RNA_lncRNA), "x", ncol(dat$RNA_lncRNA), "mRNA+lncRNA;",
    nrow(dat$smallRNA), "x", ncol(dat$smallRNA), "smallRNA\n")
