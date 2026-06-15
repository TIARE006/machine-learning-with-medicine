options(stringsAsFactors = FALSE)

.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(edgeR)
  library(GSVA)
  library(digest)
})

raw_path <- "data/raw/RNA_seq/GSE254877_raw_counts_expression.csv"
label_path <- "results/intNMF_reproducible_two_block/two_block_final_K4_labels.csv"
gene_set_path <- "scripts/fig1/01_expression_inputs/pcd_gene_sets.csv"
outdir <- "results/fig1_inputs/generated"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(raw_path), file.exists(label_path), file.exists(gene_set_path))

raw <- read_csv(
  raw_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

sample_ids <- as.character(unlist(raw[1, 3:ncol(raw)], use.names = FALSE))
dat <- raw[-1, , drop = FALSE]
genes <- toupper(as.character(dat[[1]]))
feature_type <- as.character(dat[[2]])

counts <- as.matrix(dat[, 3:ncol(raw), drop = FALSE])
counts <- apply(counts, 2, function(x) suppressWarnings(as.numeric(x)))
counts[is.na(counts)] <- 0
rownames(counts) <- genes
colnames(counts) <- sample_ids

keep <- !is.na(genes) & genes != "" & feature_type == "protein_coding"
counts <- rowsum(counts[keep, , drop = FALSE], group = genes[keep], reorder = FALSE)
storage.mode(counts) <- "numeric"
counts <- counts[rowSums(counts >= 10) >= 5, , drop = FALSE]

dge <- DGEList(counts = round(counts))
dge <- normLibSizes(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

labels <- read_csv(label_path, show_col_types = FALSE) %>%
  transmute(
    sample_id = as.character(sample_id),
    adjusted_two_block_K4 = as.character(adjusted_two_block_K4),
    geo_accession,
    title,
    subgroup
  )

if (nrow(labels) != 84 || anyDuplicated(labels$sample_id) > 0) {
  stop("Expected 84 unique final K4 labels.")
}

missing_samples <- setdiff(labels$sample_id, colnames(logcpm))
if (length(missing_samples) > 0) {
  stop("RNA matrix is missing samples: ", paste(missing_samples, collapse = ", "))
}

logcpm <- logcpm[, labels$sample_id, drop = FALSE]

gene_set_tbl <- read_csv(gene_set_path, show_col_types = FALSE) %>%
  mutate(program = as.character(program), gene = toupper(as.character(gene)))

gene_sets <- split(gene_set_tbl$gene, gene_set_tbl$program)
gene_sets <- lapply(gene_sets, unique)
present_sets <- lapply(gene_sets, intersect, y = rownames(logcpm))

gene_set_audit <- tibble(
  program = names(gene_sets),
  n_defined = lengths(gene_sets),
  n_present = lengths(present_sets),
  genes_defined = vapply(gene_sets, paste, collapse = ";", FUN.VALUE = character(1)),
  genes_present = vapply(present_sets, paste, collapse = ";", FUN.VALUE = character(1)),
  source_reference = vapply(
    names(gene_sets),
    function(x) paste(unique(gene_set_tbl$source_reference[gene_set_tbl$program == x]), collapse = ";"),
    FUN.VALUE = character(1)
  )
)

if (any(gene_set_audit$n_present < 3)) {
  stop("Every PCD program must have at least three expressed genes.")
}

write_csv(gene_set_audit, file.path(outdir, "GSE254877_14PCD_gene_set_audit.csv"))

ssgsea_param <- ssgseaParam(logcpm, present_sets, normalize = TRUE)
ssgsea <- gsva(ssgsea_param, verbose = FALSE)
pcd_z <- t(scale(t(ssgsea)))
pcd_z[!is.finite(pcd_z)] <- 0

pcd_out <- as.data.frame(t(pcd_z)) %>%
  rownames_to_column("sample_id") %>%
  left_join(labels, by = "sample_id")

write_csv(pcd_out, file.path(outdir, "GSE254877_14PCD_ssGSEA_zscores.csv"))

marker_genes <- c(
  "GSDMD", "GSDME", "GSDMB", "GSDMC", "CASP1", "CASP4", "CASP5",
  "NLRP3", "NLRP1", "NLRC4", "AIM2", "PYCARD", "IL1B", "IL18", "NOD1", "NOD2"
)

missing_markers <- setdiff(marker_genes, rownames(logcpm))
if (length(missing_markers) > 0) {
  stop("Missing pyroptosis marker genes: ", paste(missing_markers, collapse = ", "))
}

marker_z <- t(scale(t(logcpm[marker_genes, , drop = FALSE])))
marker_z[!is.finite(marker_z)] <- 0

write_csv(
  as.data.frame(marker_z) %>% rownames_to_column("gene"),
  file.path(outdir, "GSE254877_pyroptosis_marker_gene_zmatrix.csv")
)

files <- c(
  raw_path,
  label_path,
  gene_set_path,
  file.path(outdir, "GSE254877_14PCD_gene_set_audit.csv"),
  file.path(outdir, "GSE254877_14PCD_ssGSEA_zscores.csv"),
  file.path(outdir, "GSE254877_pyroptosis_marker_gene_zmatrix.csv")
)

lineage <- tibble(
  role = c("raw_expression", "final_K4_labels", "PCD_gene_set_definition", rep("generated_output", 3)),
  path = files,
  sha256 = vapply(files, digest, algo = "sha256", file = TRUE, FUN.VALUE = character(1))
)

write_csv(lineage, file.path(outdir, "expression_input_lineage_manifest.csv"))
cat("Generated traceable Fig. 1 expression inputs in:", outdir, "\n")
