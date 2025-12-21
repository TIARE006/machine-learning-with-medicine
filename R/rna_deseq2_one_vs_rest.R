#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage:\n  Rscript R/rna_deseq2_one_vs_rest.R <RUN_DIR> <RNA_COUNTS_CSV>\n")
  quit(status = 1)
}

run_dir <- args[1]
counts_csv <- args[2]

labels_csv <- file.path(run_dir, "labels", "cluster_labels.csv")
out_dir <- file.path(run_dir, "degs_deseq2")
plots_dir <- file.path(run_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat(">>> [DESeq2] run_dir:", run_dir, "\n")
cat(">>> [DESeq2] labels:", labels_csv, "\n")
cat(">>> [DESeq2] counts:", counts_csv, "\n")

# -----------------------------
# Load labels
# -----------------------------
lab <- read_csv(labels_csv, show_col_types = FALSE) %>%
  mutate(Sample_ID = as.character(Sample_ID))

# -----------------------------
# Parse your special counts CSV
# Line1: ",Unnamed: 1,<bam1>,<bam2>,..."
# Line2: "FeatureID,type,AB_3,AB_4,..."
# Data starts from line3
# -----------------------------
first_line <- readLines(counts_csv, n = 1)
parts <- strsplit(first_line, ",", fixed = TRUE)[[1]]

# parts[1] is empty, parts[2] is "Unnamed: 1", rest are bam names
bam_names <- parts[-c(1, 2)]
bam_names <- trimws(bam_names)
bam_names <- bam_names[bam_names != ""]

cat(">>> [DESeq2] bam names from line1:", length(bam_names), "\n")

# Read actual table from line2 onward
df <- read_csv(counts_csv, skip = 1, show_col_types = FALSE)

# Expect FeatureID + type + AB_* columns
if (!all(c("FeatureID", "type") %in% colnames(df))) {
  stop("Counts file format unexpected: missing FeatureID/type columns after skipping first line.")
}

gene_ids <- df[["FeatureID"]]

counts_df <- df %>%
  select(-FeatureID, -type) %>%
  as.data.frame()

# Drop obvious junk columns (in case)
drop_cols <- grepl("^Unnamed", colnames(counts_df)) | grepl("^\\.\\.\\.[0-9]+$", colnames(counts_df))
if (any(drop_cols)) {
  cat(">>> [DESeq2] drop cols:", paste(colnames(counts_df)[drop_cols], collapse = ", "), "\n")
  counts_df <- counts_df[, !drop_cols, drop = FALSE]
}

# Rename AB_* columns to bam names (must match)
if (ncol(counts_df) != length(bam_names)) {
  cat(">>> [DESeq2] ERROR: ncol(counts_df)=", ncol(counts_df),
      " but bam_names=", length(bam_names), "\n")
  stop("Counts columns do not match bam mapping line. Please re-check the CSV.")
}
colnames(counts_df) <- bam_names

# Convert to numeric (invalid -> NA)
counts_df[] <- lapply(counts_df, function(x) suppressWarnings(as.numeric(x)))

# Handle NA
na_total <- sum(is.na(as.matrix(counts_df)))
if (na_total > 0) {
  cat(">>> [DESeq2] NA found in counts:", na_total, " -> set to 0\n")
  counts_df[is.na(counts_df)] <- 0
}

counts <- as.matrix(counts_df)
rownames(counts) <- gene_ids

# Clean column names (trim)
colnames(counts) <- trimws(colnames(counts))

# Remove fake samples like 'Unnamed:*' (extra safety)
bad <- grepl("^Unnamed", colnames(counts))
if (any(bad)) counts <- counts[, !bad, drop = FALSE]

# Check non-negative
minv <- min(counts)
if (minv < 0) {
  stop(paste0("Negative values detected in counts (min=", minv, "). This is not raw counts."))
}

# -----------------------------
# Align samples to labels
# -----------------------------
common <- intersect(lab$Sample_ID, colnames(counts))
cat(">>> [DESeq2] common samples with labels:", length(common), "\n")
if (length(common) < 4) stop("Too few matched samples between labels and counts.")

lab2 <- lab %>% filter(Sample_ID %in% common)

# reorder counts to label order
counts2 <- counts[, lab2$Sample_ID, drop = FALSE]
stopifnot(all(colnames(counts2) == lab2$Sample_ID))

# DESeq2 expects integer counts
counts2 <- round(counts2)
storage.mode(counts2) <- "integer"
stopifnot(all(counts2 >= 0))

# Filter low-expression genes
keep <- rowSums(counts2) >= 10
counts2 <- counts2[keep, , drop = FALSE]
cat(">>> [DESeq2] genes after filter (rowSums>=10):", nrow(counts2), "\n")

clusters <- sort(unique(lab2$Cluster))
cat(">>> [DESeq2] clusters:", paste(clusters, collapse = ","), "\n")
cat(">>> [DESeq2] cluster sizes:\n")
print(table(lab2$Cluster))

# -----------------------------
# Save VST matrix for heatmap
# -----------------------------
dds_all <- DESeqDataSetFromMatrix(
  countData = counts2,
  colData = data.frame(row.names = lab2$Sample_ID,
                       cluster = factor(lab2$Cluster)),
  design = ~ cluster
)
dds_all <- estimateSizeFactors(dds_all)
vsd <- vst(dds_all, blind = TRUE)
vst_mat <- assay(vsd)

write_csv(
  as.data.frame(vst_mat) %>% rownames_to_column("gene"),
  file.path(out_dir, "vst_matrix.csv")
)
cat(">>> [DESeq2] saved VST matrix:", file.path(out_dir, "vst_matrix.csv"), "\n")

# -----------------------------
# one-vs-rest DESeq2 per cluster
# -----------------------------
for (c in clusters) {
  group <- ifelse(lab2$Cluster == c, paste0("C", c), "REST")
  group <- factor(group, levels = c("REST", paste0("C", c)))

  dds <- DESeqDataSetFromMatrix(
    countData = counts2,
    colData = data.frame(row.names = lab2$Sample_ID, group = group),
    design = ~ group
  )

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", paste0("C", c), "REST"))
  res <- lfcShrink(dds,
                   contrast = c("group", paste0("C", c), "REST"),
                   res = res,
                   type = "normal")

  out <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    transmute(
      gene = gene,
      log2FC = log2FoldChange,
      p_value = pvalue,
      FDR = padj,
      baseMean = baseMean,
      lfcSE = lfcSE,
      stat = stat
    )

  out_file <- file.path(out_dir, sprintf("DEG_cluster_%s_vs_rest.csv", c))
  write_csv(out, out_file)
  cat(">>> [DESeq2] saved:", out_file, " rows=", nrow(out), "\n")
}

cat(">>> [DESeq2] done.\n")

