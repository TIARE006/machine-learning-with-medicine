suppressPackageStartupMessages({
  library(GSVA)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
})

# -----------------------------
# 1) paths
# -----------------------------
project_root <- "."
expr_file  <- file.path(project_root, "data/raw/RNA_seq/GSE254877_counts_matrix_for_cps.csv")
pheno_file <- file.path(project_root, "data/raw/RNA_seq/GSE254877_pheno_ab_for_cps.csv")
out_dir    <- file.path(project_root, "results/cps_signature/GSE254877_ssgsea_r")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) gene sets
# -----------------------------
positive_genes <- c("GSDMD", "CASP1", "NLRP3", "PYCARD", "IL1B",
                    "IL18", "CASP4", "CASP5", "CCL2")
control_genes  <- c("GSDME", "CASP3", "CD68")

gene_sets <- list(
  positive = positive_genes,
  control  = control_genes
)

# -----------------------------
# 3) read expression
#    expected: FeatureID + sample columns
# -----------------------------
expr_df <- read_csv(expr_file, show_col_types = FALSE)

if (!"FeatureID" %in% colnames(expr_df)) {
  stop("Expression file must contain a 'FeatureID' column.")
}

expr_df <- expr_df %>%
  mutate(FeatureID = toupper(trimws(as.character(FeatureID)))) %>%
  filter(FeatureID != "")

sample_cols <- setdiff(colnames(expr_df), "FeatureID")

expr_mat <- expr_df %>%
  select(FeatureID, all_of(sample_cols)) %>%
  distinct(FeatureID, .keep_all = TRUE) %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()

mode(expr_mat) <- "numeric"

# -----------------------------
# 4) read pheno
#    expected: sample_id, subgroup
# -----------------------------
pheno <- read_csv(pheno_file, show_col_types = FALSE)

required_cols <- c("sample_id", "subgroup")
missing_cols <- setdiff(required_cols, colnames(pheno))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in pheno:", paste(missing_cols, collapse = ", ")))
}

pheno <- pheno %>%
  mutate(
    sample_id = as.character(sample_id),
    subgroup  = as.character(subgroup)
  )

# keep overlapping samples only
common_samples <- intersect(colnames(expr_mat), pheno$sample_id)
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between expression matrix and pheno.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
pheno <- pheno %>%
  filter(sample_id %in% common_samples) %>%
  mutate(sample_id = factor(sample_id, levels = common_samples)) %>%
  arrange(sample_id) %>%
  mutate(sample_id = as.character(sample_id))

# -----------------------------
# 5) keep present genes
# -----------------------------
positive_present <- intersect(positive_genes, rownames(expr_mat))
control_present  <- intersect(control_genes, rownames(expr_mat))

if (length(positive_present) < 2) {
  stop("Too few positive genes present in expression matrix.")
}
if (length(control_present) < 2) {
  stop("Too few control genes present in expression matrix.")
}

gene_sets_present <- list(
  positive = positive_present,
  control  = control_present
)

# -----------------------------
# 6) standard ssGSEA with GSVA
#    counts matrix -> kcdf = "Poisson"
# -----------------------------
es <- gsva(
  expr = expr_mat,
  gset.idx.list = gene_sets_present,
  method = "ssgsea",
  kcdf = "Poisson",
  verbose = FALSE
)

# es: gene_set x sample
es_df <- as.data.frame(t(es)) %>%
  rownames_to_column("sample_id")

# -----------------------------
# 7) CPS
# -----------------------------
scores <- es_df %>%
  transmute(
    sample_id   = sample_id,
    ES_positive = positive,
    ES_control  = control,
    CPS         = positive - control
  ) %>%
  mutate(
    CPS_zscore = as.numeric(scale(CPS))
  ) %>%
  left_join(pheno, by = "sample_id") %>%
  mutate(
    CPS_class = ifelse(CPS_zscore > 0.5, "CPS-High", "CPS-LowOrMid")
  )

# -----------------------------
# 8) summaries
# -----------------------------
group_summary <- scores %>%
  group_by(subgroup) %>%
  summarise(
    n = n(),
    mean_ES_positive = mean(ES_positive, na.rm = TRUE),
    mean_ES_control  = mean(ES_control, na.rm = TRUE),
    mean_CPS         = mean(CPS, na.rm = TRUE),
    median_CPS       = median(CPS, na.rm = TRUE),
    sd_CPS           = sd(CPS, na.rm = TRUE),
    mean_CPS_zscore  = mean(CPS_zscore, na.rm = TRUE),
    n_CPS_high       = sum(CPS_class == "CPS-High", na.rm = TRUE),
    frac_CPS_high    = mean(CPS_class == "CPS-High", na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 9) outputs
# -----------------------------
write_csv(scores, file.path(out_dir, "sample_scores_ssgsea_r.csv"))
write_csv(group_summary, file.path(out_dir, "group_summary_ssgsea_r.csv"))

run_log <- c(
  paste("expression_file:", expr_file),
  paste("pheno_file:", pheno_file),
  "method: GSVA::gsva with ssgseaParam()",
  "kcdf: Poisson",
  paste("n_genes_total:", nrow(expr_mat)),
  paste("n_samples_total:", ncol(expr_mat)),
  "",
  paste("positive_genes_requested:", paste(positive_genes, collapse = ", ")),
  paste("positive_genes_present:", paste(positive_present, collapse = ", ")),
  "",
  paste("control_genes_requested:", paste(control_genes, collapse = ", ")),
  paste("control_genes_present:", paste(control_present, collapse = ", ")),
  "",
  "CPS = ES_positive - ES_control",
  "CPS_high threshold: CPS_zscore > 0.5"
)
writeLines(run_log, con = file.path(out_dir, "run_log_ssgsea_r.txt"))

# -----------------------------
# 10) boxplot
# -----------------------------
png(file.path(out_dir, "cps_boxplot_ssgsea_r.png"), width = 1600, height = 1000, res = 150)
boxplot(
  CPS ~ subgroup,
  data = scores,
  main = "CPS by subgroup (standard ssGSEA, GSVA)",
  xlab = "subgroup",
  ylab = "CPS",
  las = 2
)
dev.off()

# -----------------------------
# 11) heatmap
# -----------------------------
heat_df <- scores %>%
  arrange(subgroup, desc(CPS)) %>%
  select(sample_id, subgroup, ES_positive, ES_control, CPS, CPS_zscore)

heat_mat <- heat_df %>%
  select(ES_positive, ES_control, CPS, CPS_zscore) %>%
  as.matrix()

rownames(heat_mat) <- paste(heat_df$sample_id, heat_df$subgroup, sep = " | ")

ann_row <- data.frame(subgroup = heat_df$subgroup)
rownames(ann_row) <- rownames(heat_mat)

png(file.path(out_dir, "cps_heatmap_ssgsea_r.png"), width = 1800, height = 2400, res = 180)
pheatmap(
  heat_mat,
  annotation_row = ann_row,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 7,
  main = "Sample-level CPS heatmap (standard ssGSEA, GSVA)"
)
dev.off()

message("Done.")
message("Output dir: ", out_dir)
message("sample_scores_ssgsea_r.csv")
message("group_summary_ssgsea_r.csv")
message("cps_boxplot_ssgsea_r.png")
message("cps_heatmap_ssgsea_r.png")