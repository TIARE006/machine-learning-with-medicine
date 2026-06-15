options(stringsAsFactors = FALSE, timeout = 1800)
set.seed(20260607)

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
LOCAL_LIB <- file.path(PROJECT_ROOT, "Rlibs")
dir.create(LOCAL_LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(LOCAL_LIB, winslash = "/", mustWork = TRUE), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "circlize", "cluster", "magick"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = LOCAL_LIB, dependencies = c("Depends", "Imports", "LinkingTo"))
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = LOCAL_LIB)
}

bioc_pkgs <- c("edgeR", "limma", "ConsensusClusterPlus", "ComplexHeatmap")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib = LOCAL_LIB, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(circlize)
  library(cluster)
  library(edgeR)
  library(limma)
  library(ConsensusClusterPlus)
  library(ComplexHeatmap)
  library(grid)
  library(magick)
})

OUTDIR <- file.path(
  PROJECT_ROOT,
  "results",
  "fig1_final",
  "all84_pyro_module_Euclidean_K2_main"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

RNA_PATH <- file.path(
  PROJECT_ROOT,
  "data",
  "raw",
  "RNA_seq",
  "GSE254877_raw_counts_expression.csv"
)

LABEL_PATH <- file.path(
  PROJECT_ROOT,
  "results",
  "intNMF_reproducible_two_block",
  "two_block_final_K4_labels.csv"
)

META_INFERRED_PATH <- file.path(
  PROJECT_ROOT,
  "results",
  "fig1_final",
  "reconstructed_Subtype1_pyro_K2",
  "GSE254877_sample_metadata_with_inferred_sex.csv"
)

stopifnot(file.exists(RNA_PATH))
stopifnot(file.exists(LABEL_PATH))

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

safe_factor <- function(x) {
  factor(as.character(x))
}

# ============================================================
# 0. Locked analysis design record
# ============================================================

design_lines <- c(
  "LOCKED ALL-84 PYROPTOSIS MODULE ANALYSIS",
  "",
  "Purpose:",
  "Identify pyroptosis-module-high versus pyroptosis-module-low molecular states within all 84 cancer-cachexia cohort samples.",
  "",
  "Anti-cherry-picking rules:",
  "1. Use all 84 samples. No Subtype1-like prefiltering.",
  "2. K=2 is prespecified as a mechanism-driven high/low stratification, not claimed as the global optimal K.",
  "3. K=2/3/4 consensus metrics are reported together.",
  "4. Clustering uses only prespecified module genes.",
  "5. GSDMD, GSDME, IL1B and IL18 are held out from clustering and used only for validation.",
  "6. Validation is reported with sex, cancer type and myeloid-adjusted models.",
  "7. If held-out validation fails, the groups should not be called PA/PQ."
)

writeLines(
  design_lines,
  file.path(OUTDIR, "all84_locked_analysis_design.txt")
)

# ============================================================
# 1. Read RNA and reconstruct logCPM
# ============================================================

raw <- read_csv(
  RNA_PATH,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

feature_col <- names(raw)[1]
type_col <- names(raw)[2]

sample_ids <- as.character(unlist(raw[1, 3:ncol(raw)], use.names = FALSE))
dat <- raw[-1, , drop = FALSE]

genes <- toupper(as.character(dat[[feature_col]]))
feature_type <- as.character(dat[[type_col]])

counts <- as.matrix(dat[, 3:ncol(raw), drop = FALSE])
counts <- apply(counts, 2, function(v) suppressWarnings(as.numeric(v)))
counts[is.na(counts)] <- 0

rownames(counts) <- genes
colnames(counts) <- sample_ids

keep <- !is.na(genes) & genes != "" & feature_type == "protein_coding"
counts <- counts[keep, , drop = FALSE]
counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

counts <- counts[rowSums(counts >= 10) >= 5, , drop = FALSE]

dge <- DGEList(counts = round(counts))
dge <- normLibSizes(dge)

logcpm <- cpm(dge, log = TRUE, prior.count = 1)

write_csv(
  tibble(
    item = c("genes_after_filter", "samples"),
    value = c(nrow(logcpm), ncol(logcpm))
  ),
  file.path(OUTDIR, "all84_logCPM_matrix_audit.csv")
)

# ============================================================
# 2. Read all-84 sample metadata
# ============================================================

labels84 <- read_csv(LABEL_PATH, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    subgroup = as.character(subgroup),
    adjusted_two_block_K4 = as.character(adjusted_two_block_K4)
  )

if (nrow(labels84) != 84 || anyDuplicated(labels84$sample_id) > 0) {
  stop("Expected 84 unique samples in label file.")
}

missing_samples <- setdiff(labels84$sample_id, colnames(logcpm))

if (length(missing_samples) > 0) {
  stop("Samples missing from expression matrix: ", paste(missing_samples, collapse = ", "))
}

logcpm <- logcpm[, labels84$sample_id, drop = FALSE]

# Use existing expression-inferred sex if available; otherwise infer from Y markers.
if (file.exists(META_INFERRED_PATH)) {
  sex_meta <- read_csv(META_INFERRED_PATH, show_col_types = FALSE) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    dplyr::select(sample_id, sex, sex_source) %>%
    distinct(sample_id, .keep_all = TRUE)

  meta84 <- labels84 %>%
    left_join(sex_meta, by = "sample_id")

} else {
  y_markers <- intersect(
    c("RPS4Y1", "KDM5D", "DDX3Y", "UTY", "EIF1AY", "ZFY", "USP9Y", "TMSB4Y", "NLGN4Y", "PRKY"),
    rownames(logcpm)
  )

  x_markers <- intersect("XIST", rownames(logcpm))

  if (length(y_markers) < 3) {
    stop("Cannot infer sex: fewer than three Y markers found.")
  }

  sex_gene_z <- zscore_rows(logcpm[unique(c(y_markers, x_markers)), , drop = FALSE])
  y_score <- colMeans(sex_gene_z[y_markers, , drop = FALSE])

  if (length(x_markers) > 0) {
    xist_score <- colMeans(sex_gene_z[x_markers, , drop = FALSE])
    sex_axis <- y_score - xist_score
  } else {
    sex_axis <- y_score
  }

  set.seed(20260607)
  sex_km <- kmeans(matrix(sex_axis, ncol = 1), centers = 2, nstart = 100)
  cluster_means <- tapply(sex_axis, sex_km$cluster, mean)
  male_cluster <- as.integer(names(which.max(cluster_means)))

  inferred_sex <- ifelse(sex_km$cluster == male_cluster, "Male", "Female")

  meta84 <- labels84 %>%
    mutate(
      sex = inferred_sex[match(sample_id, names(sex_axis))],
      sex_source = "Expression_inferred_in_all84_script"
    )
}

if (any(is.na(meta84$sex))) {
  stop("Missing sex values after metadata merge/inference.")
}

meta84 <- meta84 %>%
  mutate(
    sex = factor(as.character(sex), levels = c("Female", "Male")),
    subgroup = factor(as.character(subgroup)),
    adjusted_two_block_K4 = factor(as.character(adjusted_two_block_K4))
  )

write_csv(
  meta84,
  file.path(OUTDIR, "all84_sample_metadata_used.csv")
)

cat("\n===== all84 metadata audit =====\n")
print(table(meta84$sex))
print(table(meta84$subgroup))
print(table(meta84$adjusted_two_block_K4))

# ============================================================
# 3. Prespecified pyroptosis module scores
# ============================================================

pyro_modules <- list(
  Inflammasome_sensor_adapter = c(
    "NLRP1", "NLRP3", "NLRC4", "AIM2", "PYCARD"
  ),

  Inflammatory_caspase_axis = c(
    "CASP1", "CASP4", "CASP5"
  ),

  NOD_gasdermin_related_axis = c(
    "NOD1", "NOD2", "GSDMB", "GSDMC"
  )
)

heldout_genes <- c(
  "GSDMD", "GSDME", "IL1B", "IL18"
)

pyro_modules <- lapply(
  pyro_modules,
  function(g) intersect(g, rownames(logcpm))
)

heldout_present <- intersect(heldout_genes, rownames(logcpm))

module_audit <- tibble(
  module = names(pyro_modules),
  n_genes = lengths(pyro_modules),
  genes = vapply(pyro_modules, paste, collapse = ";", FUN.VALUE = character(1))
)

heldout_audit <- tibble(
  role = "Held_out_validation",
  gene = heldout_present
)

write_csv(module_audit, file.path(OUTDIR, "all84_pyro_module_gene_audit.csv"))
write_csv(heldout_audit, file.path(OUTDIR, "all84_heldout_validation_gene_audit.csv"))

if (any(module_audit$n_genes < 2)) {
  stop("At least one pyroptosis module has fewer than two expressed genes.")
}

if (length(heldout_present) < 3) {
  stop("Too few held-out validation genes are available.")
}

all_discovery_genes <- unique(unlist(pyro_modules))

gene_z <- zscore_rows(
  logcpm[all_discovery_genes, meta84$sample_id, drop = FALSE]
)

module_score <- vapply(
  pyro_modules,
  function(g) colMeans(gene_z[g, , drop = FALSE]),
  FUN.VALUE = numeric(ncol(gene_z))
)

module_score <- t(module_score)
colnames(module_score) <- meta84$sample_id

module_score_z <- zscore_rows(module_score)
integrated_score <- colMeans(module_score_z)

module_score_df <- as.data.frame(t(module_score_z)) %>%
  rownames_to_column("sample_id") %>%
  mutate(integrated_pyro_module_score = integrated_score[sample_id]) %>%
  left_join(meta84, by = "sample_id")

write_csv(
  module_score_df,
  file.path(OUTDIR, "all84_pyro_module_scores_by_sample.csv")
)

# ============================================================
# 4. Consensus clustering K=2,3,4
# ============================================================

consensus_dir <- file.path(OUTDIR, "all84_module_score_euclidean_consensus")
dir.create(consensus_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260607)

cc <- ConsensusClusterPlus(
  as.matrix(module_score_z),
  maxK = 4,
  reps = 1000,
  pItem = 0.80,
  pFeature = 1.00,
  clusterAlg = "hc",
  distance = "euclidean",
  innerLinkage = "average",
  finalLinkage = "average",
  seed = 20260607,
  plot = "pdf",
  title = consensus_dir,
  verbose = FALSE
)

calculate_consensus_metrics <- function(k) {
  cls <- cc[[k]]$consensusClass

  if (is.null(names(cls))) {
    names(cls) <- colnames(module_score_z)
  }

  cm <- cc[[k]]$consensusMatrix
  rownames(cm) <- colnames(module_score_z)
  colnames(cm) <- colnames(module_score_z)

  cm <- cm[names(cls), names(cls), drop = FALSE]

  group <- factor(cls[names(cls)])
  upper <- upper.tri(cm)
  same <- outer(group, group, "==")

  within <- cm[upper & same]
  between <- cm[upper & !same]

  pac <- mean(cm[upper] > 0.10 & cm[upper] < 0.90)

  sil <- silhouette(
    as.integer(group),
    as.dist(1 - cm)
  )

  tibble(
    K = k,
    n_clusters = length(unique(group)),
    cluster_sizes = paste(
      paste0("C", names(table(group)), "=", as.integer(table(group))),
      collapse = ";"
    ),
    min_cluster_size = min(table(group)),
    max_cluster_size = max(table(group)),
    mean_within_cluster_consensus = mean(within),
    mean_between_cluster_consensus = mean(between),
    PAC_0.1_0.9 = pac,
    mean_consensus_silhouette = mean(sil[, "sil_width"])
  )
}

consensus_metrics <- bind_rows(
  lapply(2:4, calculate_consensus_metrics)
)

write_csv(
  consensus_metrics,
  file.path(OUTDIR, "all84_K2_K4_consensus_metrics.csv")
)

# ============================================================
# 5. Prespecified K=2 labels
# ============================================================

k2 <- cc[[2]]$consensusClass

if (is.null(names(k2))) {
  names(k2) <- colnames(module_score_z)
}

cluster_mean <- tapply(
  integrated_score[names(k2)],
  k2,
  mean
)

high_cluster <- names(which.max(cluster_mean))

labels <- tibble(
  sample_id = names(k2),
  module_K2_cluster = paste0("K2_", as.integer(k2)),
  integrated_pyro_module_score = integrated_score[names(k2)],
  Pyro_module_group = ifelse(
    as.character(k2) == high_cluster,
    "PyroModuleHigh",
    "PyroModuleLow"
  )
) %>%
  mutate(
    Pyro_module_group = factor(
      Pyro_module_group,
      levels = c("PyroModuleLow", "PyroModuleHigh")
    )
  ) %>%
  left_join(meta84, by = "sample_id")

write_csv(
  labels,
  file.path(OUTDIR, "all84_pyro_module_Euclidean_K2_labels.csv")
)

cat("\n===== all84 Pyro module Euclidean K2 groups =====\n")
print(table(labels$Pyro_module_group))

cat("\nBy sex:\n")
print(table(labels$Pyro_module_group, labels$sex))

cat("\nBy cancer subgroup:\n")
print(table(labels$Pyro_module_group, labels$subgroup))

cat("\nBy final K4:\n")
print(table(labels$Pyro_module_group, labels$adjusted_two_block_K4))

# K2 consensus matrix
consensus_matrix <- cc[[2]]$consensusMatrix
rownames(consensus_matrix) <- colnames(module_score_z)
colnames(consensus_matrix) <- colnames(module_score_z)

consensus_matrix <- consensus_matrix[
  labels$sample_id,
  labels$sample_id,
  drop = FALSE
]

# ============================================================
# 6. Validation scores and models
# ============================================================

validation_meta <- labels %>%
  arrange(match(sample_id, colnames(logcpm))) %>%
  mutate(
    Pyro_module_group = factor(
      Pyro_module_group,
      levels = c("PyroModuleLow", "PyroModuleHigh")
    ),
    sex = factor(sex),
    subgroup = factor(subgroup)
  )

validation_samples <- validation_meta$sample_id

heldout_expr <- logcpm[heldout_present, validation_samples, drop = FALSE]

myeloid_genes <- intersect(
  c("CD68", "LST1", "TYROBP", "FCER1G", "CTSS", "AIF1"),
  rownames(logcpm)
)

inflammation_genes <- intersect(
  c("IL6", "TNF", "NFKB1", "RELA", "STAT3", "CXCL8", "CCL2", "S100A8", "S100A9"),
  rownames(logcpm)
)

muscle_atrophy_genes <- intersect(
  c("FBXO32", "TRIM63", "FOXO1", "FOXO3", "MSTN", "CEBPB", "SOCS3"),
  rownames(logcpm)
)

extra_gene_sets <- list(
  Heldout_pyroptosis = heldout_present,
  Inflammation = inflammation_genes,
  Myeloid_infiltration = myeloid_genes,
  Muscle_atrophy = muscle_atrophy_genes
)

extra_all_genes <- unique(unlist(extra_gene_sets))

extra_z <- zscore_rows(
  logcpm[extra_all_genes, validation_samples, drop = FALSE]
)

extra_scores <- vapply(
  extra_gene_sets,
  function(g) colMeans(extra_z[g, , drop = FALSE]),
  FUN.VALUE = numeric(length(validation_samples))
)

extra_scores_df <- as.data.frame(extra_scores) %>%
  rownames_to_column("sample_id") %>%
  left_join(
    validation_meta %>%
      dplyr::select(sample_id, Pyro_module_group, sex, subgroup, adjusted_two_block_K4),
    by = "sample_id"
  )

write_csv(
  extra_scores_df,
  file.path(OUTDIR, "all84_K2_validation_scores_by_sample.csv")
)

# Held-out genes: base model
design_base <- model.matrix(
  ~ sex + subgroup + Pyro_module_group,
  data = validation_meta
)

colnames(design_base) <- make.names(colnames(design_base))

coef_base <- grep(
  "Pyro_module_group.*PyroModuleHigh",
  colnames(design_base),
  value = TRUE
)

fit_base <- eBayes(lmFit(heldout_expr, design_base))

heldout_base <- topTable(
  fit_base,
  coef = coef_base,
  number = Inf,
  sort.by = "none"
) %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    model = "Base_adjusted_for_sex_and_cancer_type",
    adjusted_logFC = logFC,
    average_logCPM = AveExpr,
    moderated_t = t,
    P_value = P.Value,
    FDR = adj.P.Val
  )

# Held-out genes: myeloid-adjusted model
validation_meta$myeloid_score <- extra_scores_df$Myeloid_infiltration[
  match(validation_meta$sample_id, extra_scores_df$sample_id)
]

design_myeloid <- model.matrix(
  ~ sex + subgroup + myeloid_score + Pyro_module_group,
  data = validation_meta
)

colnames(design_myeloid) <- make.names(colnames(design_myeloid))

coef_myeloid <- grep(
  "Pyro_module_group.*PyroModuleHigh",
  colnames(design_myeloid),
  value = TRUE
)

fit_myeloid <- eBayes(lmFit(heldout_expr, design_myeloid))

heldout_myeloid <- topTable(
  fit_myeloid,
  coef = coef_myeloid,
  number = Inf,
  sort.by = "none"
) %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    model = "Additional_adjustment_for_myeloid_score",
    adjusted_logFC = logFC,
    average_logCPM = AveExpr,
    moderated_t = t,
    P_value = P.Value,
    FDR = adj.P.Val
  )

heldout_stats <- bind_rows(heldout_base, heldout_myeloid) %>%
  arrange(gene, model)

write_csv(
  heldout_stats,
  file.path(OUTDIR, "all84_K2_heldout_gene_validation_limma.csv")
)

# Signature-level validation
test_signature <- function(sig_name, adjust_myeloid = FALSE) {
  dat <- extra_scores_df %>%
    transmute(
      sample_id,
      score = .data[[sig_name]],
      sex = factor(sex),
      subgroup = factor(subgroup),
      Pyro_module_group = factor(
        Pyro_module_group,
        levels = c("PyroModuleLow", "PyroModuleHigh")
      ),
      myeloid_score = Myeloid_infiltration
    )

  if (adjust_myeloid && sig_name != "Myeloid_infiltration") {
    fit <- lm(score ~ sex + subgroup + myeloid_score + Pyro_module_group, data = dat)
    model_name <- "myeloid_adjusted"
  } else {
    fit <- lm(score ~ sex + subgroup + Pyro_module_group, data = dat)
    model_name <- "base"
  }

  sm <- summary(fit)$coefficients

  coef_name <- grep(
    "Pyro_module_group.*PyroModuleHigh",
    rownames(sm),
    value = TRUE
  )

  tibble(
    signature = sig_name,
    model = model_name,
    n_genes = length(extra_gene_sets[[sig_name]]),
    PyroModuleLow_mean = mean(dat$score[dat$Pyro_module_group == "PyroModuleLow"]),
    PyroModuleHigh_mean = mean(dat$score[dat$Pyro_module_group == "PyroModuleHigh"]),
    adjusted_beta = sm[coef_name, "Estimate"],
    adjusted_SE = sm[coef_name, "Std. Error"],
    P_value = sm[coef_name, "Pr(>|t|)"]
  )
}

signature_stats <- bind_rows(
  lapply(names(extra_gene_sets), function(s) test_signature(s, FALSE)),
  lapply(setdiff(names(extra_gene_sets), "Myeloid_infiltration"), function(s) test_signature(s, TRUE))
) %>%
  group_by(model) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  arrange(signature, model)

write_csv(
  signature_stats,
  file.path(OUTDIR, "all84_K2_validation_signature_statistics.csv")
)

# ============================================================
# 7. Figure panels
# ============================================================

# Panel A: analysis design
panel_A <- ggplot() +
  annotate("text", x = 0, y = 5.3, label = "84 cancer-cachexia muscle samples", size = 5.5, fontface = "bold") +
  annotate("segment", x = 0, xend = 0, y = 5.0, yend = 4.25, arrow = arrow(length = unit(0.18, "in"))) +
  annotate("text", x = 0, y = 4.0, label = "Prespecified pyroptosis modules", size = 4.6, fontface = "bold") +
  annotate("text", x = 0, y = 3.62, label = "sensor/adaptor  |  inflammatory caspase  |  NOD/gasdermin-related", size = 3.4) +
  annotate("segment", x = 0, xend = 0, y = 3.35, yend = 2.6, arrow = arrow(length = unit(0.18, "in"))) +
  annotate("text", x = 0, y = 2.35, label = "Euclidean consensus K=2", size = 4.6, fontface = "bold") +
  annotate("text", x = 0, y = 2.0, label = "PyroModuleLow vs PyroModuleHigh", size = 3.8) +
  annotate("segment", x = 0, xend = 0, y = 1.75, yend = 1.0, arrow = arrow(length = unit(0.18, "in"))) +
  annotate("text", x = 0, y = 0.75, label = "Held-out validation", size = 4.6, fontface = "bold") +
  annotate("text", x = 0, y = 0.38, label = "GSDMD, GSDME, IL1B, IL18 + myeloid-adjusted model", size = 3.4) +
  theme_void() +
  xlim(-3.6, 3.6) +
  ylim(0, 5.7)

ggsave(
  file.path(OUTDIR, "Fig1A_all84_analysis_design.png"),
  panel_A,
  width = 6,
  height = 5,
  dpi = 300
)

# Panel B: K2-K4 metrics
metrics_long <- consensus_metrics %>%
  dplyr::select(
    K,
    mean_within_cluster_consensus,
    mean_between_cluster_consensus,
    mean_consensus_silhouette,
    PAC_0.1_0.9,
    min_cluster_size
  ) %>%
  pivot_longer(
    cols = -K,
    names_to = "metric",
    values_to = "value"
  )

panel_B <- ggplot(metrics_long, aes(x = factor(K), y = value, group = metric)) +
  geom_line(aes(color = metric), linewidth = 0.8) +
  geom_point(aes(color = metric), size = 2.5) +
  theme_bw(base_size = 11) +
  labs(
    title = "Consensus diagnostics across K=2-4",
    subtitle = "K=2 is prespecified as high/low pyroptosis-module stratification",
    x = "K",
    y = "Metric value",
    color = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTDIR, "Fig1B_all84_K2_K4_consensus_metrics.png"),
  panel_B,
  width = 7,
  height = 5,
  dpi = 300
)

# Panel C: module heatmap
plot_order <- labels %>%
  arrange(desc(Pyro_module_group), desc(integrated_pyro_module_score), subgroup, sex) %>%
  pull(sample_id)

ann <- labels %>%
  dplyr::select(sample_id, Pyro_module_group, sex, subgroup, adjusted_two_block_K4) %>%
  as.data.frame()

rownames(ann) <- ann$sample_id
ann$sample_id <- NULL

make_module_heatmap <- function() {
  Heatmap(
    module_score_z[, plot_order, drop = FALSE],
    name = "Module z",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    top_annotation = HeatmapAnnotation(df = ann[plot_order, , drop = FALSE]),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_split = ann[plot_order, "Pyro_module_group"],
    column_gap = unit(2, "mm"),
    column_title = "Pyroptosis module scores across all 84 samples",
    row_names_gp = gpar(fontsize = 9)
  )
}

pdf(file.path(OUTDIR, "Fig1C_all84_module_score_heatmap.pdf"), width = 10, height = 4.8, useDingbats = FALSE)
draw(make_module_heatmap(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(file.path(OUTDIR, "Fig1C_all84_module_score_heatmap.png"), width = 3000, height = 1450, res = 300)
draw(make_module_heatmap(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Panel D: consensus heatmap
make_consensus_heatmap <- function() {
  Heatmap(
    consensus_matrix[plot_order, plot_order, drop = FALSE],
    name = "Consensus",
    col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#9ECAE1", "#08306B")),
    top_annotation = HeatmapAnnotation(df = ann[plot_order, , drop = FALSE]),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = ann[plot_order, "Pyro_module_group"],
    column_split = ann[plot_order, "Pyro_module_group"],
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),
    column_title = "Consensus matrix: all-84 pyroptosis module K=2"
  )
}

pdf(file.path(OUTDIR, "Fig1D_all84_K2_consensus_heatmap.pdf"), width = 8, height = 7, useDingbats = FALSE)
draw(make_consensus_heatmap(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(file.path(OUTDIR, "Fig1D_all84_K2_consensus_heatmap.png"), width = 2400, height = 2100, res = 300)
draw(make_consensus_heatmap(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Panel E: validation signature boxplots
validation_long <- extra_scores_df %>%
  pivot_longer(
    cols = all_of(names(extra_gene_sets)),
    names_to = "signature",
    values_to = "score"
  )

panel_E <- ggplot(
  validation_long,
  aes(Pyro_module_group, score, fill = Pyro_module_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.13, size = 1, alpha = 0.65) +
  facet_wrap(~ signature, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "Held-out and secondary validation",
    x = NULL,
    y = "Mean gene z-score"
  )

ggsave(
  file.path(OUTDIR, "Fig1E_all84_validation_signature_boxplots.png"),
  panel_E,
  width = 8,
  height = 6,
  dpi = 300
)

# Panel F: held-out gene boxplots
heldout_long <- as.data.frame(t(heldout_expr)) %>%
  rownames_to_column("sample_id") %>%
  left_join(
    validation_meta %>%
      dplyr::select(sample_id, Pyro_module_group),
    by = "sample_id"
  ) %>%
  pivot_longer(
    cols = all_of(heldout_present),
    names_to = "gene",
    values_to = "logCPM"
  )

panel_F <- ggplot(
  heldout_long,
  aes(Pyro_module_group, logCPM, fill = Pyro_module_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.13, size = 1, alpha = 0.65) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "Held-out pyroptosis genes",
    subtitle = "Not used for clustering",
    x = NULL,
    y = "logCPM"
  )

ggsave(
  file.path(OUTDIR, "Fig1F_all84_heldout_gene_boxplots.png"),
  panel_F,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 8. Decision summary
# ============================================================

k2_metric <- consensus_metrics %>% filter(K == 2)

gsdmd_base <- heldout_stats %>%
  filter(gene == "GSDMD", model == "Base_adjusted_for_sex_and_cancer_type")

gsdmd_myeloid <- heldout_stats %>%
  filter(gene == "GSDMD", model == "Additional_adjustment_for_myeloid_score")

heldout_base_score <- signature_stats %>%
  filter(signature == "Heldout_pyroptosis", model == "base")

heldout_myeloid_score <- signature_stats %>%
  filter(signature == "Heldout_pyroptosis", model == "myeloid_adjusted")

decision <- case_when(
  nrow(heldout_myeloid_score) == 1 &&
    heldout_myeloid_score$adjusted_beta > 0 &&
    heldout_myeloid_score$P_value < 0.05 &&
    k2_metric$min_cluster_size >= 10 ~
      "SUPPORTED: all-84 K2 supports a pyroptosis-module-high state by held-out, myeloid-adjusted validation.",

  nrow(heldout_myeloid_score) == 1 &&
    heldout_myeloid_score$adjusted_beta > 0 &&
    heldout_myeloid_score$P_value < 0.10 ~
      "PARTIAL: all-84 K2 shows a pyroptosis-related signal, but support is not strong enough for a definitive PA/PQ claim.",

  TRUE ~
      "NOT SUPPORTED: all-84 K2 does not provide sufficient non-circular support for pyroptosis high/low states."
)

summary_lines <- c(
  "ALL-84 PYROPTOSIS MODULE EUCLIDEAN K2 MAIN ANALYSIS",
  "",
  paste0("Samples analyzed: ", nrow(labels)),
  paste0(
    "K2 group sizes: ",
    paste(names(table(labels$Pyro_module_group)), as.integer(table(labels$Pyro_module_group)), sep = "=", collapse = ";")
  ),
  paste0(
    "Sex balance: ",
    paste(capture.output(print(table(labels$Pyro_module_group, labels$sex))), collapse = " | ")
  ),
  paste0(
    "Cancer-type balance: ",
    paste(capture.output(print(table(labels$Pyro_module_group, labels$subgroup))), collapse = " | ")
  ),
  paste0(
    "Final K4 overlap: ",
    paste(capture.output(print(table(labels$Pyro_module_group, labels$adjusted_two_block_K4))), collapse = " | ")
  ),
  "",
  "K2 consensus metrics:",
  paste0("Mean within-cluster consensus: ", sprintf("%.3f", k2_metric$mean_within_cluster_consensus)),
  paste0("Mean between-cluster consensus: ", sprintf("%.3f", k2_metric$mean_between_cluster_consensus)),
  paste0("PAC 0.1-0.9: ", sprintf("%.3f", k2_metric$PAC_0.1_0.9)),
  paste0("Mean consensus silhouette: ", sprintf("%.3f", k2_metric$mean_consensus_silhouette)),
  "",
  paste0(
    "GSDMD base logFC/FDR: ",
    ifelse(nrow(gsdmd_base) == 1, paste0(sprintf("%.3f", gsdmd_base$adjusted_logFC), "/", signif(gsdmd_base$FDR, 3)), "not available")
  ),
  paste0(
    "GSDMD myeloid-adjusted logFC/FDR: ",
    ifelse(nrow(gsdmd_myeloid) == 1, paste0(sprintf("%.3f", gsdmd_myeloid$adjusted_logFC), "/", signif(gsdmd_myeloid$FDR, 3)), "not available")
  ),
  paste0(
    "Held-out pyroptosis base beta/FDR: ",
    ifelse(nrow(heldout_base_score) == 1, paste0(sprintf("%.3f", heldout_base_score$adjusted_beta), "/", signif(heldout_base_score$FDR, 3)), "not available")
  ),
  paste0(
    "Held-out pyroptosis myeloid-adjusted beta/P: ",
    ifelse(nrow(heldout_myeloid_score) == 1, paste0(sprintf("%.3f", heldout_myeloid_score$adjusted_beta), "/", signif(heldout_myeloid_score$P_value, 3)), "not available")
  ),
  "",
  decision
)

writeLines(
  summary_lines,
  file.path(OUTDIR, "all84_module_Euclidean_K2_final_decision_summary.txt")
)

cat("\n============================================================\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nOutput directory:", OUTDIR, "\n")
cat("============================================================\n")

# ============================================================
# 9. Assemble advisor preview figure
# ============================================================

read_panel <- function(path) {
  image_read(path) %>%
    image_trim() %>%
    image_background("white", flatten = TRUE) %>%
    image_border("white", "20x20")
}

fit_box <- function(img, width, height) {
  image_resize(img, paste0(width, "x", height)) %>%
    image_extent(paste0(width, "x", height), gravity = "center", color = "white")
}

label_panel <- function(img, lab) {
  image_annotate(
    img,
    lab,
    gravity = "northwest",
    location = "+18+10",
    size = 70,
    weight = 700,
    color = "black"
  )
}

single_w <- 1500
single_h <- 1050

A <- read_panel(file.path(OUTDIR, "Fig1A_all84_analysis_design.png")) %>% fit_box(single_w, single_h) %>% label_panel("A")
B <- read_panel(file.path(OUTDIR, "Fig1B_all84_K2_K4_consensus_metrics.png")) %>% fit_box(single_w, single_h) %>% label_panel("B")
C <- read_panel(file.path(OUTDIR, "Fig1C_all84_module_score_heatmap.png")) %>% fit_box(single_w * 2, single_h) %>% label_panel("C")
D <- read_panel(file.path(OUTDIR, "Fig1D_all84_K2_consensus_heatmap.png")) %>% fit_box(single_w, single_h) %>% label_panel("D")
E <- read_panel(file.path(OUTDIR, "Fig1E_all84_validation_signature_boxplots.png")) %>% fit_box(single_w, single_h) %>% label_panel("E")
F <- read_panel(file.path(OUTDIR, "Fig1F_all84_heldout_gene_boxplots.png")) %>% fit_box(single_w, single_h) %>% label_panel("F")

row1 <- image_append(c(A, B), stack = FALSE)
row2 <- image_append(c(C, D), stack = FALSE)
row3 <- image_append(c(E, F), stack = FALSE)

title_bar <- image_blank(width = single_w * 3, height = 180, color = "white") %>%
  image_annotate(
    "Fig. 1 preview | All-84 pyroptosis-module high/low stratification",
    gravity = "northwest",
    location = "+50+35",
    size = 58,
    weight = 700,
    color = "black"
  ) %>%
  image_annotate(
    "No sample prefiltering; K=2 prespecified as module-high/module-low; held-out validation reported",
    gravity = "northwest",
    location = "+52+110",
    size = 34,
    color = "#444444"
  )

final <- image_append(c(title_bar, row1, row2, row3), stack = TRUE) %>%
  image_border("white", "35x35")

image_write(final, file.path(OUTDIR, "Fig1_all84_pyro_module_main_preview.png"))
image_write(final, file.path(OUTDIR, "Fig1_all84_pyro_module_main_preview.pdf"))
