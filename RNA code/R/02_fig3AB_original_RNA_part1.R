# Reproduce Fig. 3A and Fig. 3B from the original RNA_part1.R logic.
#
# Input/output routing uses the package-relative helpers from 00_common.R. The
# preprocessing, ssGSEA, correlation, significance, and pheatmap code below is
# copied from RNA_part1.R without changing its analytical logic.

require_packages(c(
  "limma", "edgeR", "HGNChelper", "pheatmap", "msigdbr", "GSVA", "dplyr"
))

library(limma)
library(edgeR)
library(HGNChelper)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(dplyr)

# ==============================================================================
# 1. DATA LOADING AND PREPROCESSING (copied from RNA_part1.R)
# ==============================================================================

# Load bulk RNA-seq count matrix
count_data <- read.csv(data_file("corrected_counts_tow_batch.csv"), row.names = 1)
colnames(count_data) <- gsub('\\.', '-', colnames(count_data))

# Load and filter sample metadata
sample_metadata <- read.csv(data_file("RNA_WES_match3.csv"))
sample_metadata <- sample_metadata[!is.na(sample_metadata$RNAseq_ID), ]
sample_metadata <- sample_metadata[sample_metadata$Primary.Recurrent %in% c('Primary', 'Recurrent'), ]
sample_metadata <- sample_metadata[sample_metadata$IDH %in% c('WT', 'Wt', 'wt', 'wild type'), ]
sample_metadata <- sample_metadata[sample_metadata$Annotation %in% c('CE', 'NE'), ]
sample_metadata <- sample_metadata[sample_metadata$RNAseq_ID != 'R151J', ]
rownames(sample_metadata) <- sample_metadata[, 1]
sample_metadata <- sample_metadata[!is.na(sample_metadata$WES_ID), ]
sample_metadata$patient_id2 <- sub("^([A-Za-z])([0-9]+).*", "\\1\\2", sample_metadata$WES_ID)

# Load evolutionary branches
# The submitted CSV currently carries a UTF-8 BOM; declaring it explicitly
# restores the original column name `patient_id` on this Windows R runtime.
evolutionary_branches <- read.csv(data_file("patient_group.csv"), fileEncoding = "UTF-8-BOM")
evolutionary_branches <- merge(evolutionary_branches, sample_metadata, by.x = 'patient_id', by.y = 'patient_id2')

# Intersect samples across datasets
shared_samples <- intersect(colnames(count_data), evolutionary_branches$RNAseq_ID)
count_data <- count_data[, shared_samples]
sample_metadata <- sample_metadata[shared_samples, ]

# Gene symbol standardization
check <- checkGeneSymbols(rownames(count_data), species = "human")
check <- na.omit(check)
count_data <- count_data[check$x, ]
expr_matrix <- count_data

# Load deconvolution results
deconv_result <- read.csv(data_file("theta_type_with_survival.csv"), header = TRUE)
deconv_result$X <- gsub('\\.', '-', deconv_result$X)

# Load cluster assignments
cluster_assignments <- read.csv(data_file("cluster_assignments.csv"))
evolutionary_branches <- merge(evolutionary_branches, cluster_assignments[, c('RNAseq_ID', 'GPM_cluster_by_sample')], by = 'RNAseq_ID')
colnames(evolutionary_branches)[14] <- 'Cluster'

# Final formatting and sorting
evolutionary_branches <- evolutionary_branches[evolutionary_branches$RNAseq_ID %in% colnames(expr_matrix), ]
rownames(evolutionary_branches) <- evolutionary_branches$RNAseq_ID
evolutionary_branches <- evolutionary_branches[sample_metadata$RNAseq_ID, ]
deconv_result <- deconv_result[deconv_result$X %in% colnames(expr_matrix), ]
evolutionary_branches$cluster=gsub('group','Group',evolutionary_branches$cluster)

# Filter out low expression genes
filter_low_expression <- function(expr_matrix, min_count = 10, min_samples = 0.1) {
  keep <- rowSums(expr_matrix > min_count) >= ncol(expr_matrix) * min_samples
  return(expr_matrix[keep, ])
}

# Integrate bulk RNA-seq data with evolutionary branches (voom + PCA)
integrate_bulk_evolution <- function(expr_matrix, metadata, evolutionary_branches, deconv_result) {
  # Create DGEList and calculate normalization factors
  dge <- DGEList(counts = expr_matrix)
  dge <- calcNormFactors(dge)

  # Voom transformation
  design <- model.matrix(~ evolutionary_branches$cluster)
  v <- voom(dge, design, plot = FALSE)

  # PCA analysis
  pca_result <- prcomp(t(v$E), scale. = TRUE)

  # Combine metadata
  result_list <- list(
    expr_matrix = v$E,
    pca = pca_result,
    metadata = cbind(metadata,
                     evolutionary_branch = evolutionary_branches$cluster,
                     deconv_result)
  )
  return(result_list)
}

# Calculate pathway activity (GSVA) and correlate with TME composition
analyze_pathway_activity <- function(expr_matrix, metadata, title) {
  # Extract pathways
  kegg_df <- msigdbr(species = "Homo sapiens", category = "C2")
  cell_types <- colnames(metadata)[15:20]

  pathways <- list(
    MAPK_signaling = "KEGG_MAPK_SIGNALING_PATHWAY",
    KRAS = "CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_UP",
    RTK_signaling = "KEGG_ERBB_SIGNALING_PATHWAY",
    PI3K_AKT = "REACTOME_PI3K_AKT_ACTIVATION",
    Apoptosis = "KEGG_APOPTOSIS",
    Immune_Invasion = "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"
  )

  pathway_genes <- lapply(pathways, function(p) {
    kegg_df %>% filter(gs_name == p) %>% pull(gene_symbol)
  })

  # GSVA scoring
  scores <- gsva(expr = as.matrix(expr_matrix), gset.idx.list = pathway_genes, kcdf = "Poisson", method = "ssgsea")
  pathway_scores <- t(scores)

  metadata_with_scores <- cbind(metadata, pathway_scores[metadata$RNAseq_ID, ])

  # Spearman correlation
  pathway_cell_correlations <- cor(pathway_scores, metadata[, cell_types], use = "complete.obs")
  p_value_matrix <- matrix(NA, nrow = nrow(pathway_cell_correlations), ncol = ncol(pathway_cell_correlations))

  for (i in 1:nrow(pathway_cell_correlations)) {
    for (j in 1:ncol(pathway_cell_correlations)) {
      cor_test <- cor.test(pathway_scores[, i], metadata[, cell_types][, j], method = "spearman", exact = FALSE)
      p_value_matrix[i, j] <- cor_test$p.value
    }
  }

  # Significance matrix for heatmap
  significance_matrix <- matrix("", nrow = nrow(p_value_matrix), ncol = ncol(p_value_matrix))
  significance_matrix[p_value_matrix < 0.05] <- "*"
  significance_matrix[p_value_matrix < 0.01] <- "**"
  significance_matrix[p_value_matrix < 0.001] <- "***"

  color_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

  # Save heatmap
  pheatmap(
    as.matrix(pathway_cell_correlations),
    color = color_palette,
    breaks = seq(-1, 1, length.out = length(color_palette) + 1),
    cluster_rows = FALSE, cluster_cols = FALSE,
    display_numbers = significance_matrix, number_color = "black",
    border_color = NA, cellwidth = 50, cellheight = 50,
    main = title, angle_col = 45, fontfamily = "Arial",
    fontsize = 13, fontsize_row = 12, fontsize_col = 12, fontsize_number = 10,
    filename = paste0(title, ".png"), width = 6.67, height = 6
  )

  return(list(metadata = metadata_with_scores, pathway_scores = pathway_scores, correlations = pathway_cell_correlations))
}

# ==============================================================================
# 2. ORIGINAL PREPROCESSING AND THE TWO SUBGROUP CALLS USED FOR FIG. 3A/B
# ==============================================================================

expr_matrix_filtered <- filter_low_expression(expr_matrix)
bulk_data <- integrate_bulk_evolution(expr_matrix_filtered, sample_metadata, evolutionary_branches, deconv_result)

run_original_fig3_panel <- function(wes_group, panel) {
  cluster_metadata <- bulk_data$metadata[
    bulk_data$metadata$evolutionary_branch == wes_group &
      bulk_data$metadata$Primary.Recurrent == "Primary" &
      bulk_data$metadata$Annotation == 'NE',
  ]
  cluster_mat <- bulk_data$expr_matrix[, cluster_metadata$RNAseq_ID]
  title <- paste0(wes_group, " Primary NE")
  result <- analyze_pathway_activity(cluster_mat, cluster_metadata, title)

  # Export the exact matrices produced by the original calculation without
  # changing the function or the plotted significance symbols.
  cell_types <- colnames(cluster_metadata)[15:20]
  p_value_matrix <- matrix(
    NA, nrow = nrow(result$correlations), ncol = ncol(result$correlations),
    dimnames = dimnames(result$correlations)
  )
  complete_n_matrix <- p_value_matrix
  for (i in 1:nrow(result$correlations)) {
    for (j in 1:ncol(result$correlations)) {
      complete <- complete.cases(result$pathway_scores[, i], cluster_metadata[, cell_types][, j])
      complete_n_matrix[i, j] <- sum(complete)
      cor_test <- cor.test(
        result$pathway_scores[, i], cluster_metadata[, cell_types][, j],
        method = "spearman", exact = FALSE
      )
      p_value_matrix[i, j] <- cor_test$p.value
    }
  }
  significance_matrix <- matrix(
    "", nrow = nrow(p_value_matrix), ncol = ncol(p_value_matrix),
    dimnames = dimnames(p_value_matrix)
  )
  significance_matrix[p_value_matrix < 0.05] <- "*"
  significance_matrix[p_value_matrix < 0.01] <- "**"
  significance_matrix[p_value_matrix < 0.001] <- "***"

  write_source_matrix(result$correlations, paste0(panel, "_correlation_matrix.csv"), "pathway")
  write_source_matrix(p_value_matrix, paste0(panel, "_pvalue_matrix.csv"), "pathway")
  write_source_matrix(complete_n_matrix, paste0(panel, "_complete_n_matrix.csv"), "pathway")
  write_source_matrix(t(result$pathway_scores), paste0(panel, "_ssGSEA_scores.csv"), "pathway")
  write_source_matrix(significance_matrix, paste0(panel, "_significance_symbols.csv"), "pathway")

  original_png <- file.path(getwd(), paste0(title, ".png"))
  standard_png <- panel_file(paste0(panel, "_pathway_cell_correlation.png"))
  if (!file.copy(original_png, standard_png, overwrite = TRUE)) {
    stop("Could not copy the original panel output to: ", standard_png)
  }
  invisible(result)
}

previous_working_directory <- setwd(output_file("panels"))
fig3ab_results <- tryCatch(
  list(
    Fig3A = run_original_fig3_panel("Group1.1", "Fig3A"),
    Fig3B = run_original_fig3_panel("Group2.2", "Fig3B")
  ),
  finally = setwd(previous_working_directory)
)
