# Reproduce Fig. 2A-C/F and Supplementary Fig. 3B-E.

require_packages(c("ggplot2", "pheatmap", "scales"))

# Define the categorical orders and color palettes used in the manuscript.
wes_levels <- c("Group1.1", "Group1.2", "Group2.1", "Group2.2", "Group2.3")
region_levels <- c("Primary_CE", "Primary_NE", "Recurrent_CE", "Recurrent_NE")
rna_levels <- c("GPM", "MTC", "NEU", "PPR", "Unclassified")
cell_types <- c("Oligodendrocytes", "Neuron", "Myeloid", "Malignant", "Endothelial", "Astrocyte")
wes_palette <- c(
  Group1.1 = "#2166AC", Group1.2 = "#67A9CF", Group2.1 = "#1B7837",
  Group2.2 = "#5AAE61", Group2.3 = "#A6DBA0"
)
region_palette <- c(
  Primary_CE = "#084594", Primary_NE = "#A94F1D",
  Recurrent_CE = "#9ECAE1", Recurrent_NE = "#FE9929"
)
rna_palette <- c(
  GPM = "#3C5488", MTC = "#4DBBD5", NEU = "#E64B35",
  PPR = "#00A087", Unclassified = "#DADADA"
)
cell_palette <- c(
  Oligodendrocytes = "#009E73", Neuron = "#56B4E9", Myeloid = "#E69F00",
  Malignant = "#B81830", Endothelial = "#888888", Astrocyte = "#CC79A7"
)

# Load the WES-matched metadata used for evolutionary-group panels.
metadata <- utils::read.csv(data_file("metadata_WES_RNA_match.csv"), check.names = TRUE)
metadata$WES_group <- factor(format_group(metadata$evolutionary_branch), levels = wes_levels)
metadata$Clinical_region <- factor(
  paste(metadata$Primary.Recurrent, metadata$Annotation, sep = "_"),
  levels = region_levels
)
metadata <- metadata[!is.na(metadata$WES_group), ]

# Fig. 2A: calculate clinical-region counts and proportions within each WES group.
fig2a_counts <- as.data.frame(xtabs(~ WES_group + Clinical_region, data = metadata))
names(fig2a_counts) <- c("WES_group", "Clinical_region", "count")
fig2a_counts <- fig2a_counts[fig2a_counts$count > 0, ]
fig2a_counts$proportion <- ave(
  fig2a_counts$count, fig2a_counts$WES_group,
  FUN = function(x) x / sum(x)
)
write_source_table(fig2a_counts, "Fig2A_clinical_region_distribution.csv")

p_fig2a <- ggplot2::ggplot(
  fig2a_counts,
  ggplot2::aes(x = WES_group, y = proportion, fill = Clinical_region)
) +
  ggplot2::geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(proportion, accuracy = 1)),
    position = ggplot2::position_stack(vjust = 0.5), color = "white", size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = region_palette) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = "WES cluster", y = "Sample distribution", fill = "Clinical region") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_fig2a, "Fig2A_clinical_region_distribution.png", 5.8, 5.2)

# Fig. 2B: join RNA-subtype assignments to WES groups and calculate proportions.
rna_assignments <- utils::read.csv(data_file("cluster_assignments_RNA_300.csv"), check.names = FALSE)
rna_assignments <- rna_assignments[, c("RNAseq_ID", "Cluster")]
rna_assignments$RNA_subtype <- factor(rna_assignments$Cluster, levels = rna_levels)
rna_wes <- merge(
  metadata[, c("RNAseq_ID", "WES_group")],
  rna_assignments[, c("RNAseq_ID", "RNA_subtype")],
  by = "RNAseq_ID"
)
fig2b_counts <- as.data.frame(xtabs(~ WES_group + RNA_subtype, data = rna_wes))
names(fig2b_counts) <- c("WES_group", "RNA_subtype", "count")
fig2b_counts <- fig2b_counts[fig2b_counts$count > 0, ]
fig2b_counts$proportion <- ave(
  fig2b_counts$count, fig2b_counts$WES_group,
  FUN = function(x) x / sum(x)
)
write_source_table(fig2b_counts, "Fig2B_RNA_subtype_distribution.csv")

p_fig2b <- ggplot2::ggplot(
  fig2b_counts,
  ggplot2::aes(x = WES_group, y = proportion, fill = RNA_subtype)
) +
  ggplot2::geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(proportion >= 0.06, scales::percent(proportion, accuracy = 0.1), "")),
    position = ggplot2::position_stack(vjust = 0.5), color = "white", size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = rna_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = "WES cluster", y = "Sample distribution", fill = "RNA subtype") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_fig2b, "Fig2B_RNA_subtype_distribution.png", 5.8, 5.2)

# Fig. 2C: reproduce the displayed row-scaled evolutionary expression heatmap.
fig2c_expression <- utils::read.csv(
  data_file("Fig.2C.plot_matrix.csv"), row.names = 1, check.names = FALSE
)
fig2c_z <- t(scale(t(as.matrix(fig2c_expression))))
write_source_matrix(fig2c_expression, "Fig2C_expression_voom.csv", "gene")
write_source_matrix(fig2c_z, "Fig2C_expression_zscore.csv", "gene")

heatmap_meta <- utils::read.csv(
  data_file("evolutionary_branches_WES_RNA_match.csv"), check.names = FALSE
)
rownames(heatmap_meta) <- heatmap_meta$RNAseq_ID
heatmap_meta <- heatmap_meta[colnames(fig2c_expression), , drop = FALSE]
fig2c_annotation <- data.frame(
  Subtype = heatmap_meta$Cluster,
  `MRI Region` = heatmap_meta$Annotation,
  `Pathology Status` = heatmap_meta$Primary.Recurrent,
  `Evolution Group` = format_group(heatmap_meta$cluster),
  check.names = FALSE,
  row.names = colnames(fig2c_expression)
)
annotation_colors <- list(
  Subtype = rna_palette,
  `MRI Region` = c(CE = "#0072B2", NE = "#D55E00"),
  `Pathology Status` = c(Primary = "#004D40", Recurrent = "#80CBC4"),
  `Evolution Group` = wes_palette
)
pheatmap::pheatmap(
  fig2c_z, cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = fig2c_annotation, annotation_colors = annotation_colors,
  color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101), show_colnames = FALSE,
  border_color = NA, fontsize_row = 8,
  filename = panel_file("Fig2C_evolutionary_expression_heatmap.png"),
  width = 8.0, height = 6.5
)

# Fig. 2F: average the six deconvolved cell fractions within each WES group.
fig2f_means <- stats::aggregate(
  metadata[, cell_types], list(WES_group = metadata$WES_group), mean, na.rm = TRUE
)
fig2f_long <- reshape(
  fig2f_means, varying = cell_types, v.names = "mean_proportion",
  timevar = "Cell_type", times = cell_types, direction = "long"
)
rownames(fig2f_long) <- NULL
write_source_table(fig2f_long[, c("WES_group", "Cell_type", "mean_proportion")],
                   "Fig2F_mean_cell_composition.csv")
p_fig2f <- ggplot2::ggplot(
  fig2f_long,
  ggplot2::aes(x = WES_group, y = mean_proportion, fill = Cell_type)
) +
  ggplot2::geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  ggplot2::scale_fill_manual(values = cell_palette) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Cell composition", fill = "Cell type") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_fig2f, "Fig2F_mean_cell_composition.png", 6.6, 5.2)

# Supplementary Fig. 3B: calculate clinical-region enrichment within WES groups.
region_counts <- as.data.frame(xtabs(~ WES_group + Clinical_region, data = metadata))
names(region_counts) <- c("WES_group", "Clinical_region", "count")
region_enrichment <- calculate_enrichment(region_counts, "WES_group", "Clinical_region")
write_source_table(region_enrichment, "FigS3B_region_enrichment_long.csv")
region_matrix <- enrichment_matrix(
  region_enrichment, "WES_group", "Clinical_region", wes_levels, region_levels
)
write_source_matrix(region_matrix, "FigS3B_region_enrichment_matrix.csv", "WES_group")
pheatmap::pheatmap(
  region_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
  color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101), border_color = NA,
  main = "Region Preference Within WES Subtypes log2(O/E)",
  filename = panel_file("FigS3B_region_enrichment.png"), width = 7.0, height = 5.5
)

# Supplementary Fig. 3C: calculate RNA-subtype enrichment within WES groups.
# Count only combinations observed in the data, matching the original script.
# Unobserved combinations therefore remain 0 (the white midpoint) in the matrix
# instead of being treated as count-zero observations with negative enrichment.
rna_counts <- stats::aggregate(
  list(count = rep.int(1L, nrow(rna_wes))),
  by = list(
    WES_group = as.character(rna_wes$WES_group),
    RNA_subtype = as.character(rna_wes$RNA_subtype)
  ),
  FUN = sum
)
rna_enrichment <- calculate_enrichment(rna_counts, "WES_group", "RNA_subtype")
write_source_table(rna_enrichment, "FigS3C_RNA_subtype_enrichment_long.csv")
present_rna_levels <- rna_levels[rna_levels %in% rna_enrichment$RNA_subtype]
rna_matrix <- enrichment_matrix(
  rna_enrichment, "WES_group", "RNA_subtype", wes_levels, present_rna_levels
)
write_source_matrix(rna_matrix, "FigS3C_RNA_subtype_enrichment_matrix.csv", "WES_group")
pheatmap::pheatmap(
  rna_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
  color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101), border_color = NA,
  main = "RNA Cluster Preference Within WES Subtypes log2(O/E)",
  filename = panel_file("FigS3C_RNA_subtype_enrichment.png"), width = 7.0, height = 5.5
)

# Supplementary Fig. 3D: primary/recurrent proportions within WES groups.
figs3d_counts <- as.data.frame(xtabs(~ WES_group + Primary.Recurrent, data = metadata))
names(figs3d_counts) <- c("WES_group", "Clinical_status", "count")
figs3d_counts$proportion <- ave(
  figs3d_counts$count, figs3d_counts$WES_group, FUN = function(x) x / sum(x)
)
write_source_table(figs3d_counts, "FigS3D_primary_recurrent_distribution.csv")
p_figs3d <- ggplot2::ggplot(
  figs3d_counts,
  ggplot2::aes(x = WES_group, y = proportion, fill = Clinical_status)
) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(proportion, accuracy = 1)),
    position = ggplot2::position_stack(vjust = 0.5), color = "white"
  ) +
  ggplot2::scale_fill_manual(values = c(Primary = "#236A5D", Recurrent = "#8DD3C7")) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = "WES cluster", y = "Sample distribution", fill = "Status") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_figs3d, "FigS3D_primary_recurrent_distribution.png", 5.8, 5.2)

# Supplementary Fig. 3E: CE/NE proportions within WES groups.
figs3e_counts <- as.data.frame(xtabs(~ WES_group + Annotation, data = metadata))
names(figs3e_counts) <- c("WES_group", "MRI_region", "count")
figs3e_counts$proportion <- ave(
  figs3e_counts$count, figs3e_counts$WES_group, FUN = function(x) x / sum(x)
)
write_source_table(figs3e_counts, "FigS3E_CE_NE_distribution.csv")
p_figs3e <- ggplot2::ggplot(
  figs3e_counts,
  ggplot2::aes(x = WES_group, y = proportion, fill = MRI_region)
) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(proportion, accuracy = 1)),
    position = ggplot2::position_stack(vjust = 0.5), color = "white"
  ) +
  ggplot2::scale_fill_manual(values = c(CE = "#0072B2", NE = "#D55E00")) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = "WES cluster", y = "Sample distribution", fill = "MRI region") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_figs3e, "FigS3E_CE_NE_distribution.png", 5.8, 5.2)
