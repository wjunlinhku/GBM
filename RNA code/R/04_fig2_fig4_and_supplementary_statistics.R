# Reproduce differential-expression, enrichment, clinical-composition, and cell-composition panels.

require_packages(c("ggplot2", "ggrepel", "scales", "cowplot"))

# Shared palettes and helper functions for volcano and enrichment panels.
region_levels <- c("Primary_CE", "Primary_NE", "Recurrent_CE", "Recurrent_NE")
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

# Draw a volcano plot using the significance labels supplied in each DEG table.
make_volcano <- function(data, title, output_name) {
  plot_data <- data[!is.na(data$padj), ]
  significant <- plot_data[plot_data$Significance != "NotSig", ]
  labels <- head(significant[order(significant$padj), ], 10)
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)), color = Significance)
  ) +
    ggplot2::geom_point(alpha = 0.75, size = 1.4) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666") +
    ggplot2::scale_color_manual(values = c(Up = "#E64B35", Down = "#3182BD", NotSig = "#BDBDBD")) +
    ggrepel::geom_text_repel(
      data = labels, ggplot2::aes(label = Gene), color = "black",
      size = 3.2, fontface = "italic", max.overlaps = 20
    ) +
    ggplot2::labs(x = "Log2 fold change", y = "-Log10 adjusted P value", title = title) +
    theme_publication()
  save_panel(p, output_name, 6.0, 6.0)
}

# Combine up- and down-regulated enrichment tables into signed -log10(FDR) values.
prepare_mirror_data <- function(up, down, pathways = NULL) {
  up$Direction <- "Up-regulated"
  down$Direction <- "Down-regulated"
  combined <- rbind(up, down)
  if (!is.null(pathways)) {
    combined <- combined[combined$Description %in% pathways | combined$ID %in% pathways, ]
  }
  combined$minus_log10_FDR <- -log10(pmax(combined$p.adjust, .Machine$double.xmin))
  combined$signed_minus_log10_FDR <- ifelse(
    combined$Direction == "Up-regulated",
    combined$minus_log10_FDR, -combined$minus_log10_FDR
  )
  combined
}

# Draw a horizontal signed enrichment bar chart.
make_mirror_plot <- function(data, title, output_name, width, height) {
  ordered <- data[order(data$signed_minus_log10_FDR), ]
  ordered$Description <- factor(ordered$Description, levels = ordered$Description)
  p <- ggplot2::ggplot(
    ordered,
    ggplot2::aes(x = signed_minus_log10_FDR, y = Description, fill = Direction)
  ) +
    ggplot2::geom_col(color = "#444444", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = c(`Up-regulated` = "#E64B35", `Down-regulated` = "#4DBBD5")) +
    ggplot2::labs(x = "Signed -Log10 adjusted P value", y = NULL, title = title) +
    theme_publication() +
    ggplot2::theme(legend.position = "top", axis.text.y = ggplot2::element_text(size = 8))
  save_panel(p, output_name, width, height)
}

# Fig. 2D: EGFR-amplified versus non-EGFR differential expression.
fig2d <- utils::read.csv(data_file("EGFR vs none_EGFR.csv"), check.names = FALSE)
write_source_table(fig2d, "Fig2D_EGFR_vs_nonEGFR_DEG.csv")
make_volcano(fig2d, "EGFR versus non-EGFR", "Fig2D_EGFR_vs_nonEGFR_volcano.png")

# Fig. 2E: KEGG enrichment of genes upregulated in EGFR-amplified samples.
fig2e <- utils::read.csv(
  data_file("KEGG_ALL_none_EGFR vs EGFR up.csv"), check.names = FALSE
)
fig2e$GeneRatio_numeric <- ratio_to_numeric(fig2e$GeneRatio)
fig2e <- fig2e[order(fig2e$GeneRatio_numeric), ]
fig2e$Description <- factor(fig2e$Description, levels = fig2e$Description)
write_source_table(fig2e, "Fig2E_EGFR_upregulated_KEGG.csv")
p_fig2e <- ggplot2::ggplot(fig2e, ggplot2::aes(x = GeneRatio_numeric, y = Description)) +
  ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
  ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35", limits = c(0, 0.05),
                                oob = scales::squish) +
  ggplot2::labs(x = "Gene Ratio", y = NULL, color = "Adjusted p-value") +
  theme_publication()
save_panel(p_fig2e, "Fig2E_EGFR_upregulated_KEGG.png", 7.0, 6.5)

# Load clinical subtype assignments and deconvolution fractions for Fig. 4A-B/D and Fig. S2C.
clinical <- utils::read.csv(data_file("cluster_assignments_RNA_300.csv"), check.names = FALSE)
clinical$Clinical_region <- factor(
  paste(clinical$Primary.Recurrent, clinical$Annotation, sep = "_"),
  levels = region_levels
)

# Fig. 4A: RNA subtype proportions across the four clinical regions.
fig4a_counts <- as.data.frame(xtabs(~ Clinical_region + GPM_cluster_by_sample, data = clinical))
names(fig4a_counts) <- c("Clinical_region", "RNA_subtype", "count")
fig4a_counts <- fig4a_counts[fig4a_counts$count > 0, ]
fig4a_counts$proportion <- ave(
  fig4a_counts$count, fig4a_counts$Clinical_region, FUN = function(x) x / sum(x)
)
write_source_table(fig4a_counts, "Fig4A_RNA_subtype_by_clinical_region.csv")
p_fig4a <- ggplot2::ggplot(
  fig4a_counts,
  ggplot2::aes(x = Clinical_region, y = proportion, fill = RNA_subtype)
) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::scale_fill_manual(values = rna_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Proportion", fill = "Subtype") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_fig4a, "Fig4A_RNA_subtype_by_clinical_region.png", 5.5, 5.2)

# Merge deconvolution values after normalizing R-converted sample identifiers.
deconvolution <- utils::read.csv(data_file("deconv_result_RNA_300.csv"), check.names = TRUE)
sample_column <- if ("X" %in% names(deconvolution)) "X" else names(deconvolution)[2]
deconvolution[[sample_column]] <- gsub("\\.", "-", deconvolution[[sample_column]])
clinical_deconv <- merge(
  clinical, deconvolution,
  by.x = "RNAseq_ID", by.y = sample_column
)
cell_types <- c("Oligodendrocytes", "Neuron", "Myeloid", "Malignant", "Endothelial", "Astrocyte")

# Fig. 4B: average deconvolved composition in each clinical region.
fig4b_means <- stats::aggregate(
  clinical_deconv[, cell_types],
  list(Clinical_region = clinical_deconv$Clinical_region), mean, na.rm = TRUE
)
fig4b_long <- reshape(
  fig4b_means, varying = cell_types, v.names = "mean_proportion",
  timevar = "Cell_type", times = cell_types, direction = "long"
)
rownames(fig4b_long) <- NULL
fig4b_long <- fig4b_long[, c("Clinical_region", "Cell_type", "mean_proportion")]
write_source_table(fig4b_long, "Fig4B_mean_cell_composition.csv")
p_fig4b <- ggplot2::ggplot(
  fig4b_long,
  ggplot2::aes(x = Clinical_region, y = mean_proportion, fill = Cell_type)
) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::scale_fill_manual(values = cell_palette) +
  ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Cell composition", fill = "Cell type") +
  theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_panel(p_fig4b, "Fig4B_mean_cell_composition.png", 6.0, 5.2)

# Fig. 4C: recurrent-NE versus primary-NE differential expression.
fig4c <- utils::read.csv(data_file("Recurrent_NE vs Primary_NE.csv"), check.names = FALSE)
write_source_table(fig4c, "Fig4C_Recurrent_NE_vs_Primary_NE_DEG.csv")
make_volcano(fig4c, "Recurrent NE versus Primary NE", "Fig4C_Recurrent_NE_vs_Primary_NE_volcano.png")

# Fig. 4D and Supplementary Fig. 2C: sample-level fractions and Welch tests.
ne_data <- clinical_deconv[
  clinical_deconv$Clinical_region %in% c("Primary_NE", "Recurrent_NE"),
  c("RNAseq_ID", "Clinical_region", "Myeloid", "Malignant", "Oligodendrocytes")
]
ne_data$Clinical_region <- droplevels(ne_data$Clinical_region)
write_source_table(ne_data[, c("RNAseq_ID", "Clinical_region", "Myeloid", "Malignant")],
                   "Fig4D_NE_cell_fraction_values.csv")
write_source_table(pairwise_welch_tests(ne_data, "Clinical_region", "Myeloid"),
                   "Fig4D_Myeloid_Welch_test.csv")
write_source_table(pairwise_welch_tests(ne_data, "Clinical_region", "Malignant"),
                   "Fig4D_Malignant_Welch_test.csv")
write_source_table(ne_data[, c("RNAseq_ID", "Clinical_region", "Oligodendrocytes")],
                   "FigS2C_Oligodendrocyte_values.csv")
write_source_table(pairwise_welch_tests(ne_data, "Clinical_region", "Oligodendrocytes"),
                   "FigS2C_Oligodendrocyte_Welch_test.csv")

p_myeloid <- cell_boxplot(ne_data, "Clinical_region", "Myeloid", "Myeloid", region_palette)
p_malignant <- cell_boxplot(ne_data, "Clinical_region", "Malignant", "Malignant", region_palette)
save_panel(p_myeloid, "Fig4D_Myeloid.png", 4.6, 4.6)
save_panel(p_malignant, "Fig4D_Malignant.png", 4.6, 4.6)
p_figs2c <- cell_boxplot(
  ne_data, "Clinical_region", "Oligodendrocytes", "Oligodendrocytes", region_palette
)
save_panel(p_figs2c, "FigS2C_Oligodendrocytes.png", 4.8, 4.8)

# Fig. 4E and Supplementary Fig. 5A: KEGG pathways for recurrent NE versus primary NE.
kegg_up <- utils::read.csv(
  data_file("KEGG_ALL_Recurrent_NE vs Primary_NE_up.csv"), check.names = FALSE
)
kegg_down <- utils::read.csv(
  data_file("KEGG_ALL_Recurrent_NE vs Primary_NE_down.csv"), check.names = FALSE
)
selected_pathways <- c(
  "Cytokine-cytokine receptor interaction", "Chemokine signaling pathway",
  "TNF signaling pathway", "NF-kappa B signaling pathway", "Hematopoietic cell lineage",
  "Neuroactive ligand-receptor interaction", "Neuroactive ligand signaling",
  "Calcium signaling pathway", "Serotonergic synapse", "Dopaminergic synapse",
  "Cell cycle", "Notch signaling pathway", "Homologous recombination",
  "Fanconi anemia pathway", "p53 signaling pathway"
)
fig4e <- prepare_mirror_data(kegg_up, kegg_down, selected_pathways)
write_source_table(fig4e, "Fig4E_selected_KEGG_mirror.csv")
make_mirror_plot(fig4e, NULL, "Fig4E_selected_KEGG_mirror.png", 7.0, 6.0)

figs5a <- prepare_mirror_data(kegg_up, kegg_down)
write_source_table(figs5a, "FigS5A_all_KEGG_mirror.csv")
make_mirror_plot(figs5a, NULL, "FigS5A_all_KEGG_mirror.png", 8.0, 11.0)

# Supplementary Fig. 3A: split the eight supplied GO terms into BP and MF panels.
figs3a <- utils::read.csv(
  data_file("GO_ALL_none_EGFR vs EGFR down.csv"), check.names = FALSE
)
figs3a$GeneRatio_numeric <- ratio_to_numeric(figs3a$GeneRatio)
write_source_table(figs3a, "FigS3A_EGFR_downregulated_GO.csv")
make_go_panel <- function(data, title) {
  data <- data[order(data$GeneRatio_numeric), ]
  data$Description <- factor(data$Description, levels = data$Description)
  ggplot2::ggplot(data, ggplot2::aes(x = GeneRatio_numeric, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35", limits = c(0, 0.05),
                                  oob = scales::squish) +
    ggplot2::labs(x = "Gene Ratio", y = NULL, title = title, color = "Adjusted p-value") +
    theme_publication()
}
p_bp <- make_go_panel(figs3a[figs3a$ONTOLOGY == "BP", ], "BP")
p_mf <- make_go_panel(figs3a[figs3a$ONTOLOGY == "MF", ], "MF")
p_figs3a <- cowplot::plot_grid(p_bp, p_mf, nrow = 1, align = "h")
save_panel(p_figs3a, "FigS3A_EGFR_downregulated_GO.png", 13.0, 6.2)
