# Reproduce Fig. 3C-E and export the requested hub-gene matrices.
# Fig. 3A/B are reproduced independently by 02_fig3AB_original_RNA_part1.R.

require_packages(c("pheatmap", "ggplot2", "WGCNA", "scales"))

# Load matched metadata for Fig. 3E.
metadata <- utils::read.csv(data_file("metadata_WES_RNA_match.csv"), check.names = TRUE)
metadata$WES_group <- format_group(metadata$evolutionary_branch)

# Fig. 3C: reproduce purple-module hub selection from the saved WGCNA intermediate object.
wgcna_object <- readRDS(data_file("kmt2c_only_vs_no_egfr_kmt2c_WGCNA_intermediate.rds"))
dat_expr <- as.data.frame(wgcna_object$datExpr)
dat_traits <- as.data.frame(wgcna_object$datTraits)
module_colors <- wgcna_object$moduleColors
module_eigengenes <- as.data.frame(wgcna_object$MEs)
purple_index <- which(module_colors == "purple")

gene_significance <- stats::cor(dat_expr, dat_traits[, "KMT2C_mut"], use = "pairwise.complete.obs")
module_membership <- stats::cor(dat_expr, module_eigengenes, use = "pairwise.complete.obs")
purple_mm_column <- grep("purple$", colnames(module_membership), value = TRUE)
purple_statistics <- data.frame(
  gene = colnames(dat_expr)[purple_index],
  MM = module_membership[purple_index, purple_mm_column],
  GS = gene_significance[purple_index, 1],
  stringsAsFactors = FALSE
)
purple_statistics$MM_rank <- rank(-abs(purple_statistics$MM), ties.method = "first")
purple_statistics$GS_rank <- rank(-abs(purple_statistics$GS), ties.method = "first")
purple_statistics$sum_rank <- purple_statistics$MM_rank + purple_statistics$GS_rank
purple_statistics <- purple_statistics[order(purple_statistics$sum_rank), ]
hub_genes <- head(purple_statistics$gene, 20)
write_source_table(purple_statistics, "Fig3C_purple_module_MM_GS.csv")

# Scale each hub gene across samples and order samples by mutation status.
sample_order <- rownames(dat_traits)[order(dat_traits$KMT2C_mut, decreasing = TRUE)]
hub_raw <- t(as.matrix(dat_expr[sample_order, hub_genes, drop = FALSE]))
hub_z <- t(scale(as.matrix(dat_expr[, hub_genes, drop = FALSE])))
hub_z <- hub_z[, sample_order, drop = FALSE]
annotation <- data.frame(
  Truncal_mutation = ifelse(
    dat_traits[sample_order, "KMT2C_mut"] == 1,
    "KMT2C+ group", "EGFR- and KMT2C- group"
  ),
  row.names = sample_order
)

# Obtain the row order used by hierarchical clustering before exporting the plotted matrix.
heatmap_probe <- pheatmap::pheatmap(
  hub_z, cluster_rows = TRUE, cluster_cols = FALSE, silent = TRUE
)
visual_order <- rownames(hub_z)[heatmap_probe$tree_row$order]
hub_raw <- hub_raw[visual_order, , drop = FALSE]
hub_z <- hub_z[visual_order, , drop = FALSE]
write_source_matrix(hub_raw, "Fig3C_hub_gene_expression_raw.csv", "gene")
write_source_matrix(hub_z, "Fig3C_hub_gene_expression_zscore.csv", "gene")
write_source_table(
  data.frame(sample = sample_order, KMT2C_mut = dat_traits[sample_order, "KMT2C_mut"],
             group = annotation$Truncal_mutation),
  "Fig3C_sample_annotation.csv"
)
pheatmap::pheatmap(
  hub_z, cluster_rows = TRUE, cluster_cols = FALSE,
  annotation_col = annotation,
  annotation_colors = list(Truncal_mutation = c(
    `KMT2C+ group` = "#3C5488", `EGFR- and KMT2C- group` = "#D3D3D3"
  )),
  color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-4, 4, length.out = 101), show_colnames = FALSE,
  main = "Module purple hub genes (n = 20)", border_color = NA,
  filename = panel_file("Fig3C_purple_hub_gene_heatmap.png"),
  width = 10, height = 8
)

# Fig. 3D: plot the supplied GO Biological Process enrichment for the purple module.
fig3d <- utils::read.csv(
  data_file("kmt2c_only_vs_no_egfr_kmt2c_purple_GO_BP.csv"), check.names = FALSE
)
fig3d$GeneRatio_numeric <- ratio_to_numeric(fig3d$GeneRatio)
fig3d <- fig3d[order(fig3d$GeneRatio_numeric), ]
fig3d$Description <- factor(fig3d$Description, levels = fig3d$Description)
write_source_table(fig3d, "Fig3D_purple_module_GO_BP.csv")
p_fig3d <- ggplot2::ggplot(
  fig3d, ggplot2::aes(x = GeneRatio_numeric, y = Description)
) +
  ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
  ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35", limits = c(0, 0.05),
                                oob = scales::squish) +
  ggplot2::labs(x = "Gene Ratio", y = NULL, color = "Adjusted p-value") +
  theme_publication()
save_panel(p_fig3d, "Fig3D_purple_module_GO_BP.png", 7.0, 6.5)

# Fig. 3E: export primary-NE myeloid fractions and all legacy pairwise Welch tests.
fig3e <- metadata[
  metadata$Primary.Recurrent == "Primary" & metadata$Annotation == "NE" &
    metadata$WES_group %in% c("Group1.1", "Group1.2", "Group2.1", "Group2.2", "Group2.3"),
  c("RNAseq_ID", "WES_group", "Myeloid")
]
fig3e$WES_group <- factor(fig3e$WES_group,
                          levels = c("Group1.1", "Group1.2", "Group2.1", "Group2.2", "Group2.3"))
fig3e_tests <- pairwise_welch_tests(fig3e, "WES_group", "Myeloid")
write_source_table(fig3e, "Fig3E_primary_NE_myeloid_values.csv")
write_source_table(fig3e_tests, "Fig3E_pairwise_Welch_tests.csv")
p_fig3e <- cell_boxplot(
  fig3e, "WES_group", "Myeloid", "Myeloid",
  c(Group1.1 = "#2166AC", Group1.2 = "#67A9CF", Group2.1 = "#1B7837",
    Group2.2 = "#5AAE61", Group2.3 = "#A6DBA0")
)
save_panel(p_fig3e, "Fig3E_primary_NE_myeloid.png", 5.5, 5.0)
