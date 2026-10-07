# Reproduce Fig. 3F-G and derive Table 2 from enhancer targets and RNA-seq DEGs.

require_packages(c("ggplot2", "VennDiagram", "scales", "grid"))

# Fig. 3F: retain the 20 GO terms explicitly selected in the original plotting script.
target_terms <- c(
  "extracellular matrix organization",
  "extracellular structure organization",
  "external encapsulating structure organization",
  "camera-type eye development",
  "eye development",
  "visual system development",
  "sensory system development",
  "eye morphogenesis",
  "sensory organ morphogenesis",
  "regulation of transmembrane receptor protein serine/threonine kinase signaling pathway",
  "regulation of morphogenesis of an epithelium",
  "proteoglycan metabolic process",
  "negative regulation of cellular response to growth factor stimulus",
  "camera-type eye morphogenesis",
  "lateral sprouting from an epithelium",
  "morphogenesis of an epithelial bud",
  "morphogenesis of an epithelial fold",
  "prostate gland epithelium morphogenesis",
  "regulation of animal organ formation",
  "prostate gland morphogenesis"
)
fig3f_all <- utils::read.csv(data_file("GO_ALL_mouse gain_enhancers2.csv"), check.names = FALSE)
fig3f <- fig3f_all[fig3f_all$Description %in% target_terms, , drop = FALSE]
fig3f$GeneRatio_numeric <- ratio_to_numeric(fig3f$GeneRatio)
fig3f <- fig3f[order(fig3f$GeneRatio_numeric), ]
fig3f$Description <- factor(fig3f$Description, levels = fig3f$Description)
write_source_table(fig3f, "Fig3F_gain_enhancer_GO_selected.csv")

p_fig3f <- ggplot2::ggplot(
  fig3f, ggplot2::aes(x = GeneRatio_numeric, y = Description)
) +
  ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
  ggplot2::scale_color_gradient(
    low = "#4DBBD5", high = "#E64B35", limits = c(0, 0.05),
    oob = scales::squish, name = "Adjusted p-value"
  ) +
  ggplot2::scale_size(range = c(2, 7)) +
  ggplot2::labs(x = "Gene Ratio", y = NULL) +
  theme_publication()
save_panel(p_fig3f, "Fig3F_gain_enhancer_GO.png", 9.5, 7.2)

# Fig. 3G: apply the documented DEG thresholds and compare with gain-enhancer targets.
deg <- utils::read.csv(data_file("KMT2C_mut_DESeq2_all_results.csv"), check.names = FALSE)
enhancer_links <- utils::read.csv(
  data_file("KMT2C_target_genes_rewriting gain enhancer.csv"), check.names = FALSE
)
# Normalize UTF-8-BOM-prefixed first-column names found in the supplied CSV files.
names(deg)[1] <- "gene"
names(enhancer_links)[1] <- "entrez"
enhancer_links$SYMBOL <- toupper(trimws(as.character(enhancer_links$SYMBOL)))
deg$gene <- toupper(trimws(as.character(deg$gene)))
chip_genes <- sort(unique(enhancer_links$SYMBOL[!is.na(enhancer_links$SYMBOL) & enhancer_links$SYMBOL != ""]))
upregulated_genes <- sort(unique(deg$gene[
  !is.na(deg$padj) & deg$padj < 0.05 & deg$log2FoldChange > 1 &
    !is.na(deg$gene) & deg$gene != ""
]))
overlap_genes <- intersect(chip_genes, upregulated_genes)

venn_membership <- data.frame(
  gene = sort(union(chip_genes, upregulated_genes)),
  stringsAsFactors = FALSE
)
venn_membership$gain_enhancer_target <- venn_membership$gene %in% chip_genes
venn_membership$upregulated_DEG <- venn_membership$gene %in% upregulated_genes
venn_membership$overlap <- venn_membership$gain_enhancer_target & venn_membership$upregulated_DEG
write_source_table(venn_membership, "Fig3G_venn_membership.csv")
write_source_table(
  data.frame(
    set = c("Gain-enhancer targets", "Upregulated DEGs", "Intersection"),
    count = c(length(chip_genes), length(upregulated_genes), length(overlap_genes))
  ),
  "Fig3G_venn_counts.csv"
)

grDevices::png(panel_file("Fig3G_gain_enhancer_DEG_venn.png"), width = 3500, height = 3000, res = 600)
grid::grid.newpage()
grid::grid.draw(VennDiagram::draw.pairwise.venn(
  area1 = length(chip_genes), area2 = length(upregulated_genes),
  cross.area = length(overlap_genes),
  category = c("KMT2C-KO gain enhancers", "Upregulated DEGs"),
  fill = c("#2C7BB6", "#F0A202"), alpha = c(0.55, 0.55),
  scaled = FALSE, cat.cex = 1.2, cex = 1.5
))
grDevices::dev.off()

# Table 2: count detailed gained-enhancer associations for each overlapping human gene.
enhancer_associations <- utils::read.csv(
  data_file("KMT2C_gain_enhancer_gene_associations.csv"), check.names = FALSE
)
enhancer_associations$Human_Gene <- toupper(trimws(enhancer_associations$Human_Gene))
write_source_table(enhancer_associations, "Table2_enhancer_gene_associations.csv")
enhancer_counts <- as.data.frame(table(enhancer_associations$Human_Gene), stringsAsFactors = FALSE)
names(enhancer_counts) <- c("Human_Gene", "No_of_gained_enhancers")
enhancer_counts <- enhancer_counts[enhancer_counts$Human_Gene %in% overlap_genes, ]
deg_selected <- deg[deg$gene %in% overlap_genes, c("gene", "log2FoldChange", "padj")]
names(deg_selected) <- c("Human_Gene", "RNA_Log2FC", "RNA_FDR")
table2 <- merge(deg_selected, enhancer_counts, by = "Human_Gene")
table2 <- table2[order(-table2$No_of_gained_enhancers, table2$Human_Gene), ]
write_source_table(table2, "Table2_gain_enhancer_upregulated_genes.csv")
