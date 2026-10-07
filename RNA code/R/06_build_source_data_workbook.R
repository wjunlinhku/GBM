# Assemble all plotted values and statistical results into SourceData.xlsx.

require_packages(c("openxlsx", "zip"))

# Remove dangling drawing relationships produced by openxlsx 4.2.5.2 on this platform.
# These relationships point to files that do not exist and can trigger repair prompts in
# non-R spreadsheet readers even though openxlsx itself can read the workbook.
repair_dangling_drawing_relationships <- function(path) {
  extract_dir <- tempfile("xlsx_repair_")
  dir.create(extract_dir)
  utils::unzip(path, exdir = extract_dir)
  relationship_dir <- file.path(extract_dir, "xl", "worksheets", "_rels")
  relationship_files <- list.files(
    relationship_dir, pattern = "[.]rels$", full.names = TRUE
  )
  for (relationship_file in relationship_files) {
    xml <- paste(readLines(relationship_file, warn = FALSE), collapse = "")
    xml <- gsub(
      "<Relationship[^>]+/relationships/(drawing|vmlDrawing)\"[^>]*/>",
      "", xml, perl = TRUE
    )
    writeLines(xml, relationship_file, useBytes = TRUE)
  }
  repaired_path <- tempfile(fileext = ".xlsx")
  zip::zipr(
    repaired_path,
    files = list.files(extract_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    root = extract_dir, include_directories = FALSE, mode = "mirror"
  )
  file.copy(repaired_path, path, overwrite = TRUE)
  unlink(extract_dir, recursive = TRUE)
  unlink(repaired_path)
  invisible(path)
}

# Define one worksheet per main or supplementary panel.
sheet_blocks <- list(
  Fig2A = list(c("Plotted counts and proportions", "Fig2A_clinical_region_distribution.csv")),
  Fig2B = list(c("Plotted counts and proportions", "Fig2B_RNA_subtype_distribution.csv")),
  Fig2C = list(
    c("Plotted row z-scores", "Fig2C_expression_zscore.csv"),
    c("Underlying voom expression", "Fig2C_expression_voom.csv")
  ),
  Fig2D = list(c("All plotted genes", "Fig2D_EGFR_vs_nonEGFR_DEG.csv")),
  Fig2E = list(c("All plotted KEGG terms", "Fig2E_EGFR_upregulated_KEGG.csv")),
  Fig2F = list(c("Plotted group means", "Fig2F_mean_cell_composition.csv")),
  Fig3A = list(
    c("Displayed coefficient matrix (Pearson r)", "Fig3A_correlation_matrix.csv"),
    c("Displayed significance symbols from original code", "Fig3A_significance_symbols.csv"),
    c("Significance P-value matrix (Spearman test)", "Fig3A_pvalue_matrix.csv"),
    c("Complete-pair sample sizes", "Fig3A_complete_n_matrix.csv"),
    c("ssGSEA scores", "Fig3A_ssGSEA_scores.csv")
  ),
  Fig3B = list(
    c("Displayed coefficient matrix (Pearson r)", "Fig3B_correlation_matrix.csv"),
    c("Displayed significance symbols from original code", "Fig3B_significance_symbols.csv"),
    c("Significance P-value matrix (Spearman test)", "Fig3B_pvalue_matrix.csv"),
    c("Complete-pair sample sizes", "Fig3B_complete_n_matrix.csv"),
    c("ssGSEA scores", "Fig3B_ssGSEA_scores.csv")
  ),
  Fig3C = list(
    c("Plotted hub-gene z-scores in heatmap order", "Fig3C_hub_gene_expression_zscore.csv"),
    c("Underlying hub-gene expression", "Fig3C_hub_gene_expression_raw.csv"),
    c("Sample annotation", "Fig3C_sample_annotation.csv"),
    c("Purple-module membership and gene significance", "Fig3C_purple_module_MM_GS.csv")
  ),
  Fig3D = list(c("Purple-module GO enrichment", "Fig3D_purple_module_GO_BP.csv")),
  Fig3E = list(
    c("Sample-level plotted values", "Fig3E_primary_NE_myeloid_values.csv"),
    c("Pairwise Welch tests", "Fig3E_pairwise_Welch_tests.csv")
  ),
  Fig3F = list(c("Selected gain-enhancer GO terms", "Fig3F_gain_enhancer_GO_selected.csv")),
  Fig3G = list(
    c("Venn counts", "Fig3G_venn_counts.csv"),
    c("Gene-level set membership", "Fig3G_venn_membership.csv")
  ),
  Fig4A = list(c("Plotted subtype counts and proportions", "Fig4A_RNA_subtype_by_clinical_region.csv")),
  Fig4B = list(c("Plotted group means", "Fig4B_mean_cell_composition.csv")),
  Fig4C = list(c("All plotted genes", "Fig4C_Recurrent_NE_vs_Primary_NE_DEG.csv")),
  Fig4D = list(
    c("Sample-level Myeloid and Malignant values", "Fig4D_NE_cell_fraction_values.csv"),
    c("Myeloid Welch test", "Fig4D_Myeloid_Welch_test.csv"),
    c("Malignant Welch test", "Fig4D_Malignant_Welch_test.csv")
  ),
  Fig4E = list(c("Selected signed KEGG values", "Fig4E_selected_KEGG_mirror.csv")),
  Fig4F = list(c("Selected aggregated LR pairs", "Fig4F_selected_ligand_receptor_pairs.csv")),
  FigS2C = list(
    c("Sample-level Oligodendrocyte values", "FigS2C_Oligodendrocyte_values.csv"),
    c("Welch test", "FigS2C_Oligodendrocyte_Welch_test.csv")
  ),
  FigS3A = list(c("BP and MF enrichment values", "FigS3A_EGFR_downregulated_GO.csv")),
  FigS3B = list(
    c("Plotted enrichment matrix", "FigS3B_region_enrichment_matrix.csv"),
    c("Counts, expected counts, and enrichment", "FigS3B_region_enrichment_long.csv")
  ),
  FigS3C = list(
    c("Plotted enrichment matrix", "FigS3C_RNA_subtype_enrichment_matrix.csv"),
    c("Counts, expected counts, and enrichment", "FigS3C_RNA_subtype_enrichment_long.csv")
  ),
  FigS3D = list(c("Plotted counts and proportions", "FigS3D_primary_recurrent_distribution.csv")),
  FigS3E = list(c("Plotted counts and proportions", "FigS3E_CE_NE_distribution.csv")),
  FigS5A = list(c("All signed KEGG values", "FigS5A_all_KEGG_mirror.csv")),
  FigS5B = list(c("Top-40 LR-pair aggregated values", "FigS5B_top40_ligand_receptor_pairs.csv")),
  Table2 = list(
    c("Gain-enhancer targets that are upregulated DEGs", "Table2_gain_enhancer_upregulated_genes.csv"),
    c("Underlying enhancer-gene associations", "Table2_enhancer_gene_associations.csv")
  )
)

# Initialize workbook styles.
workbook <- openxlsx::createWorkbook(creator = "Wang et al.")
title_style <- openxlsx::createStyle(
  fontSize = 12, textDecoration = "bold", fgFill = "#D9EAF7"
)
header_style <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = "#E7E6E6", border = "Bottom"
)
note_style <- openxlsx::createStyle(
  fontColour = "#666666", textDecoration = "italic", wrapText = TRUE
)

# Add workbook-level documentation and the key statistical caveat.
openxlsx::addWorksheet(workbook, "README")
readme_rows <- data.frame(
  Field = c(
    "Purpose", "Input location", "Code location", "Contrasts",
    "DEG threshold", "Enrichment correction", "Cell-composition tests",
    "Fig3A/B coefficient and significance calculation",
    "Fig3C scaling", "Missing values", "Reproduction command"
  ),
  Description = c(
    "Numeric source data underlying Fig. 2-4, Fig. S2C, Fig. S3, Fig. S5, and Table 2.",
    "All inputs are read from the project-relative source data directory.",
    "All scripts are in R/.",
    "Positive log2 fold change is the first named group relative to the second: EGFR versus non-EGFR; Recurrent_NE versus Primary_NE; KMT2C-mutant versus comparator.",
    "Differentially expressed genes used for the enhancer overlap require log2FC > 1 and BH-adjusted P < 0.05.",
    "Enrichment tables use their supplied BH-adjusted P values; plotted limits are 0-0.05.",
    "Two-sided Welch t-tests, matching the original plotting scripts; no multiplicity adjustment was applied to displayed pairwise tests.",
    "Fig3A/B are generated by the original RNA_part1.R logic: heatmap colors show default Pearson correlation coefficients, while displayed symbols come from two-sided asymptotic Spearman tests (exact = FALSE). The calculated P-value and significance-symbol matrices are supplied separately.",
    "Each gene was z-scored across the 44 WGCNA samples; samples were ordered by KMT2C mutation status and rows were hierarchically clustered.",
    "Blank cells represent missing values, not zero.",
    "From the package directory run: Rscript run_all.R"
  ),
  stringsAsFactors = FALSE
)
openxlsx::writeData(workbook, "README", readme_rows, withFilter = FALSE)
openxlsx::addStyle(workbook, "README", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
openxlsx::setColWidths(workbook, "README", cols = 1, widths = 32)
openxlsx::setColWidths(workbook, "README", cols = 2, widths = 110)
openxlsx::addStyle(workbook, "README", note_style, rows = 2:(nrow(readme_rows) + 1), cols = 2,
                   gridExpand = TRUE)
openxlsx::freezePane(workbook, "README", firstRow = TRUE)

# Add all panel tables, separating multiple blocks with a title and source filename.
for (sheet_name in names(sheet_blocks)) {
  openxlsx::addWorksheet(workbook, sheet_name)
  current_row <- 1L
  maximum_columns <- 1L
  for (block in sheet_blocks[[sheet_name]]) {
    title <- block[[1]]
    filename <- block[[2]]
    path <- table_file(filename)
    if (!file.exists(path)) stop("Missing generated source table: ", path)
    data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    maximum_columns <- max(maximum_columns, ncol(data))

    openxlsx::writeData(workbook, sheet_name, title, startRow = current_row, startCol = 1)
    openxlsx::addStyle(workbook, sheet_name, title_style, rows = current_row, cols = 1,
                       gridExpand = TRUE)
    current_row <- current_row + 1L
    openxlsx::writeData(workbook, sheet_name, paste0("Generated from: ", filename),
                        startRow = current_row, startCol = 1)
    openxlsx::addStyle(workbook, sheet_name, note_style, rows = current_row, cols = 1)
    current_row <- current_row + 1L
    openxlsx::writeData(workbook, sheet_name, data, startRow = current_row, startCol = 1,
                        withFilter = FALSE, keepNA = FALSE)
    openxlsx::addStyle(workbook, sheet_name, header_style, rows = current_row,
                       cols = seq_len(ncol(data)), gridExpand = TRUE)
    current_row <- current_row + nrow(data) + 3L
  }
  openxlsx::freezePane(workbook, sheet_name, firstActiveRow = 4)
  openxlsx::setColWidths(workbook, sheet_name, cols = seq_len(maximum_columns), widths = "auto")
}

# Save both inside the reproducibility package and at the project root for submission.
workbook_path <- output_file("SourceData.xlsx")
openxlsx::saveWorkbook(workbook, workbook_path, overwrite = TRUE)
repair_dangling_drawing_relationships(workbook_path)
file.copy(workbook_path, file.path(get_project_dir(), "SourceData.xlsx"), overwrite = TRUE)
