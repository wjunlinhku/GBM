# Reproduce Fig. 4F and Supplementary Fig. 5B from the supplied ligand-receptor table.

require_packages(c("ggplot2", "scales"))

# Load and filter the two pathways excluded by the original script.
lr_raw <- utils::read.csv(
  data_file("group Recurrent_NE vs Primary_NE.csv"), check.names = FALSE
)
lr_filtered <- lr_raw[
  !lr_raw$pw.name %in% c(
    "Anti-inflammatory response favouring Leishmania parasite infection",
    "Surfactant metabolism"
  ),
]
lr_filtered$LR_pair <- paste(lr_filtered$L, lr_filtered$R, sep = "-")
lr_filtered$mean_log2FC <- (lr_filtered$L.logFC + lr_filtered$R.logFC) / 2

# Aggregate repeated pathway/pair records exactly as in the plotting code.
lr_aggregated <- stats::aggregate(
  cbind(LR.score, mean_log2FC) ~ pw.name + LR_pair,
  data = lr_filtered, FUN = mean, na.rm = TRUE
)
names(lr_aggregated)[3:4] <- c("mean_LR_score", "mean_log2FC")

# Build a reusable dot plot with pathway and LR-pair ordering based on the displayed values.
make_lr_plot <- function(data, title, output_name, width, height) {
  pathway_order <- names(sort(tapply(data$mean_log2FC, data$pw.name, mean)))
  pair_order <- names(sort(tapply(data$mean_LR_score, data$LR_pair, mean)))
  data$pw.name <- factor(data$pw.name, levels = pathway_order)
  data$LR_pair <- factor(data$LR_pair, levels = pair_order)
  p <- ggplot2::ggplot(
    data, ggplot2::aes(x = LR_pair, y = pw.name)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = mean_LR_score, color = mean_log2FC), alpha = 0.9
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-2, 2), oob = scales::squish,
      name = "Mean log2FC"
    ) +
    ggplot2::scale_size(range = c(1.5, 6), name = "LR score") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    theme_publication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold", size = 8)
    )
  save_panel(p, output_name, width, height)
}

# Fig. 4F: retain the narrative-selected pathway/LR-pair combinations.
selected_targets <- data.frame(
  LR_pair = c(
    "B2M-KLRD1", "HDC-HRH2", "EFNB3-EPHB2", "NTN4-DCC", "SLIT2-ROBO1",
    "GAS6-TYRO3", "VIM-CD44", "LRPAP1-LRP8", "TGS1-RXRA", "COL2A1-ITGA2B"
  ),
  pw.name = c(
    "DAP12 interactions",
    "ADORA2B mediated anti-inflammatory cytokines production",
    "Ephrin signaling",
    "Netrin mediated repulsion signals",
    "Regulation of commissural axon pathfinding by SLIT and ROBO",
    "nervous system development",
    "extracellular matrix organization",
    "ECM proteoglycans", "ECM proteoglycans", "ECM proteoglycans"
  ),
  stringsAsFactors = FALSE
)
fig4f <- merge(lr_aggregated, selected_targets, by = c("pw.name", "LR_pair"))
write_source_table(fig4f, "Fig4F_selected_ligand_receptor_pairs.csv")
make_lr_plot(
  fig4f,
  "Recurrent NE activated pathway-LR pairs compared with Primary NE",
  "Fig4F_selected_ligand_receptor_pairs.png", 12.0, 6.0
)

# Supplementary Fig. 5B: retain all rows belonging to the top 40 LR pairs by maximum score.
pair_max <- stats::aggregate(mean_LR_score ~ LR_pair, lr_aggregated, max, na.rm = TRUE)
top_pairs <- head(pair_max$LR_pair[order(pair_max$mean_LR_score, decreasing = TRUE)], 40)
figs5b <- lr_aggregated[lr_aggregated$LR_pair %in% top_pairs, ]
write_source_table(figs5b, "FigS5B_top40_ligand_receptor_pairs.csv")
make_lr_plot(
  figs5b,
  "Recurrent NE activated pathway-LR pairs compared with Primary NE",
  "FigS5B_top40_ligand_receptor_pairs.png", 12.0, 18.0
)
