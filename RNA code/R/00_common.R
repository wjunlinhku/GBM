# Shared utilities for all figure-reproduction scripts.

# Resolve the package and project roots from run_all.R rather than from a user-specific path.
get_package_dir <- function() {
  value <- getOption("nc.package_dir")
  if (is.null(value)) {
    stop("Option 'nc.package_dir' is not set. Run the workflow with run_all.R.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

get_project_dir <- function() {
  value <- getOption("nc.project_dir")
  if (is.null(value)) {
    stop("Option 'nc.project_dir' is not set. Run the workflow with run_all.R.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

data_file <- function(name) {
  path <- file.path(get_project_dir(), "source data", name)
  if (!file.exists(path)) stop("Missing input file: ", path)
  path
}

output_file <- function(...) {
  file.path(get_package_dir(), "output", ...)
}

panel_file <- function(name) output_file("panels", name)
table_file <- function(name) output_file("source_tables", name)

# Stop early when a required package is unavailable.
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Install the following R packages before running the workflow: ",
         paste(missing, collapse = ", "))
  }
}

# Convert cluster labels to the capitalization used in the figures.
format_group <- function(x) sub("^group", "Group", as.character(x))

# Convert a ratio such as "12/150" into a numeric value without evaluating code.
ratio_to_numeric <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[[1]]) / as.numeric(z[[2]])
  }, numeric(1))
}

# Use a consistent publication-style ggplot theme.
theme_publication <- function(base_size = 11) {
  # Use the portable sans alias (Arial on standard Windows installations).
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

# Save a figure panel with deterministic dimensions and high resolution.
save_panel <- function(plot, name, width, height, dpi = 600) {
  ggplot2::ggsave(
    filename = panel_file(name), plot = plot, width = width, height = height,
    units = "in", dpi = dpi, bg = "white"
  )
}

# Save a rectangular source-data table while preserving numeric cells.
write_source_table <- function(data, name, row_names = FALSE) {
  utils::write.csv(data, table_file(name), row.names = row_names, na = "")
}

# Save a matrix with its row identifiers in the first column.
write_source_matrix <- function(x, name, row_label = "row") {
  out <- data.frame(row_id = rownames(x), x, check.names = FALSE)
  names(out)[1] <- row_label
  write_source_table(out, name)
}

# Calculate expected counts and log2 observed-to-expected enrichment.
calculate_enrichment <- function(count_data, row_variable, column_variable) {
  row_totals <- stats::aggregate(
    count_data$count,
    list(count_data[[row_variable]]), sum
  )
  names(row_totals) <- c(row_variable, "row_total")
  column_totals <- stats::aggregate(
    count_data$count,
    list(count_data[[column_variable]]), sum
  )
  names(column_totals) <- c(column_variable, "column_total")
  result <- merge(count_data, row_totals, by = row_variable, all.x = TRUE)
  result <- merge(result, column_totals, by = column_variable, all.x = TRUE)
  result$expected <- result$row_total * result$column_total / sum(count_data$count)
  result$log2_observed_expected <- log2((result$count + 1) / (result$expected + 1))
  result
}

# Convert long enrichment data into a matrix with an explicit row/column order.
enrichment_matrix <- function(data, row_variable, column_variable,
                              row_levels = NULL, column_levels = NULL) {
  if (is.null(row_levels)) row_levels <- unique(data[[row_variable]])
  if (is.null(column_levels)) column_levels <- unique(data[[column_variable]])
  ans <- matrix(
    0, nrow = length(row_levels), ncol = length(column_levels),
    dimnames = list(row_levels, column_levels)
  )
  for (i in seq_len(nrow(data))) {
    ans[as.character(data[[row_variable]][i]),
        as.character(data[[column_variable]][i])] <- data$log2_observed_expected[i]
  }
  ans
}

# Return all pairwise Welch t-tests used by the legacy plotting scripts.
pairwise_welch_tests <- function(data, group, value) {
  groups <- unique(as.character(data[[group]]))
  combinations <- utils::combn(groups, 2, simplify = FALSE)
  rows <- lapply(combinations, function(z) {
    x <- data[data[[group]] == z[[1]], value]
    y <- data[data[[group]] == z[[2]], value]
    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) return(NULL)
    test <- stats::t.test(x, y)
    data.frame(
      group_1 = z[[1]], group_2 = z[[2]],
      n_1 = sum(!is.na(x)), n_2 = sum(!is.na(y)),
      mean_1 = mean(x, na.rm = TRUE), mean_2 = mean(y, na.rm = TRUE),
      p_value = test$p.value, test = "two-sided Welch t-test",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Build a standard box-and-jitter panel for cell-composition comparisons.
cell_boxplot <- function(data, group, value, title, palette = NULL) {
  p <- ggplot2::ggplot(
    data, ggplot2::aes(x = .data[[group]], y = .data[[value]], fill = .data[[group]])
  ) +
    ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, color = "black") +
    ggplot2::geom_jitter(width = 0.14, size = 1.4, alpha = 0.75) +
    ggplot2::labs(x = NULL, y = "Cell proportion", title = title) +
    theme_publication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_manual(values = palette)
  p
}
