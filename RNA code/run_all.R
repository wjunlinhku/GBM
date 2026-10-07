# Reproduce all target panels and assemble SourceData.xlsx from a clean R session.

# Locate this script even when Rscript is called from another working directory.
command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this workflow with: Rscript run_all.R")
}
package_dir <- normalizePath(
  dirname(sub("^--file=", "", file_argument)),
  winslash = "/", mustWork = TRUE
)
project_dir <- package_dir
options(nc.package_dir = package_dir, nc.project_dir = project_dir)

# Create deterministic output folders without changing the caller's working directory.
dir.create(file.path(package_dir, "output", "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(package_dir, "output", "source_tables"), recursive = TRUE, showWarnings = FALSE)

# Run each analysis module in panel order.
scripts <- c(
  "00_common.R",
  "01_fig2_and_figS3_distributions.R",
  "02_fig3AB_original_RNA_part1.R",
  "02_fig3_pathway_and_wgcna.R",
  "03_fig3_enhancer_analysis.R",
  "04_fig2_fig4_and_supplementary_statistics.R",
  "05_fig4_ligand_receptor.R",
  "06_build_source_data_workbook.R"
)
for (script in scripts) {
  message("Running ", script)
  sys.source(file.path(package_dir, "R", script), envir = globalenv())
}

# Record the exact R and package versions used for this run.
session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(package_dir, "sessionInfo.txt"), useBytes = TRUE)
message("Workflow completed. Outputs: ", file.path(package_dir, "output"))
