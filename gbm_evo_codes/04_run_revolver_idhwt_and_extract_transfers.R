#!/usr/bin/env Rscript

# Stage 4: run REVOLVER for the IDH-WT cohort and export information transfers.

library(ggplot2)

input_file <- Sys.getenv("REVOLVER_IDHWT_INPUT_FILE", "/path/to/project/revolver/by_idh_status/revolver_input_idhWT.tsv")
output_dir <- Sys.getenv("REVOLVER_IDHWT_OUTPUT_DIR", "/path/to/project/revolver/revolver_output_idhWT")

dirs <- c(
  output_dir,
  file.path(output_dir, "1_plot"),
  file.path(output_dir, "patient_full_trees"),
  file.path(output_dir, "patient_only_trees"),
  file.path(output_dir, "clusters")
)
lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

input_data <- readr::read_tsv(input_file, locale = locale(encoding = "UTF-8"))

input_data$is.driver <- as.logical(input_data$is.driver)
input_data$is.clonal <- as.logical(input_data$is.clonal)
input_data$cluster   <- as.character(input_data$cluster)
input_data$CCF       <- as.character(input_data$CCF)

input_data <- input_data %>%
  filter(
    !is.na(patientID), patientID != "",
    !is.na(variantID), variantID != "",
    !is.na(cluster), cluster != "",
    !is.na(is.driver),
    !is.na(is.clonal),
    !is.na(CCF), CCF != ""
  )

cat("Running REVOLVER for", length(unique(input_data$patientID)), "patients.\n")
rev_obj <- revolver_cohort(
  dataset = input_data,
  CCF_parser = revolver::CCF_parser,
  ONLY.DRIVER = FALSE,
  MIN.CLUSTER.SIZE = 3,
  annotation = "GBM REVOLVER - "
)

patients <- unique(input_data$patientID)
successful_patients <- c()

if(is.null(rev_obj$phylogenies)) rev_obj$phylogenies <- list()

for(p in patients) {
  tryCatch({
    temp_obj <- compute_clone_trees(rev_obj, patients = p, overwrite = TRUE)
    if(!is.null(temp_obj$phylogenies) && p %in% names(temp_obj$phylogenies)) {
      rev_obj$phylogenies[[p]] <- temp_obj$phylogenies[[p]]
      successful_patients <- c(successful_patients, p)
    }
  }, error = function(e) {
    warning("Clone-tree inference failed for ", p, ": ", conditionMessage(e))
  })
}

if(is.null(rev_obj$variantIDs_by_cluster) || nrow(rev_obj$variantIDs_by_cluster) == 0) {
  warning("REVOLVER did not return variant-cluster assignments; using driver variants from input.")
  rev_obj$variantIDs_by_cluster <- input_data[input_data$is.driver == TRUE, 
                                            c("patientID", "variantID", "cluster", "is.driver", "is.clonal")]
}

saveRDS(rev_obj, file.path(output_dir, "revolver_trees_processed.rds"))

# Driver Occurrence Top 50
driver_counts <- as.data.frame(rev_obj$variantIDs_by_cluster) %>%
  group_by(variantID) %>%
  summarise(occurrence = n()) %>%
  arrange(desc(occurrence))

top50_genes <- head(driver_counts$variantID, 50)
rev_obj_top50 <- rev_obj
rev_obj_top50$variantIDs_by_cluster <- rev_obj$variantIDs_by_cluster[rev_obj$variantIDs_by_cluster$variantID %in% top50_genes, ]

pdf(file.path(output_dir, "1_plot/drivers_occurrence_top50.pdf"), width=15, height=12)
print(plot_drivers_occurrence(rev_obj_top50))
dev.off()

patients_with_trees <- names(rev_obj$phylogenies)

variant_ids <- rev_obj$variantIDs_by_cluster$variantID
driver_counts <- table(variant_ids)
common_drivers <- names(driver_counts[driver_counts >= 2])

if(length(common_drivers) < 2) {
  warning("Fewer than two recurrent drivers were found; using the ten most frequent variants.")
  all_counts <- table(input_data$variantID)
  top_vars <- names(sort(all_counts, decreasing = TRUE))[1:10]
  rev_obj$variantIDs_by_cluster <- input_data[input_data$variantID %in% top_vars, 
                                            c("patientID", "variantID", "cluster", "is.driver", "is.clonal")]
  rev_obj$variantIDs_by_cluster$is.driver <- TRUE
}

tryCatch({
  rev_obj_fit <- revolver_fit(rev_obj, 
                             patients = patients_with_trees,
                             n.attempts = 1, 
                             auto.merge = TRUE,
                             force = TRUE)
  
  saveRDS(rev_obj_fit, file.path(output_dir, "revolver_fit_result.rds"))
  
  pdf(file.path(output_dir, "fit_information.pdf"), width = 10, height = 8)
  print(plot_fit(rev_obj_fit))
  dev.off()
}, error = function(e) {
  warning("REVOLVER fit failed: ", conditionMessage(e))
})

if(exists("rev_obj_fit")) {
  for(k in 2:4) {
    tryCatch({
      rev_obj_cluster <- revolver_cluster(rev_obj_fit, k = k)
      saveRDS(rev_obj_cluster, file.path(output_dir, "clusters", paste0("cluster_k", k, ".rds")))
      
      pdf(file.path(output_dir, "clusters", paste0("cluster_trajectories_k", k, ".pdf")), width = 12, height = 10)
      print(plot_trajectories(rev_obj_cluster))
      dev.off()
    }, error = function(e) {
      warning("REVOLVER clustering failed for k = ", k, ": ", conditionMessage(e))
    })
  }
}


# Information-transfer extraction stage.


output_file <- Sys.getenv("REVOLVER_ITRANSFER_FILE", "/path/to/project/revolver/revolver_output_idhWT/information_transfer_results.txt")

rev_obj_path <- Sys.getenv("REVOLVER_TREES_RDS", "/path/to/project/revolver/revolver_output_idhWT/revolver_trees_processed.rds")
rev_obj <- readRDS(rev_obj_path)

patients <- names(rev_obj$phylogenies)

file.create(output_file)

for (patient_id in patients) {
  tryCatch({
    if (!is.null(rev_obj$phylogenies[[patient_id]]) && length(rev_obj$phylogenies[[patient_id]]) > 0) {
      transfer_result <- ITransfer(rev_obj, patient_id, rank = 1, type = 'drivers')
      
      write(paste0("patientID: ", patient_id), file = output_file, append = TRUE)
      
      write(" from  to   ", file = output_file, append = TRUE)
      
      for (i in 1:nrow(transfer_result)) {
        line <- sprintf("%-2d %-5s %-5s", i, transfer_result$from[i], transfer_result$to[i])
        write(line, file = output_file, append = TRUE)
      }
      
      write("", file = output_file, append = TRUE)
    } else {
      warning("No phylogeny was available for patient: ", patient_id)
    }
  }, error = function(e) {
    warning("Information-transfer extraction failed for ", patient_id, ": ", conditionMessage(e))
  })
}

if (file.exists(output_file)) {
  lines <- readLines(output_file)
  patient_count <- sum(grepl("^patientID:", lines))
  cat("Information transfers written for", patient_count, "patients:", output_file, "\n")
} else {
  stop("Information-transfer output was not created: ", output_file)
}
