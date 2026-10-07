#!/usr/bin/env Rscript

# Stage 3: convert annotated PyClone-VI results to REVOLVER input.

# Input and output configuration.

library(tidyverse)
library(data.table)

pyclone_annotated_dir <- Sys.getenv("PYCLONE_ANNOTATED_DIR", "/path/to/project/pyclone_output/pyclone_results_with_gene_annotation")
output_file <- Sys.getenv("REVOLVER_INPUT_FILE", "/path/to/project/revolver/revolver_input.tsv")
driver_gene_file <- Sys.getenv("DRIVER_GENE_FILE", "/path/to/driver_gene_list.txt")

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

HANDLE_DUPLICATES <- TRUE

drivers <- read_tsv(driver_gene_file, col_names = FALSE) %>%
  pull(X1) %>%
  unique()

annotated_files <- list.files(pyclone_annotated_dir, 
                              pattern = "_results_with_gene_annotation\\.tsv$", 
                              full.names = TRUE)

if (length(annotated_files) == 0) {
  stop(" *_results_with_gene_annotation.tsv ")
}

revolver_data <- tibble()

for (file in annotated_files) {
  patient_id <- basename(file) %>% 
    str_replace("_results_with_gene_annotation.tsv", "")
  
  pyclone_results <- read_tsv(file, col_types = cols(.default = col_character()))
  
  pyclone_results <- pyclone_results %>%
    mutate(
      cellular_prevalence = as.numeric(cellular_prevalence),
      cluster_id = as.character(as.numeric(cluster_id))
    )
  
  clonal_info <- pyclone_results %>%
    group_by(cluster_id) %>%
    summarise(mean_cp = mean(cellular_prevalence, na.rm = TRUE)) %>%
    arrange(desc(mean_cp))
  
  clonal_cluster <- clonal_info$cluster_id[1]
  patient_data <- pyclone_results %>%
    distinct(mutation_id, cluster_id, .keep_all = TRUE) %>%
    mutate(
      cluster = as.character(as.integer(cluster_id) + 1),
      Misc = paste0(patient_id, ":", mutation_id),
      patientID = patient_id,
      variantID = ifelse(is.na(gene) | gene == "Unknown", mutation_id, gene),
      is.driver = !is.na(gene) & gene %in% drivers,
      is.clonal = (cluster_id == clonal_cluster),
      CCF = map_chr(mutation_id, function(mut) {
        samples <- pyclone_results %>% filter(mutation_id == mut)
        paste(paste0(samples$sample_id, ":", round(as.numeric(samples$cellular_prevalence), 4)),
              collapse = ";")
      }),
      mean_ccf = map_dbl(mutation_id, function(mut) {
        samples <- pyclone_results %>% filter(mutation_id == mut)
        mean(as.numeric(samples$cellular_prevalence), na.rm = TRUE)
      })
    )
  
  if (HANDLE_DUPLICATES) {
    duplicate_check <- patient_data %>%
      group_by(patientID, variantID) %>%
      summarise(count = n(), .groups = "drop") %>%
      filter(count > 1)
    
    if (nrow(duplicate_check) > 0) {
      patient_data <- patient_data %>%
        group_by(patientID, variantID) %>%
        slice_max(order_by = mean_ccf, n = 1, with_ties = FALSE) %>%
        ungroup()
      
    }
  }
  
  patient_output <- patient_data %>%
    select(Misc, patientID, variantID, cluster, is.driver, is.clonal, CCF,
           gene, func_type = func, exonic_func_type = exonic_func, aa_change_info = aa_change, mean_ccf)
  
  revolver_data <- bind_rows(revolver_data, patient_output)
}

if (HANDLE_DUPLICATES) {
  global_duplicates <- revolver_data %>%
    group_by(patientID, variantID) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)
  
  if (nrow(global_duplicates) > 0) {
    revolver_data <- revolver_data %>%
      group_by(patientID, variantID) %>%
      slice_max(order_by = mean_ccf, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
}

revolver_output <- revolver_data %>%
  select(Misc, patientID, variantID, cluster, is.driver, is.clonal, CCF)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write_tsv(revolver_output, output_file)

detailed_output_file <- file.path(output_dir, "revolver_input_with_gene_annotations.tsv")
write_tsv(revolver_data, detailed_output_file)

final_check <- revolver_output %>%
  group_by(patientID, variantID) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1)

if (nrow(final_check) == 0) {
  cat("REVOLVER input written:", output_file, "\n")
} else {
  warning("Duplicate patientID-variantID pairs remain in REVOLVER input: ", nrow(final_check))
}


patient_info_file <- Sys.getenv("PATIENT_INFO_FILE", "/path/to/patient_info.csv")

if (file.exists(patient_info_file)) {
  patient_info <- read_csv(patient_info_file, col_types = cols(.default = col_character())) %>%
    select(patientID, IDH_MUTation_status) %>%
    distinct(patientID, .keep_all = TRUE) %>%
    group_by(patientID) %>%
    summarise(IDH_MUTation_status = first(IDH_MUTation_status), 
              .groups = "drop") %>%
    mutate(
      IDH_MUTation_status = case_when(
        is.na(IDH_MUTation_status) | IDH_MUTation_status == "" ~ "NA",
        toupper(IDH_MUTation_status) == "WT" ~ "WT",
        toupper(IDH_MUTation_status) == "MUT" ~ "MUT",
        TRUE ~ "NA"  # NA
      )
    )
  
  revolver_input <- read_tsv(output_file, col_types = cols(.default = col_character()))
  
  revolver_with_idh <- revolver_input %>%
    left_join(patient_info, by = "patientID") %>%
    mutate(
      IDH_MUTation_status = ifelse(is.na(IDH_MUTation_status), "NA", IDH_MUTation_status)
    )
  
  missing_patients <- revolver_input %>%
    anti_join(patient_info, by = "patientID") %>%
    pull(patientID) %>%
    unique()
  
  if (length(missing_patients) > 0) {
    warning("Patients missing IDH status: ", paste(missing_patients, collapse = ", "))
  }
  
  output_dir_idh <- file.path(output_dir, "by_idh_status")
  if (!dir.exists(output_dir_idh)) {
    dir.create(output_dir_idh, recursive = TRUE)
  }
  
  wt_output <- revolver_with_idh %>%
    filter(IDH_MUTation_status == "WT") %>%
    select(-IDH_MUTation_status)
  
  wt_file <- file.path(output_dir_idh, "revolver_input_idhWT.tsv")
  write_tsv(wt_output, wt_file)
  
  mut_output <- revolver_with_idh %>%
    filter(IDH_MUTation_status == "MUT") %>%
    select(-IDH_MUTation_status)
  
  mut_file <- file.path(output_dir_idh, "revolver_input_idhMUT.tsv")
  write_tsv(mut_output, mut_file)
  
  na_output <- revolver_with_idh %>%
    filter(IDH_MUTation_status == "NA") %>%
    select(-IDH_MUTation_status)
  
  na_file <- file.path(output_dir_idh, "revolver_input_idhNA.tsv")
  write_tsv(na_output, na_file)
  
  patient_classification <- revolver_with_idh %>%
    distinct(patientID, IDH_MUTation_status) %>%
    group_by(patientID) %>%
    summarise(n_status = n_distinct(IDH_MUTation_status), .groups = "drop")
  
  if (max(patient_classification$n_status) > 1) {
    multi_status <- patient_classification %>% filter(n_status > 1)
    warning("Patients with multiple IDH classifications: ", paste(multi_status$patientID, collapse = ", "))
  }
  cat("IDH-stratified REVOLVER inputs written to:", output_dir_idh, "\n")
} else {
  warning("Patient information file not found; IDH-stratified outputs were not created: ", patient_info_file)
}
