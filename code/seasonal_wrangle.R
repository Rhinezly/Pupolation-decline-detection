library(tidyverse)
library(reticulate)
np <- import("numpy")

stat_specs <- tribble(
  ~orig_name,                     ~output_name,
  "Number of trees",              "Number_of_trees",
  "Number of mutations",          "Number_of_mutations",
  "Number of polymorphic sites",  "Number_of_polymorphic_sites",
  "Density of segregating sites", "Density_of_segregating_sites",
  "Nucleotide diversity (pi)",    "Nucleotide_diversity_pi",
  "Tajima's D",                   "Tajimas_D"
)

# Function to process NPZ files with Python
process_npz_files <- function(input_dir, output_dir) {
  npz_files <- list.files(input_dir, pattern = "\\.npz$", full.names = TRUE)
  if (length(npz_files) == 0) return()
  
  py_run_string("
import os
import numpy as np

def save_npz_data(files, model_type, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    stats = {k: [] for k in ['fst', 'dxy', 'gr']}
    
    for file in files:
        data = np.load(file)
        for stat in stats.keys():
            stats[stat].append(data[stat])
    
    for stat, arr in stats.items():
        np.save(f'{output_dir}/{model_type}_{stat}.npy', np.stack(arr))
")

  # Process files by model type
  list(
    seasonal = grep("seasonal_", npz_files, value = TRUE),
    decline = grep("decline_", npz_files, value = TRUE)
  ) %>% 
    walk2(names(.), function(files, model_type) {
      if (length(files) > 0) py$save_npz_data(files, model_type, output_dir)
    })
}


##### Main function to process a single directory #####

process_single_directory <- function(input_dir, output_root_dir) {
  drytime <- basename(input_dir) %>% str_extract("\\d+")
  drytype <- basename(dirname(input_dir))
  
  # Create output directory
  output_dir <- file.path(output_root_dir, drytype, paste0("DryTime_", drytime))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process CSV files
  list.files(input_dir, pattern = "\\.csv$", full.names = TRUE) %>%
    map_dfr(function(file) {
      seed_num <- str_extract(basename(file), "\\d+(?=\\.csv)") %>% as.integer()
      model_type <- ifelse(str_detect(file, "seasonal_"), "seasonal", "decline")
      
      read_csv(file, show_col_types = FALSE) %>%
        rename(Time = 1) %>%
        mutate(Seed = seed_num, Model = model_type)
    }) %>%
    pivot_longer(-c(Time, Seed, Model), names_to = "Statistic", values_to = "Value") %>%
    filter(Statistic %in% stat_specs$orig_name) %>%
    left_join(stat_specs, by = c("Statistic" = "orig_name")) %>%
    select(-Statistic) %>%
    pivot_wider(names_from = Time, values_from = Value) %>%
    group_by(output_name) %>%
    group_split() %>%
    walk(~{
      write_csv(select(.x, -output_name), 
               file.path(output_dir, paste0(first(.x$output_name), ".csv")))
    })
  
  # Process NPZ files
  process_npz_files(input_dir, output_dir)
}