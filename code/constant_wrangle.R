library(tidyverse)
library(reticulate)  # Import reticulate for Python integration
np <- import("numpy")

cat("Python configuration:\n")
print(reticulate::py_config())
cat("NumPy available:", reticulate::py_module_available("numpy"), "\n")

# Define statistic name mappings
stat_specs <- tribble(
  ~orig_name,                     ~output_name,
  "Number of trees",              "Number_of_trees",
  "Number of mutations",          "Number_of_mutations",
  "Density of segregating sites", "Density_of_segregating_sites",
  "Nucleotide diversity (pi)",    "Nucleotide_diversity_pi",
  "Tajima's D",                   "Tajimas_D"
)

# Function to process NPZ files with Python
process_npz_files <- function(input_dir, output_root_dir) {
  npz_files <- list.files(input_dir, pattern = "\\.npz$", full.names = TRUE)
  cat("Found", length(npz_files), ".npz files in", input_dir, "\n")
  
  if (length(npz_files) == 0) return()
  
  py_run_string("
import numpy as np
import glob

def process_npz_group(files, model_type, output_dir):
    fst_list = []
    dxy_list = []
    gr_list = []
    
    for file in files:
        data = np.load(file)
        fst_list.append(data['fst'])
        dxy_list.append(data['dxy'])
        gr_list.append(data['gr'])
    
    np.save(f'{output_dir}/{model_type}_fst.npy', np.stack(fst_list))
    np.save(f'{output_dir}/{model_type}_dxy.npy', np.stack(dxy_list))
    np.save(f'{output_dir}/{model_type}_gr.npy', np.stack(gr_list))
")
  
  # Separate files by model type
  constant_files <- npz_files[grep("constant_", npz_files)]
  decline_files <- npz_files[grep("decline_", npz_files)]
  
  drytime <- basename(input_dir) %>% str_extract("\\d+")
  output_dir <- file.path(output_root_dir, paste0("Constant", drytime))
  
  if (length(constant_files) > 0) {
    py$process_npz_group(constant_files, "constant", output_dir)
    cat("Processed", length(constant_files), "constant files\n")
  }
  
  if (length(decline_files) > 0) {
    py$process_npz_group(decline_files, "decline", output_dir)
    cat("Processed", length(decline_files), "decline files\n")
  }
}



##### Main function to process a single directory #####

process_single_directory <- function(input_dir, output_root_dir) {
  # Process CSV files
  list.files(input_dir, pattern = "\\.csv$", full.names = TRUE) %>%
    map_dfr(function(file) {
      seed_num <- str_extract(basename(file), "\\d+(?=\\.csv)") %>% as.integer()
      model_type <- ifelse(str_detect(file, "constant_"), "constant", "decline")
      
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
      stat_name <- first(.x$output_name)
      output_dir <- file.path(output_root_dir, paste0("Constant", basename(input_dir) %>% str_extract("\\d+")))
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      write_csv(select(.x, -output_name), file.path(output_dir, paste0(stat_name, ".csv")))
    })
  
  # Process NPZ files
  process_npz_files(input_dir, output_root_dir)
}