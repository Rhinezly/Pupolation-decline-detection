library(tidyverse)

# Read in CSV files
constant_files <- Sys.glob("constant_*.csv")
decline_files <- Sys.glob("decline_*.csv") 

# Define the mapping of original names to output names
stat_specs <- tribble(
  ~orig_name,                     ~output_name,
  "Number of trees",              "Number_of_trees",
  "Number of mutations",          "Number_of_mutations",
  "Density of segregating sites", "Density_of_segregating_sites",
  "Nucleotide diversity (pi)",    "Nucleotide_diversity_pi",
  "Tajima's D",                   "Tajimas_D"
)

# Function to process files of a specific type (constant or decline)
process_file_type <- function(files, type) {
  map_dfr(files, function(file) {
    # Extract numbers from filenames（e.g. "constant_1.csv" → 1）
    file_num <- str_extract(file, "(?<=_)\\d+") %>% as.integer()
    
    # Read the CSV file and process it
    read_csv(file, show_col_types = FALSE) %>%
      rename(Time = 1) %>%
      mutate(
        File = tools::file_path_sans_ext(basename(file)),
        FileNum = file_num,
        FileType = type  # annotate file types（constant/decline）
      )
  }) %>%
  arrange(FileNum)  # arrange by FileNum
}

# Process both file types and combine the results
combined_data <- bind_rows(
  process_file_type(constant_files, "constant"),
  process_file_type(decline_files, "decline")
) %>%
  # Transform data to long format
  pivot_longer(
    cols = -c(Time, File, FileNum, FileType),
    names_to = "Statistic",
    values_to = "Value"
  ) %>%
  filter(Statistic %in% stat_specs$orig_name) %>%
  # Join with stat_specs to get output names
  left_join(stat_specs, by = c("Statistic" = "orig_name")) %>%
  select(-Statistic, -FileNum, -FileType) %>%  # remove temporary columns
  # Transform to wide format
  pivot_wider(
    names_from = Time,
    values_from = Value
  ) %>%
  # Group by output_name and summarize
  group_by(output_name) %>%
  group_split()

# Save each group to a separate CSV file
walk(combined_data, ~{
  output_file <- paste0(first(.x$output_name), ".csv")
  write_csv(.x %>% select(-output_name), output_file)
  message("Generated: ", output_file)
})

message("All files processed successfully")