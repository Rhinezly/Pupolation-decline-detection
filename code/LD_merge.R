args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1){
  stop("Usage: Rscript merge_model_csv.R <model_name>")
}

model_name <- args[1]

# Find all CSV files in the current directory
files <- list.files(pattern = "\\.csv$", full.names = TRUE)

if(length(files) == 0){
  stop("No CSV files found for model ", model_name)
}

# Read and merge CSV files
df_list <- lapply(files, read.csv)
df_all <- do.call(rbind, df_list)

# Write merged data to a new CSV file
output_file <- paste0(model_name, "_merged.csv")
write.csv(df_all, output_file, row.names = FALSE)
cat("Merged file saved to:", output_file, "\n")