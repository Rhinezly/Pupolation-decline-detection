#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -j oe
#PBS -N wrangle


module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Set working directories
INPUT_DIR="$HOME/outputs"
OUTPUT_DIR="$INPUT_DIR/statistics" 

# Process each directory
for Ne_DIR in "$INPUT_DIR"/Constant_*; do
  Rscript --vanilla -e "
    library(reticulate);
    use_python('$(which python)', required = TRUE);
    source('$HOME/constant_wrangle.R');
    process_single_directory('$Ne_DIR', '$OUTPUT_DIR')
  "
done

# Move log files to logs/ directory
mv "$PBS_O_WORKDIR"/*.o* "$HOME/logs/"

echo "Job completed. Results in $OUTPUT_DIR"