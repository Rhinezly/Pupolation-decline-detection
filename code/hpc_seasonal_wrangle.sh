#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -j oe
#PBS -N wrangle

module load anaconda3/personal
source activate msprime_env

BASE_DIR="$HOME/outputs"
OUTPUT_DIR="$BASE_DIR/statistics"

mkdir -p "$OUTPUT_DIR"/{Seasonal300,Seasonal3000}


# Process directories Seasonal3000 and Seasonal300
for crashtype in Seasonal300 Seasonal3000; do
  for drytime_dir in "$BASE_DIR/$crashtype"/DryTime_*; do
    echo "Processing: $drytime_dir"
    Rscript --vanilla -e "
      library(reticulate)
      use_python('$(which python)', required=TRUE)
      source('$HOME/seasonal_wrangle.R')
      process_single_directory('$drytime_dir', '$OUTPUT_DIR')
    "
  done
done


# Move log files to logs/ directory
mv "$PBS_O_WORKDIR"/*.o* "$HOME/logs/"

echo "Job completed. Results in $OUTPUT_DIR"