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

# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/wrangle.R $TMPDIR

# Copy input data to temporary directory
time (
  echo "Starting parallel copy..."
  cp "$HOME/outputs/constant_"*.csv $TMPDIR &
  cp "$HOME/outputs/decline_"*.csv $TMPDIR &
  wait
)
echo "Copy completed in $SECONDS seconds"

# Run R script in temporary directory
R --vanilla < wrangle.R

# Move results and clean up
for stat in "Number_of_trees" "Number_of_mutations" \
             "Number_of_polymorphic_sites" "Density_of_segregating_sites" \
             "Nucleotide_diversity_pi" "Tajimas_D"; do
  mv "$TMPDIR/$stat.csv" "$OUTPUT_DIR/" 2>/dev/null
done

# Move log files to $HOME/logs directory
mv "$HOME/"*.o* $HOME/logs/

echo "Job completed successfully. Results saved to $OUTPUT_DIR"