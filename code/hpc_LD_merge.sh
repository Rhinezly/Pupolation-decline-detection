#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -j oe
#PBS -N merge

module load anaconda3/personal

OUTPUT_DIR="$HOME/outputs/LD"

# List of scenarios to process
MODELS=("con" "dec300" "dec3000" "s3" "s6" "sd3_300" "sd3_3000" "sd6_300" "sd6_3000")

for MODEL in "${MODELS[@]}"; do
    MODEL_DIR="$OUTPUT_DIR/$MODEL"
    if [ -d "$MODEL_DIR" ]; then
        echo "Merging CSV files for model $MODEL ..."
        cd "$MODEL_DIR"
        Rscript $HOME/LD_merge.R "$MODEL"
        mv *_merged.csv $OUTPUT_DIR
    else
        echo "Directory $MODEL_DIR does not exist, skipping."
    fi
done

# Move log files to logs directory
mkdir -p $HOME/logs
mv "$HOME/"*.o* $HOME/logs/