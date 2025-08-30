#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -J 1-1000
#PBS -j oe
#PBS -N sim_all_models

module load anaconda3/personal

# # List of scenarios to process
MODELS=(
  "con" "dec300" "dec3000" "s3" "s6" 
  "sd3_300" "sd3_3000" "sd6_300" "sd6_3000"
)

cp $HOME/{LD.R,aestivation3.so} $TMPDIR/

# Run all models
for MODEL in "${MODELS[@]}"; do
    echo "Running model $MODEL replicate ${PBS_ARRAY_INDEX}"

    Rscript LD.R ${PBS_ARRAY_INDEX} ${MODEL}

    OUTPUT_DIR=$HOME/outputs/LD/$MODEL
    mkdir -p $OUTPUT_DIR
    mv "${MODEL}_${PBS_ARRAY_INDEX}.csv" "$OUTPUT_DIR/"
done


mkdir -p $HOME/logs
mv "$HOME/"*.o* $HOME/logs/