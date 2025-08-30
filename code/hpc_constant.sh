#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -J 1-1000
#PBS -j oe
#PBS -N simulation

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Define Ne values after crash
Ne_current=(300 3000)

# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/{constant.py,decline.py,statistics.py} $TMPDIR

# Run the simulations for each Ne value
for Ne_current in "${Ne_current[@]}"; do
    python constant.py
    python decline.py $Ne_current
    python statistics.py
    
    # Copy the output files back to home directory with Ne value in filename
    mkdir -p $HOME/outputs/Constant_$Ne_current
    mv "$TMPDIR/constant_"* "$TMPDIR/decline_"* $HOME/outputs/Constant_$Ne_current/
done

# Move log files to logs directory
mkdir -p $HOME/logs
mv "$HOME/"*.o* $HOME/logs/

echo "msprime and statistics calculation has finished running"