#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -J 1-1000
#PBS -j oe
#PBS -N simulation

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/{constant.py,decline.py,statistics.py} $TMPDIR

# Run the simulations
python constant.py
python decline.py
python statistics.py

# Copy the output files back to home directory
mv "$TMPDIR/constant_"* "$TMPDIR/decline_"* $HOME/outputs/

# Move log files to $HOME/logs directory
mv "$HOME/"*.o* $HOME/logs/

echo "msprime and statistics calculation has finished running"