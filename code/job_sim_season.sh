#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -J 1-1000
#PBS -j oe
#PBS -N simulation_seasonal

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Define Ne values to test
time_dry=(6 3)

# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/{seasonal.py, statistics.py} $TMPDIR

# Run the simulations for each Ne value
for time_dry in "${time_dry[@]}"; do
    # Run with current Ne value
    python seasonal.py $time_dry
    python statistics.py
    
    # Copy the output files back to home directory with Ne in filename
    mkdir -p $HOME/outputs/DryTime_$time_dry
    mv "$TMPDIR/seasonal_"* $HOME/outputs/DryTime_$time_dry/
    mv "$TMPDIR/"*.csv "$TMPDIR/"*.npz $HOME/outputs/DryTime_$time_dry/ 2>/dev/null || true
done

# Move log files to $HOME/logs directory
mv "$HOME/"*.o* $HOME/logs/

echo "msprime and statistics calculation has finished running"