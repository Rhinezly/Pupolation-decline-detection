#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -J 1-1000
#PBS -j oe
#PBS -N simulation_seasonal

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Define dry/wet season time
half_period=(6 3)
# Define initial population sizes
initial_sizes=(300 3000)


# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/{seasonal.py,decline_seasonal.py,statistics.py} $TMPDIR


# Run the simulations for each half period value and each initial size
for half_period in "${half_period[@]}"; do
    for initial_size in "${initial_sizes[@]}"; do
        python seasonal.py $half_period
        python decline_seasonal.py $half_period $initial_size
        python statistics.py
        
        # Copy the output files back to outputs directory with parameters in filename
        mkdir -p $HOME/outputs/Seasonal${initial_size}/DryTime_$half_period
        mv "$TMPDIR/seasonal_"* "$TMPDIR/decline_season_"* $HOME/outputs/Seasonal${initial_size}/DryTime_$half_period/
    done
done



# Move log files to $HOME/logs directory
mv "$HOME/"*.o* $HOME/logs/

echo "msprime and statistics calculation has finished running"