import msprime
import os
import sys

# If the shell script set the half_period values，use them; or else use the default
if len(sys.argv) > 1:
    half_period = int(sys.argv[1])
else:
    half_period = 6

seed_number = int(os.getenv("PBS_ARRAY_INDEX"))


### Parameters ###
Ne_wet = 30000  # Effective population size at wet season
Ne_dry = 300  # Population size at dry season
max_time = 1e6  # Maximum generations of simulation (for the while loop)
seqlength = 1e6  
mu = 3.75e-8  
recombination_rate = 10 * mu  
timepoints = list(range(0, 73, 3)) 


### Define demographic model ###
demography = msprime.Demography()
demography.add_population(name="pop", initial_size=Ne_wet)  # The most recent population size

# Add periodic population size changes
sim_time = half_period
sim_Ne = Ne_dry

while sim_time <= max_time:
    demography.add_population_parameters_change(
        time=sim_time,
        initial_size=sim_Ne,
        population="pop",
    )
    # Switch Ne
    sim_Ne = Ne_wet if sim_Ne == Ne_dry else Ne_dry
    sim_time += half_period

samples = [msprime.SampleSet(50, population="pop",time=t) for t in timepoints]


### Simulate ancestry ###
ts = msprime.sim_ancestry(
    samples=samples,
    demography=demography,
    sequence_length=seqlength, 
    recombination_rate=recombination_rate,
    model=[
        msprime.DiscreteTimeWrightFisher(duration=200),
        msprime.StandardCoalescent()
        ],
    random_seed=seed_number
)

# Add mutations
mts = msprime.sim_mutations(
    ts, 
    rate=mu, 
    random_seed=seed_number
)


### Save tree sequence ###
mts.dump(f"seasonal_{seed_number}.trees")

