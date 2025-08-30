import msprime
import os

seed_number = int(os.getenv("PBS_ARRAY_INDEX"))  # Find out the job number

### Parameters ###
Ne = 30000  # Effective population size
seqlength = 1e6  # Sequence length
mu = 3.75e-8  # Mutation rate
recombination_rate = 10 * mu  # Recombination rate
timepoints = list(range(0, 73, 3))  # Time points from 0 to 72 in steps of 3
samples = [msprime.SampleSet(50, time=t) for t in timepoints]


### Simulate ancestry ###
ts = msprime.sim_ancestry(
    samples=samples,
    population_size=Ne,
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
mts.dump(f"constant_{seed_number}.trees")