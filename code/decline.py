import msprime
import os

seed_number = int(os.getenv("PBS_ARRAY_INDEX"))  # Find out the job number

### Parameters ###
Ne = 30000  # Effective population size
Ne_current = 300  # Population size at time 0 (the most recent time)
seqlength = 1e6  # Sequence length
mu = 3.75e-8  # Mutation rate
recombination_rate = 10 * mu  # Recombination rate
timepoints = list(range(0, 73, 3))  # Time points from 0 to 72 in steps of 3

### Define demographic model ###
demography = msprime.Demography()
demography.add_population(name="pop", initial_size=Ne_current)  # The most recent population size
demography.add_population_parameters_change(
    time=36,  # Time of population decline (36 generations ago)
    initial_size=Ne,  # Population size before crash
    population="pop"  # Apply to the population
)

samples = [msprime.SampleSet(50, population="pop",time=t) for t in timepoints]

# Generate individual names in vcf file
individual_names = []
for t in timepoints:
    if t >= 36:
        prefix = "post"
    else:
        prefix = "pre"
    for i in range(50):
        individual_names.append(f"{prefix}{t}_{i}")


### Simulate ancestry ###
ts = msprime.sim_ancestry(
    samples=samples,
    demography=demography,  # Use the defined demographic model
    sequence_length=seqlength,  # Genome length
    recombination_rate=recombination_rate,
    model=[
        msprime.StandardCoalescent(),
    ],
    random_seed=seed_number
)

# Add mutations
mts = msprime.sim_mutations(
    ts,  # Input tree sequence
    rate=mu,  # Mutation rate
    random_seed=seed_number
)



### Save tree sequence to file ###
mts.dump(f"decline_{seed_number}.trees")  # Save the tree sequence to a file

### Save to VCF file ###
with open(f"decline_{seed_number}.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)  # Write the tree sequence to a VCF file
