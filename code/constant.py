import msprime
import os

seed_number = int(os.getenv("PBS_ARRAY_INDEX"))  # Find out the job number

### Parameters ###
Ne = 594  # Effective population size (harmonic mean of 300 and 30000)
seqlength = 1e6  # Sequence length
mu = 3.75e-8  # Mutation rate
recombination_rate = 10 * mu  # Recombination rate
timepoints = list(range(0, 73, 3))  # Time points from 0 to 72 in steps of 3
samples = [msprime.SampleSet(50, time=t) for t in timepoints]

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
    population_size=Ne,
    sequence_length=seqlength,  # Genome length
    recombination_rate=recombination_rate,
    model=[
        msprime.StandardCoalescent()
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
mts.dump(f"constant_{seed_number}.trees")  # Save the tree sequence to a file

### Save to VCF file ###
with open(f"constant_{seed_number}.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)  # Write the tree sequence to a VCF file