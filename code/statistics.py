import tskit
import pandas as pd
import numpy as np
import glob


# Time points from 0 to 72 in steps of 3
timepoints = list(range(0, 73, 3))

# Get the prefix of the tree file
def get_tree_prefixes():
    prefixes = []
    for f in glob.glob("*.trees"):
        prefixes.append(f.replace(".trees", ""))
    return prefixes

prefix = get_tree_prefixes()

# Calculate summary statistics
def oneway_stats(ts):
    stats = {
        "Number of trees": ts.num_trees,
        "Number of mutations": ts.num_mutations,
        "Density of segregating sites": ts.segregating_sites(),
        "Nucleotide diversity (pi)": ts.diversity(),
        "Tajima's D": ts.Tajimas_D(),
    }
    return stats

def oneway_output(ts, prefix):
    results = {}
    for t in timepoints:
        # Extract samples for the current time point
        sample_ids = [n.id for n in ts.nodes() if n.time == t]
        sub_ts = ts.simplify(samples=sample_ids)

        # Calculate statistics
        stats = oneway_stats(sub_ts)
        results[t] = stats

    # Convert results to a pandas DataFrame
    oneway_df = pd.DataFrame(results).T  # Transpose to have time points as rows
    oneway_df.index.name = "t"
    oneway_df.to_csv(f"{prefix}.csv")

def multiway_stats(ts, sample_id_dict):
    fst_array = []
    dxy_array = []
    gr_array = []

    for t1 in sample_id_dict.keys():
        if 36 < t1 <= 72:
            fst_row = []
            dxy_row = []
            gr_row = []
            for t2 in sample_id_dict.keys():
                if 0 <= t2 < 36:
                    sample_set_1 = sample_id_dict[t1]
                    sample_set_2 = sample_id_dict[t2]
                    fst = ts.Fst([sample_set_1, sample_set_2])
                    dxy = ts.divergence([sample_set_1, sample_set_2])
                    gr = ts.genetic_relatedness([sample_set_1, sample_set_2])
                    fst_row.append(fst)
                    dxy_row.append(dxy)
                    gr_row.append(gr)
            fst_array.append(fst_row)
            dxy_array.append(dxy_row)
            gr_array.append(gr_row)
            
    fst_array = np.array(fst_array)
    dxy_array = np.array(dxy_array)
    gr_array = np.array(gr_array)
    return fst_array, dxy_array, gr_array


def process_trees_file():
    for prefix in get_tree_prefixes():
        ts = tskit.load(f"{prefix}.trees")
        sample_id_dict = {t: ts.samples(time=t) for t in timepoints}
        
        oneway_output(ts, prefix)
        
        fst, dxy, gr = multiway_stats(ts, sample_id_dict)
        np.savez(f"{prefix}.npz", fst=fst, dxy=dxy, gr=gr)


if __name__ == "__main__":
    process_trees_file()