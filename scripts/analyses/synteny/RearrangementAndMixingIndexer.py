import pandas as pd
import os
import argparse
from math import log

"""
Extend the RearrangementIndexer from Lewin et al. (see https://github.com/symgenoevolab/RearrangementIndexer/blob/main/RearrangementIndexer.py)
to compute a Mixing Index based on a normalized Shannon entropy definition.
"""


def compute_shannon_entropy(probabilities):
    # Compute the Shannon entropy from a list of probabilities of k events
    logs = [log(p) if p > 0 else 0 for p in probabilities]
    return -sum([p * l for p, l in zip(probabilities, logs)])


def compute_region_probabilities(df, most_genes_chrom, alg, max_gap=0, freq_max=1.0, chrom_freqs=None):  
    # Compute the probabilites of a gene belonging to an ALG=alg to be present
    # in a region i on the chromomsome=most_genes_chrom.
    df_sorted = df[df['Chromosome'] == most_genes_chrom].sort_values(by=['Start'])

    if chrom_freqs is None:
        counts = df_sorted['ALG'].value_counts()
        chrom_freqs = (counts / len(df_sorted)).to_dict()

    all_regions, cur_region, gap = [], 0, 0
    for row in df_sorted.itertuples(index=False):
        if row.ALG == alg:
            cur_region += 1
            gap = 0
        elif cur_region > 0:
            gap += 1
            class_freq = chrom_freqs.get(row.ALG, 0.0)
            if gap > max_gap or class_freq >= freq_max:
                all_regions.append(cur_region)
                cur_region = 0
                gap = 0
    if cur_region > 0:
        all_regions.append(cur_region)
    tot_algs = sum(all_regions)
    return tot_algs, [r / tot_algs for r in all_regions]
    
    
# Define functions
def process_tsv_file(file_path, max_gap=0, freq_max=1.0):
    print(f"Processing file: {file_path} (max_gap={max_gap}, freq_max={freq_max})")
    
    # Load the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=None)
    
    # Rename columns for clarity
    df.columns = ["Index", "Status", "Chromosome", "Start", "End", "ALG"]
    
    # Combine ALG sub-parts into single ALGs
    df['ALG'] = df['ALG'].replace({'A1a': 'A1', 'A1b': 'A1', 'Ea': "E", 'Eb': "E", 'Qa': 'Q', 'Qb': 'Q', 'Qc': 'Q', 'Qd': 'Q'})
    
    # Calculate proportion of genes belonging to each ALG for each chromosome
    chrom_alg_counts = df.groupby(["Chromosome", "ALG"]).size().unstack(fill_value=0)
    
    # Calculate the total count of genes for each chromosome
    chrom_total_counts = chrom_alg_counts.sum(axis=1)
    
    # Calculate the proportion of genes belonging to each ALG for each chromosome
    chrom_alg_proportions = chrom_alg_counts.div(chrom_total_counts, axis=0)
    
    # Determine all unique ALGs
    all_algs = chrom_alg_counts.columns
    
    # Determine chromosome with the most genes from each ALG
    chrom_most_genes = chrom_alg_counts.idxmax(axis=0)
    
    # Calculate the rearrangement index for each ALG
    rearrangement_indices = {}
    proportion_on_chrom_values = {}
    total_proportion_on_chrom_values = {}
    mixing_on_chrom_values = {}
    
    for alg in all_algs:
        most_genes_chrom = chrom_most_genes[alg]
        print(f"Calculating metrics for ALG: {alg}")
        
        # Find the proportion of genes from the ALG on the chromosome with the most genes from that ALG
        if alg in chrom_alg_counts.columns:
            proportion_on_chrom = chrom_alg_counts.loc[most_genes_chrom, alg] / chrom_alg_counts[alg].sum()
        else:
            proportion_on_chrom = pd.NA
        proportion_on_chrom_values[alg] = proportion_on_chrom
        
        # Find the proportion of total genes on that chromosome that are from that ALG
        if alg in chrom_alg_counts.columns:
            total_proportion_on_chrom = chrom_alg_counts.loc[most_genes_chrom, alg] / chrom_total_counts[most_genes_chrom]
        else:
            total_proportion_on_chrom = pd.NA
        total_proportion_on_chrom_values[alg] = total_proportion_on_chrom
        
        # Compute the mixing of that ALG on that chromosome.
        chrom_freqs = chrom_alg_proportions.loc[most_genes_chrom].to_dict()
        tot_algs, probabilities = compute_region_probabilities(
            df, most_genes_chrom, alg, max_gap=max_gap, freq_max=freq_max, chrom_freqs=chrom_freqs
        )
        if tot_algs <= 1:
            mixing_on_chrom_values[alg] = 0.0
        else:
            mixing_on_chrom_values[alg] = abs(compute_shannon_entropy(probabilities) / log(tot_algs))
        
        # Calculate the rearrangement index for this ALG
        if alg in chrom_alg_counts.columns:
            rearrangement_index = 1 - (proportion_on_chrom * total_proportion_on_chrom)
        else:
            rearrangement_index = pd.NA
        rearrangement_indices[alg] = rearrangement_index
    
    # Convert dictionaries to pandas Series
    rearrangement_series = pd.Series(rearrangement_indices, name=os.path.basename(file_path))
    proportion_on_chrom_series = pd.Series(proportion_on_chrom_values, name=os.path.basename(file_path))
    total_proportion_on_chrom_series = pd.Series(total_proportion_on_chrom_values, name=os.path.basename(file_path))
    mixing_on_chrom_series = pd.Series(mixing_on_chrom_values, name=os.path.basename(file_path))
    
    return rearrangement_series, proportion_on_chrom_series, total_proportion_on_chrom_series, mixing_on_chrom_series

def main(input_dir, output_dir, max_gap=0, freq_max=1.0):
    # Get list of TSV files in the input directory
    tsv_files = [file for file in os.listdir(input_dir) if file.endswith(".tsv")]

    # Initialize DataFrames to store the outputs
    rearrangement_df = pd.DataFrame()
    proportion_on_chrom_df = pd.DataFrame()
    total_proportion_on_chrom_df = pd.DataFrame()
    mixing_on_chrom_df = pd.DataFrame()

    # Process each TSV file
    for tsv_file in tsv_files:
        tsv_path = os.path.join(input_dir, tsv_file)
        rearrangement_series, proportion_series, total_proportion_series, mixing_series = process_tsv_file(tsv_path, max_gap=max_gap, freq_max=freq_max)

        # Append to respective DataFrames
        rearrangement_df = pd.concat([rearrangement_df, rearrangement_series], axis=1)
        proportion_on_chrom_df = pd.concat([proportion_on_chrom_df, proportion_series], axis=1)
        total_proportion_on_chrom_df = pd.concat([total_proportion_on_chrom_df, total_proportion_series], axis=1)
        mixing_on_chrom_df = pd.concat([mixing_on_chrom_df, mixing_series], axis=1)

    # Output files suffix
    suffix = f"_gap{max_gap}_freq{freq_max}"

    # Save the rearrangement index to a table format
    rearrangement_df.to_csv(f"{output_dir}/Rearrangement_index{suffix}.tsv", sep="\t")
    print(f"Rearrangement indices saved to Rearrangement_index{suffix}.tsv")

    # Save the proportion_on_chrom values to a table format
    proportion_on_chrom_df.to_csv(f"{output_dir}/Splitting_parameter{suffix}.tsv", sep="\t")
    print(f"Fission parameters saved to Splitting_parameter{suffix}.tsv")

    # Save the total_proportion_on_chrom values to a table format
    total_proportion_on_chrom_df.to_csv(f"{output_dir}/Combining_parameter{suffix}.tsv", sep="\t")
    print(f"Fusion parameters values saved to Combining_parameter{suffix}.tsv")

    # Save the mixing_on_chrom_values tot a table format
    mixing_on_chrom_df.to_csv(f"{output_dir}/Mixing_index{suffix}.tsv", sep="\t")
    print(f"Mixing indices values saved to Mixing_index{suffix}.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Rearrangement and Mixing indices.")
    parser.add_argument('--in_dir', required=True, help='Path to the directory containing coordinates files with information about genes genomic location and ALG.')
    parser.add_argument('--out_dir', required=False, default='.', help='Path to output directory.')
    parser.add_argument('--max_gap', required=False, default=0, type=int, help='Skip gap(s) <= max_gap between 2 consecutive regions of an ALG.')
    parser.add_argument('--max_freq', required=False, default=1.0, type=float, help='Skip gap(s) between consecutive regions of an ALG x if each ALG y != x present in a gap has a frequency < max_freq.')
    args = parser.parse_args() 
    main(args.in_dir, args.out_dir, args.max_gap, args.max_freq)
    
