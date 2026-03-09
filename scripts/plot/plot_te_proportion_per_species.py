import sys
import re
import json
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot_repeat_landscape import NAME_MAP, GRAPH_LABELS


"""
Plot the repeat type proportion per species.
"""


#--------------------- HELPERS ---------------------#

RE_TYPES = ['DNA', 'RC', 'PLE', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple_repeat', 'Low_complexity', 'Unclassified']
TE_REPEATS = ['DNA', 'RC', 'PLE', 'LTR', 'LINE', 'SINE']

TE_REGEX = re.compile(r'^#{3}\w+')

GRAPH_LABELS_TE = GRAPH_LABELS.copy()
del GRAPH_LABELS_TE['Unknown']
del GRAPH_LABELS_TE['Other']

GRAPH_LABELS.update({
    'Satellite': '#7f7f7f',
    'Simple_repeat': '#999999',
    'Low_complexity': '#b3b3b3',
    'Unclassified': '#cccccc'})    


def get_chr_sizes(fai: str) -> dict:
    """Return a dictionary of chromosome sizes from a .fai file.

    Args:
        fai (str): Path to the .fai file containing chromosome sizes.

    Returns:
        dict: A dictionary with chromosome names as keys and their sizes as values.
    """
    with open(fai, 'r') as f:
        lines = f.readlines()
    return {line.split('\t')[0]: int(line.split('\t')[1]) for line in lines}


def parse_rm2b(fpath: str, chr_size: dict) -> tuple:
    """Parse the RepeatMasker *.BED file to compute the proportion of each repeat type.

    Args:
        fpath (str): Path to the RM *.BED file.
        chr_size (dict): Dictionary of chromosome sizes.
    Returns:
        tuple: Two dictionaries containing the proportion of each repeat type
               relative to the genome size and relative to total TE base pairs.
    """
    print(f'Parsing {fpath}...')
    
    summary_re, summary_te = {}, {}
    df = pd.read_csv(fpath, delimiter='\t', header=None)
    df = df[[0, 1, 2, 4, 6, 7]]
    df.columns = ['contig', 'start', 'end', 'length', 'class', 'family']
    
    # format class/family names
    df['class'] = df['class'].replace('Unknown', 'Unclassified')
    df.loc[df['family'].str.contains('Penelope', na=False), 'class'] = 'PLE'
    
    # Loop through chromosomes, keep assigned positions to prevent overlaps
    for chr in chr_size.keys():
        assigned = [False] * chr_size[chr]
        for _, row in df[df['contig'] == chr].sort_values(by='start').iterrows():
            # Compute real length base on not assigned bases
            length = 0
            for i in range(row['start'], row['end']):
                if not assigned[i]:
                    length += 1
                    assigned[i] = True
            if row['class'] in TE_REPEATS:
                summary_te[row['class']] = summary_te.get(row['class'], 0) + length
            if row['class'] in RE_TYPES:
                summary_re[row['class']] = summary_re.get(row['class'], 0) + length
            
    # Compute percentage for repetitive elements
    genome_size = sum(chr_size.values())
    for re, bp in summary_re.items():
        summary_re[re] = (bp / genome_size) * 100
    # Compute percentage for TEs only
    total_te_bp = sum([bp for _, bp in summary_te.items()])
    for te, bp in summary_te.items():
        summary_te[te] = (bp / total_te_bp) * 100
    return summary_re, summary_te
             
            
#----------------- PLOTTING FUNCTION ------------------#    
    
def plot_te_proportion(sp_details_files: list, rm2b_files: list, fai_files: list, output_file: str) -> None:
    """Plot the proportion of different TE types per species.
    
    Args:
        sp_details_files (list): List of paths to species details JSON files.
        rm2b_files (list): List of paths to RepeatMasker .bed files.
        fai_files (list): List of paths to genome .fai files.
        output_file (str): Path to save the output plot.
    """
    
    # Get species names and genome sizes
    samples = []
    for sp_details_file in sp_details_files:
        with open(sp_details_file, 'rt') as f:
            sp_details = json.load(f)
        samples.append(sp_details['scientific_name'])
        
    # Get repeat proportions
    all_data, te_data = [], []
    for rm2b, fai in zip(rm2b_files, fai_files):
        ad, td = parse_rm2b(rm2b, get_chr_sizes(fai))
        all_data.append([ad.get(k, 0) for k in RE_TYPES])
        te_data.append([td.get(k, 0) for k in GRAPH_LABELS_TE])
            
    all_data = np.array(all_data)
    te_data = np.array(te_data)    
        
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 0.33 * len(samples)))

    # Left panel: Stacked horizontal bar chart
    left = np.zeros(len(samples))
    for i, repeat_type in enumerate(RE_TYPES):
        ax1.barh(samples, all_data[:, i], left=left, 
                color=GRAPH_LABELS[repeat_type], label=repeat_type, edgecolor='none')
        left += all_data[:, i]

    ax1.set_xlabel('Proportion of genome (%)', fontsize=12)
    ax1.set_xlim(0, 80)
    ax1.invert_yaxis()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(left=False)

    # Right panel: Stacked horizontal bar chart (normalized to 100%)
    # Only show TE repeats
    left = np.zeros(len(samples))
    for i, repeat_type in enumerate(GRAPH_LABELS_TE):
        ax2.barh(samples, te_data[:, i], left=left, 
                color=GRAPH_LABELS[repeat_type], edgecolor='none')
        left += te_data[:, i]

    ax2.set_xlabel('Relative proportion of total TEs (%)', fontsize=12)
    ax2.set_xlim(0, 100)
    ax2.set_yticklabels([])
    ax2.invert_yaxis()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(left=False)

    # Create legend 
    handles = [plt.Rectangle((0,0),1,1, color=GRAPH_LABELS[rt]) for rt in RE_TYPES]
    ax2.legend(handles, RE_TYPES, loc='center left', bbox_to_anchor=(1.02, 0.5),
            frameon=True, title='Repeat type', fontsize=10, title_fontsize=11, 
            facecolor='#f2f2f2', edgecolor='none')
    
    fig.suptitle('Genomic repeat content and diversity', fontsize=16)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()    

    

# Test  
species = ['a_vaga', 'a_sp_wild', 'a_ricciae', 'r_rotatoria', 'm_quadricornifera', 'h_sp_wild', 'p_roseola', 's_nebaliae', 'n_agilis', 'b_calyciflorus', 'f_enflata']
sp_details = ['/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/genomes/{}/{}.details.json'.format(sp, sp) for sp in species]
rm2b_files = ['/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/genomes/{}/te_annotation/{}_rm.bed'.format(sp, sp) for sp in species]
fai_files = ['/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/genomes/{}/{}.assembly.fa.fai'.format(sp, sp) for sp in species]
output_file = '/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/analyses/te/te_proportion_per_species.svg'
plot_te_proportion(sp_details, rm2b_files, fai_files, output_file)