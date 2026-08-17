import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.colors import LogNorm

""" 
Compute and plot an inter-species colinearity similarity Jaccard-like index. 
For each species pair, this is calculated as the sum of their genes 
anchored in colinear regions divided by the total number of genes in both species.
"""

plt.rcParams['svg.fonttype'] = 'none'

#------------------- CONSTANTS -----------------#

BLOCK_RE = re.compile(r'^## Alignment\s+\d+:.*?N=(\d+)\s+(\S+)&(\S+)')
GENE_RE  = re.compile(r'^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)')
MATCH_SIZE = 5 # MCScanX parameter for 'Number of genes required to call a collinear block'
SPECIES = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}
MCSCANX_GFF = 'rotifers_mcscanx.gff'
MCSCANX_COLINEARITY = 'rotifers_mcscanx.collinearity'

OUT_JACCARD = 'bdelloid_jaccard_indices.tsv'
OUT_FRAG = 'bdelloid_fragmentation_indices.tsv'
OUT_PLOT = 'bdelloid_jaccard_and_frag_indices.svg'

#------------------- PARSERS -------------------#

def get_sp_genes(mcscanx_gff, species):
    """
    Parse a MCScanX gff input file and extract all species genes.
    """
    sp_genes = {k: set() for k in species}
    df = pd.read_csv(mcscanx_gff, sep='\t', header=None, 
                     names=['chr', 'gene', 'start', 'end'])
    for row in df.itertuples(index=False):
        chr = row.chr[:2]
        if chr in sp_genes:
            sp_genes[chr].add(row.gene)
    return sp_genes


def parse_collinearity(mcscanx_colinearity):
    """
    Parse a MCScanX colinearity file to extract anchored genes
    between pair of species and the the sizes of the colinear regions.
    """
    anchored = defaultdict(lambda: {'spA': set(), 'spB': set()})
    block_sizes = defaultdict(list)
    current_pair = None
    
    with open(mcscanx_colinearity) as f:
        for line in f:
            m = BLOCK_RE.match(line)
            
            # --- block header ---
            if m:
                n_pairs = int(m.group(1))
                chrA = m.group(2)
                chrB = m.group(3)
                spA = chrA[:2]
                spB = chrB[:2]
                
                if spA == spB:
                    # skip intra-species colinear regions
                    current_pair = None
                elif chrA[2:] != chrB[2:]:
                    # skip regions between non-orthologous chromosomes
                    current_pair = None
                else:
                    pair_key     = tuple(sorted([spA, spB]))
                    current_pair = pair_key
                    block_sizes[pair_key].append(n_pairs)
                 
            # --- gene pair line ---
            else:
                m = GENE_RE.match(line)
                if m and current_pair is not None:
                    geneA = m.group(1)
                    geneB = m.group(2)
                    anchored[current_pair]['spA'].add(geneA)
                    anchored[current_pair]['spB'].add(geneB)

    return anchored, dict(block_sizes)


#---------------- JACCARD INDEX ----------------#

def compute_jaccard(anchored, block_sizes, total_genes):
    """
    Compute the Jaccard-like index for each species pair and
    a fragmentation index defined as the min number of gene in a block
    divided by the block sizes average.
    """
    species = sorted(total_genes.keys())
    j_matrix = pd.DataFrame(np.nan, index=species, columns=species)
    frag_matrix = pd.DataFrame(np.nan, index=species, columns=species)
    
    for pair_key, gene_sets in anchored.items():
        spA, spB = pair_key

        anch_A = len(gene_sets['spA'])
        anch_B = len(gene_sets['spB'])
        tot_A = len(total_genes.get(spA))
        tot_B = len(total_genes.get(spB))

        if tot_A is None or tot_B is None:
            print(f'Missing total for {spA} or {spB}, skipping...')
            continue
        
        sizes = block_sizes.get(pair_key, [])
        n_blocks = len(sizes)
        total_pairs = sum(sizes)

        mean_block_size = total_pairs / n_blocks if n_blocks > 0 else np.nan
        frag = MATCH_SIZE / mean_block_size if mean_block_size else np.nan

        jaccard = (anch_A + anch_B) / (tot_A + tot_B)

        j_matrix.loc[spA, spB] = jaccard
        j_matrix.loc[spB, spA] = jaccard
        frag_matrix.loc[spA, spB] = frag
        frag_matrix.loc[spB, spA] = frag

    return j_matrix, frag_matrix


def plot_matrices(j_matrix, frag_matrix, species, out_file):
    """
    Plot the Jaccard and Fragementation index matrices.
    """
    order  = [sp for sp in species if sp in j_matrix.index]
    order += [sp for sp in j_matrix.index if sp not in order]
    labels = [species.get(sp, sp) for sp in order]

    j_matrix    = j_matrix.loc[order, order]
    frag_matrix = frag_matrix.loc[order, order]

    sns.set_theme(style='white')
    _, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(len(order) * 2.4 + 3, len(order) * 1.2 + 0.5)
    )

    # --- Left: Jaccard ---
    sns.heatmap(
        j_matrix.astype(float), ax=ax1,
        annot=True, fmt='.3f', cmap='YlOrRd',
        vmin=0, vmax=1,
        linewidths=0.5, linecolor='white',
        xticklabels=labels, yticklabels=labels,
        cbar_kws={'label': 'Jaccard synteny index', 'shrink': 0.7},
    )
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right',
                        fontstyle='italic', fontsize=12)
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0,
                        fontstyle='italic', fontsize=12)
    ax1.set_title('Pairwise genome-wide Jaccard synteny index',
                  fontsize=15, pad=12)

    # --- Right: fragmentation ---
    sns.heatmap(
        frag_matrix.astype(float), ax=ax2,
        annot=True, fmt='.2f', cmap='Blues',
        linewidths=0.5, linecolor='white',
        xticklabels=labels, yticklabels=labels,
        cbar_kws={'label': 'Fragmentation index', 'shrink': 0.7},
    ) 
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right',
                        fontstyle='italic', fontsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0,
                        fontstyle='italic', fontsize=12)
    ax2.set_title('Synteny fragmentation\n(lower = larger, more contiguous blocks)',
                  fontsize=15, pad=12)

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    
    
#-------------------- RUN ----------------------#

anchored, block_sizes = parse_collinearity(MCSCANX_COLINEARITY)
j_matrix, f_matrix = compute_jaccard(anchored, block_sizes, get_sp_genes(MCSCANX_GFF, SPECIES))
j_matrix.to_csv(OUT_JACCARD, sep='\t')
f_matrix.to_csv(OUT_FRAG, sep='\t')
plot_matrices(j_matrix, f_matrix, SPECIES, OUT_PLOT)
