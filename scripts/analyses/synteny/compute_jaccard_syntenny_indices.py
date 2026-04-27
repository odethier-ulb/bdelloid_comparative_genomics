# DISCLAIMER:
# This script was initially generated with assistance from Claude (Anthropic)
# and subsequently reviewed, modified, and adapted by the authors.
# The authors take full responsibility for the final code, its accuracy,
# implementation, and any resulting analyses.

"""
Compute pairwise genome-wide Jaccard synteny indices from a single
all-vs-all MCScanX .collinearity file.

For each species pair (A, B):

    J(A, B) = (|anchored_A| + |anchored_B|) / (|total_A| + |total_B|)

where anchored_X = genes of species X present in at least one syntenic
block with species Y (deduplicated).

Outputs
-------
- TSV matrix of pairwise Jaccard indices
- TSV of per-species-pair anchoring statistics (raw numbers)
- Heatmap SVG/PDF
"""

import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'


# ------------------------------------------------------------------ #
#                         PARSER
# ------------------------------------------------------------------ #

HOMOELOGS = set([(1, 4), (2, 5), (3, 6)])

BLOCK_RE = re.compile(
    r'^## Alignment\s+\d+:.*?N=(\d+)\s+(\S+)&(\S+)'
)
GENE_RE = re.compile(
    r'^\s+\d+-\s*\d+:\s+(\S+)\s+(\S+)'
)

GENE_RE = re.compile(r'^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)')


def species_from_chrom(chrom: str) -> str:
    """Return the 2-letter species code from a chromosome ID (e.g. 'ar1' → 'ar')."""
    return chrom[:2]


def parse_collinearity(collinearity_file: str):
    """
    Parse an all-vs-all MCScanX .collinearity file.

    Returns
    -------
    anchored : dict  {species_pair_tuple: {'spA': set_of_genes, 'spB': set_of_genes}}
        species_pair_tuple is always (min, max) alphabetically so A-B == B-A.
    block_stats : list of dicts
        One entry per block with keys: chrA, chrB, spA, spB, n_pairs.
    """
    anchored    = defaultdict(lambda: {'spA': set(), 'spB': set()})
    block_stats = []
    current_pair = None
    
    with open(collinearity_file) as fh:
        for line in fh:
            m = BLOCK_RE.match(line)
            # --- block header ---
            if m:
                n_pairs = int(m.group(1))
                chrA    = m.group(2)
                chrB    = m.group(3)
                spA     = species_from_chrom(chrA)
                spB     = species_from_chrom(chrB)
                
                if spA == spB:
                    current_pair = None  # skip self-vs-self blocks       
                elif tuple(sorted((int(chrA[2:]), int(chrB[2:])))) in HOMOELOGS:
                    current_pair = None  # skip known homoeologous pairs
                else:     
                    pair_key     = tuple(sorted([spA, spB]))
                    current_pair = pair_key

                    block_stats.append({
                        'chrA': chrA, 'chrB': chrB,
                        'spA': spA,   'spB': spB,
                        'n_pairs': n_pairs,
                    })
            # --- gene pair line ---
            else:
                m = GENE_RE.match(line)
                if m and current_pair is not None:
                    geneA = m.group(1)
                    geneB = m.group(2)
                    anchored[current_pair]['spA'].add(geneA)
                    anchored[current_pair]['spB'].add(geneB)
         
    return anchored, block_stats


# ------------------------------------------------------------------ #
#                       JACCARD CALCULATION
# ------------------------------------------------------------------ #

def compute_jaccard(anchored: dict, total_genes: dict) -> pd.DataFrame:
    """
    Parameters
    ----------
    anchored    : output of parse_collinearity()
    total_genes : dict {species_code: int}  total gene count per species

    Returns
    -------
    matrix : pd.DataFrame  symmetric NxN Jaccard matrix (diagonal = within-species)
    stats  : pd.DataFrame  full per-pair statistics table
    """
    species = sorted(total_genes.keys())
    matrix = pd.DataFrame(np.nan, index=species, columns=species)
    rows   = []

    for pair_key, gene_sets in anchored.items():
        spA, spB = pair_key          # already sorted

        anch_A = len(gene_sets['spA'])
        anch_B = len(gene_sets['spB'])
        tot_A  = total_genes.get(spA, None)
        tot_B  = total_genes.get(spB, None)

        if tot_A is None or tot_B is None:
            print(f'[WARN] Missing total gene count for {spA} or {spB}, skipping pair.')
            continue

        jaccard = (anch_A + anch_B) / (tot_A + tot_B)

        matrix.loc[spA, spB] = jaccard
        matrix.loc[spB, spA] = jaccard

        rows.append({
            'species_A'  : spA,
            'species_B'  : spB,
            'anchored_A' : anch_A,
            'anchored_B' : anch_B,
            'total_A'    : tot_A,
            'total_B'    : tot_B,
            'pct_anch_A' : round(anch_A / tot_A * 100, 2),
            'pct_anch_B' : round(anch_B / tot_B * 100, 2),
            'jaccard'    : round(jaccard, 4),
        })

    stats = pd.DataFrame(rows).sort_values(['species_A', 'species_B'])
    return matrix, stats


# ------------------------------------------------------------------ #
#                           HEATMAP
# ------------------------------------------------------------------ #

def plot_jaccard_heatmap(matrix: pd.DataFrame, species_labels: dict,
                          out_file: str):
    """
    Parameters
    ----------
    matrix        : NxN Jaccard DataFrame (species codes as index/columns)
    species_labels: dict {code: full_name} for axis labels (optional)
    out_file      : output path (.svg or .pdf)
    """
    # reorder rows/cols to match a meaningful species order if provided
    order  = [sp for sp in species_labels if sp in matrix.index]
    order += [sp for sp in matrix.index if sp not in order]
    matrix = matrix.loc[order, order]

    labels = [species_labels.get(sp, sp) for sp in order]

    sns.set_theme(style='white')
    _, ax = plt.subplots(figsize=(len(order) * 1.2 + 1.5, len(order) * 1.2 + 0.5))

    sns.heatmap(
        matrix.astype(float),
        ax=ax,
        annot=True, fmt='.3f',
        cmap='YlOrRd',
        vmin=0, vmax=1,
        linewidths=0.5, linecolor='white',
        xticklabels=labels,
        yticklabels=labels,
        cbar_kws={'label': 'Jaccard synteny index', 'shrink': 0.7},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right',
                        fontstyle='italic', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,
                        fontstyle='italic', fontsize=10)
    ax.set_title('Pairwise genome-wide Jaccard synteny index', fontsize=13, pad=12)

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    print(f'Saved → {out_file}')


# ------------------------------------------------------------------ #
#                           EXECUTION
# ------------------------------------------------------------------ #

gff_file = 'rotifers_mcscanx.gff'
collinearity_file = 'rotifers_mcscanx.collinearity'
out_matrix  = 'jaccard_matrix.tsv'
out_stats   = 'jaccard_stats.tsv'
out_heatmap = 'jaccard_heatmap.svg'

# Compute gene count per species from the GFF
total_genes = {
    'av': 0,
    'as': 0,
    'ar': 0,
    'rr': 0,
    'mq': 0,
    'hs': 0,   
    'pr': 0, 
}
with open(gff_file) as f:
    for line in f:
        lsplt = line.strip().split('\t')
        if len(lsplt) < 3:
            continue
        chrom = lsplt[0]
        sp = species_from_chrom(chrom)
        if sp in total_genes:
            total_genes[sp] += 1
            
# Full species names for heatmap labels (order also controls row/col order)
species_labels = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}

print('Parsing collinearity file...')
anchored, block_stats = parse_collinearity(collinearity_file)
print(f'  {len(block_stats)} blocks parsed across {len(anchored)} species pairs.')

matrix, stats = compute_jaccard(anchored, total_genes)

matrix.to_csv(out_matrix, sep='\t')
stats.to_csv(out_stats, sep='\t', index=False)
print(f'Saved → {out_matrix}')
print(f'Saved → {out_stats}')

print('\nJaccard matrix:')
print(matrix.to_string())

plot_jaccard_heatmap(matrix, species_labels, out_heatmap)

