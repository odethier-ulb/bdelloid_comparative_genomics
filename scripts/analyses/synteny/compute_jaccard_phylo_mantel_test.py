# DISCLAIMER:
# This script was initially generated with assistance from Claude (Anthropic)
# and subsequently reviewed, modified, and adapted by the authors.
# The authors take full responsibility for the final code, its accuracy,
# implementation, and any resulting analyses.

"""
Test whether the pairwise Jaccard synteny index correlates with phylogenetic
(patristic) distance across bdelloid species.
"""

import itertools

import numpy as np
import pandas as pd
import dendropy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import squareform

plt.rcParams['svg.fonttype'] = 'none'


# ------------------------------------------------------------------ #
#                    1. PATRISTIC DISTANCES
# ------------------------------------------------------------------ #

def extract_patristic(treefile, label_map):
    """Return a DataFrame of pairwise patristic distances for focal taxa."""
    tree = dendropy.Tree.get(path=treefile, schema='newick')
    pdm  = tree.phylogenetic_distance_matrix()

    # index taxa by label
    taxa = {t.label: t for t in tree.taxon_namespace}

    codes   = sorted(label_map.values())
    code2label = {v: k for k, v in label_map.items()}

    matrix = pd.DataFrame(0.0, index=codes, columns=codes)

    for cA, cB in itertools.combinations(codes, 2):
        tA = taxa[code2label[cA]]
        tB = taxa[code2label[cB]]
        d  = pdm(tA, tB)
        matrix.loc[cA, cB] = d
        matrix.loc[cB, cA] = d

    return matrix


# ------------------------------------------------------------------ #
#                       2. MANTEL TEST
# ------------------------------------------------------------------ #

def mantel_test(mat1, mat2, n_permutations=9999, seed=42):
    """
    Mantel test between two symmetric distance matrices.

    Parameters
    ----------
    mat1, mat2 : pd.DataFrame  symmetric NxN matrices with same index/columns
    n_permutations : int

    Returns
    -------
    r : float   Pearson correlation between upper-triangle vectors
    p : float   permutation p-value (one-tailed, r >= observed)
    """
    rng = np.random.default_rng(seed)

    # align
    order = sorted(mat1.index)
    m1 = mat1.loc[order, order].values
    m2 = mat2.loc[order, order].values

    # upper-triangle vectors
    idx = np.triu_indices(len(order), k=1)
    v1  = m1[idx]
    v2  = m2[idx]

    observed_r, _ = stats.pearsonr(v1, v2)

    # permute rows/cols of m1 simultaneously
    count = 0
    for _ in range(n_permutations):
        perm = rng.permutation(len(order))
        m1_p = m1[np.ix_(perm, perm)]
        r_p, _ = stats.pearsonr(m1_p[idx], v2)
        if r_p >= observed_r:
            count += 1

    p_value = (count + 1) / (n_permutations + 1)
    return observed_r, p_value


# ------------------------------------------------------------------ #
#                       3. SCATTER PLOT
# ------------------------------------------------------------------ #

def abbrev(sp_code):
    """Return a short italic-friendly label: 'A. vaga' → 'A. vaga', etc."""
    return SPECIES_NAMES[sp_code]


def plot_scatter(patristic_df, jaccard_df, out_file):
    codes = sorted(set(patristic_df.index) & set(jaccard_df.index))

    rows = []
    for cA, cB in itertools.combinations(codes, 2):
        rows.append({
            'spA'     : cA,
            'spB'     : cB,
            'distance': patristic_df.loc[cA, cB],
            'jaccard' : jaccard_df.loc[cA, cB],
        })
    df = pd.DataFrame(rows).dropna(subset=['jaccard'])

    sns.set_theme(style='white')
    _, ax = plt.subplots(figsize=(9, 6.5))

    # single colour for all points
    ax.scatter(
        df['distance'], df['jaccard'],
        color='#2f5597', s=70, alpha=0.85,
        edgecolors='white', linewidths=0.5, zorder=3,
    )

    # regression line
    slope, intercept, r, p, _ = stats.linregress(df['distance'], df['jaccard'])
    x_line = np.linspace(df['distance'].min(), df['distance'].max(), 200)
    ax.plot(x_line, slope * x_line + intercept,
            color='#555555', linewidth=1.2, linestyle='--', zorder=2,
            label=f'Linear fit (r = {r:.3f}, p = {p:.4f})')

    # label each point: "A. vaga – A. ricciae" on two lines
    for _, row in df.iterrows():
        label = f"{abbrev(row['spA'])} – {abbrev(row['spB'])}"
        ax.annotate(
            label,
            xy=(row['distance'], row['jaccard']),
            xytext=(5, 3), textcoords='offset points',
            fontsize=6, color='#333333',
            fontstyle='italic',
        )

    ax.set_xlabel('Patristic distance (substitutions per site)', fontsize=11)
    ax.set_ylabel('Jaccard synteny index', fontsize=11)
    ax.set_title('Synteny conservation vs. phylogenetic distance', fontsize=12)
    ax.legend(frameon=True, facecolor='#f2f2f2', edgecolor='none', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    print(f'Saved → {out_file}')
    return df


# ------------------------------------------------------------------ #
#                           EXECUTION
# ------------------------------------------------------------------ #


treefile    = 'partitions.txt.treefile'
jaccard_tsv = 'jaccard_matrix.tsv'
out_scatter = 'jaccard_vs_distance.svg'
out_mantel  = 'mantel_results.tsv'

# Mapping: exact taxon label in tree → 2-letter Jaccard matrix code
LABEL_MAP = {
    'a ricciae'        : 'ar',
    'a sp'             : 'as',
    'a vaga'           : 'av',
    'h sp'             : 'hs',
    'm quadricornifera': 'mq',
    'r rotatoria'      : 'rr',
    'p roseola'        : 'pr',
}

# Full species names for plot labels
SPECIES_NAMES = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}


print('Extracting patristic distances...')
patristic_df = extract_patristic(treefile, LABEL_MAP)
print(patristic_df.to_string())

print('\nLoading Jaccard matrix...')
jaccard_df = pd.read_csv(jaccard_tsv, sep='\t', index_col=0)
# keep only focal species present in both matrices
common = sorted(set(patristic_df.index) & set(jaccard_df.index))
patristic_df = patristic_df.loc[common, common]
jaccard_df   = jaccard_df.loc[common, common]

# convert Jaccard to dissimilarity for Mantel test
dissim_df = 1 - jaccard_df

print('\nRunning Mantel test (Jaccard dissimilarity vs patristic distance)...')
r, p = mantel_test(patristic_df, dissim_df, n_permutations=9999)
print(f'  Mantel r = {r:.4f}')
print(f'  p-value  = {p:.4f}  (9999 permutations)')

mantel_results = pd.DataFrame([{
    'mantel_r'      : round(r, 4),
    'p_value'       : round(p, 4),
    'n_permutations': 9999,
    'metric'        : 'Jaccard dissimilarity (1 - Jaccard) vs patristic distance',
}])
mantel_results.to_csv(out_mantel, sep='\t', index=False)
print(f'Saved → {out_mantel}')

print('\nPlotting...')
pair_df = plot_scatter(patristic_df, jaccard_df, out_scatter)
