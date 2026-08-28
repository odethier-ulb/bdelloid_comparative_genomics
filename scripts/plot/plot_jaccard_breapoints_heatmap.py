import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
Plot heatmaps of Jaccard synteny indices and breakpoints per
average genome size for each pair of species.
"""

plt.rcParams['svg.fonttype'] = 'none'

# ------------------- CONSTANTS ----------------- #

SPECIES = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}
GENOME_SIZE_MB = {
    'av': 101,   
    'as': 123,  
    'ar': 105,  
    'rr': 183,   
    'mq': 278,   
    'hs': 111,   
    'pr': 90,   
}

IN_DIR = './jaccard'
JACCARD_MATRIX_FILE     = f'{IN_DIR}/bdelloid_jaccard_matrix.tsv'
BREAKPOINTS_MATRIX_FILE = f'{IN_DIR}/bdelloid_n_breakpoints_matrix.tsv'

OUT_PLOT = f'{IN_DIR}/bdelloid_jaccard_breakpoint_heatmaps.svg'
J_MIN = 0.65
J_MAX = 1
B_MIN = 0.35
B_MAX = 4

# ------------------- LOADERS -------------------- #

def get_matrix(matrix_file):
    """Load a symmetric species x species matrix TSV (diagonal = NaN)."""
    return pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)


def compute_breakpoint_rate(breakpoints_matrix, genome_size_mb):
    """
    Convert the number of breakpoints between pair of species to
    a rate per Mb by dividing the value by their average genome size.
    """
    rate_matrix = pd.DataFrame(
        np.nan, index=breakpoints_matrix.index, columns=breakpoints_matrix.columns
    )
    for sp1 in breakpoints_matrix.index:
        for sp2 in breakpoints_matrix.columns:
            if sp1 == sp2:
                continue
            n_bp = breakpoints_matrix.loc[sp1, sp2]
            if pd.isna(n_bp):
                continue
            combined_size = (genome_size_mb[sp1] + genome_size_mb[sp2]) / 2
            rate_matrix.loc[sp1, sp2] = n_bp / combined_size
    return rate_matrix


# -------------------- PLOTTING ------------------- #

def plot_matrices(j_matrix, rate_matrix, species, out_file):
    """
    Plot the Jaccard index and breakpoint rate matrices side by side.
    """
    order  = [sp for sp in species if sp in j_matrix.index]
    order += [sp for sp in j_matrix.index if sp not in order]
    labels = [species.get(sp, sp) for sp in order]

    j_matrix    = j_matrix.loc[order, order]
    rate_matrix = rate_matrix.loc[order, order]

    sns.set_theme(style='white')
    _, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(len(order) * 2.4 + 3, len(order) * 1.2 + 0.5)
    )

    # --- Left: Jaccard --- 
    sns.heatmap(
        j_matrix.astype(float), ax=ax1,
        annot=True, fmt='.3f', cmap='YlOrRd',
        vmin=J_MIN, vmax=J_MAX,
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

    # --- Right: breakpoint rate ---
    sns.heatmap(
        rate_matrix.astype(float), ax=ax2,
        annot=True, fmt='.3f', cmap='Blues',
        linewidths=0.5, linecolor='white',
        vmin=B_MIN, vmax=B_MAX,
        xticklabels=labels, yticklabels=labels,
        cbar_kws={'label': 'Breakpoints per Mbp', 'shrink': 0.7},
    )
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right',
                        fontstyle='italic', fontsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0,
                        fontstyle='italic', fontsize=12)
    ax2.set_title('Breakpoints per Mbp',
                  fontsize=15, pad=12)

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')


# ----------------------- RUN --------------------- #

j_matrix = get_matrix(JACCARD_MATRIX_FILE)
breakpoints_matrix = get_matrix(BREAKPOINTS_MATRIX_FILE)
rate_matrix = compute_breakpoint_rate(breakpoints_matrix, GENOME_SIZE_MB)

plot_matrices(j_matrix, rate_matrix, SPECIES, OUT_PLOT)
