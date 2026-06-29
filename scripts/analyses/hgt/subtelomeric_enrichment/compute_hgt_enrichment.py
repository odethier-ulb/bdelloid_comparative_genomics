import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from scipy.stats import fisher_exact

plt.rcParams['svg.fonttype'] = 'none'
warnings.filterwarnings('ignore', category=RuntimeWarning)


"""
Check if horizontaly transferred genes are enriched in subtelomeric regions using
Fisher's exact test. The subtelomeric regions are defined by a split ratio of the chromosome length or gene positions.
"""

# Enrichment constants
SPLIT_RATIOS = [0.1, 0.15, 0.2, 0.25]
CHR_NAME_PATTERN = r'(?:(?<=\D)|(?<=^))(\d+)$'
SPECIES = [
    'a_vaga',
    'a_sp_wild',
    'a_ricciae',
    'r_rotatoria',
    'm_quadricornifera',
    'h_sp_wild',
    'p_roseola',
]

SPECIES_LABELS = {
    'a_vaga'           : 'A. vaga',
    'a_sp_wild'        : 'Adineta sp.',
    'a_ricciae'        : 'A. ricciae',
    'r_rotatoria'      : 'R. rotatoria',
    'm_quadricornifera': 'M. quadricornifera',
    'h_sp_wild'        : 'Habrotrocha sp.',
    'p_roseola'        : 'P. roseola',
}
# Plot constants
WINDOW_SIZE = 100_000
Y_MAX = 35
COLORS = {
    'subtelomere_pos': '#c00000',
    'subtelomere_idx': "#5a4cef",
    'core'       : '#2f5597',
    'density'    : '#333333',
    'hgt'        : '#c00000',
}


#-----------------------------------------------------------------------------#

def get_chr_size(fai):
    """
    Get chromosome sizes from fai file.
    
    The chromosome is identified by its number.
    """
    chr_sizes = {}
    with open(fai, 'rt') as f:
        for line in f:
            lsplt = line.strip().split('\t')
            if len(lsplt) > 0:
                match = re.search(CHR_NAME_PATTERN, lsplt[0])
                chr_sizes[int(match.group(1))] = int(lsplt[1]) 
    return chr_sizes            


def get_htgs(hgt_results):
    """
    Get list of horizontally transferred genes from hgt results file.
    """
    df = pd.read_csv(hgt_results, header=0, sep='\t')
    return set(df.loc[df['hgt_class'] >= 3, 'prot_id'])


def get_gene_pos(chrom_file, htgs):
    """
    Retrieve gene positions.
    """
    df = pd.read_csv(chrom_file, header=None, sep='\t', names=['gene_id', 'chrom', 'orientation', 'start', 'end'])
    df['position'] = (df['start'] + df['end']) / 2
    df['is_hgt'] = df['gene_id'].isin(htgs)
    return df.sort_values(by=['chrom', 'position'])


def compute_fisher_exact_test(subtelHTGs, subtelSelfs, coreHTGs, coreSelfs):
    """
    Compute Fisher's exact test for enrichment of HTGs in subtelomeric regions.
    """
    table = [[subtelHTGs, subtelSelfs],
             [coreHTGs, coreSelfs]]
    return fisher_exact(table, alternative='greater')
    

def compute_sp_hgt_enrichment(species, chrom_file, fai, hgt_results, split_ratios):
    """
    Compute HGT enrichment in subtelomeric regions
    """
    chr_sizes = get_chr_size(fai)
    gene_pos = get_gene_pos(chrom_file, get_htgs(hgt_results))
    results = []
    
    for chrom in gene_pos['chrom'].unique():
        for split_ratio in split_ratios:
            chr_size = chr_sizes[int(chrom[-1])]
            chrom_genes = gene_pos[gene_pos['chrom'] == chrom]
            
            # retrieve chromosome 'core' location
            core_pos_start, core_pos_end = split_ratio * chr_size, (1 - split_ratio) * chr_size
            core_idx_start, core_idx_end = split_ratio * len(chrom_genes), (1 - split_ratio) * len(chrom_genes)
            
            # compute stats
            # -- position based --#
            subL_pos = chrom_genes[chrom_genes['position'] < core_pos_start]
            subR_pos = chrom_genes[chrom_genes['position'] > core_pos_end]
            core_pos = chrom_genes[(chrom_genes['position'] > core_pos_start) & (chrom_genes['position'] < core_pos_end)]
            or_pos, p_pos = compute_fisher_exact_test(len(subL_pos[subL_pos['is_hgt'] == True]) + len(subR_pos[subR_pos['is_hgt'] == True]),
                                    len(subL_pos[subL_pos['is_hgt'] == False]) + len(subR_pos[subR_pos['is_hgt'] == False]),
                                    len(core_pos[core_pos['is_hgt'] == True]),
                                    len(core_pos[core_pos['is_hgt'] == False]))
            # -- index based --#
            subL_idx = chrom_genes.iloc[0:int(core_idx_start)]
            subR_idx = chrom_genes.iloc[int(core_idx_end):]
            core_idx = chrom_genes.iloc[int(core_idx_start):int(core_idx_end)] 
            or_idx, p_idx = compute_fisher_exact_test(len(subL_idx[subL_idx['is_hgt'] == True]) + len(subR_idx[subR_idx['is_hgt'] == True]),
                                    len(subL_idx[subL_idx['is_hgt'] == False]) + len(subR_idx[subR_idx['is_hgt'] == False]),
                                    len(core_idx[core_idx['is_hgt'] == True]),
                                    len(core_idx[core_idx['is_hgt'] == False]))
            # append results
            results.append({
                'species': species,
                'chr': int(chrom[-1]),
                'split_ratio': split_ratio,
                'p_pos': round(p_pos, 4),
                'or_pos': round(or_pos, 4),
                'p_idx': round(p_idx, 4),
                'or_idx': round(or_idx, 4),
            })          
    return results


def compute_enrichment_summary(df, out_file):
    """
    Summarize the enrichment results by counting for each chromosome how many time it has found 
    a significant (p-value < 0.05) enrichment of HTGs in subtelomeric regions.
    """
    summary = df.groupby(['chr', 'split_ratio']).agg(
        n_chr_p_pos = ('p_pos', lambda x: (x < 0.05).sum()),
        n_chr_p_idx = ('p_idx', lambda x: (x < 0.05).sum()),
    ).reset_index()
    summary.to_csv(out_file, sep='\t', index=False)
    
    
def plot_hgt_density(species, chrom_file, fai, htgs, split_ratio, out_file):
    """
    Plot HGT density along the chromosomes with an highlight of the subtelomeric regions.
    """
    # retrieve data prior to plot
    genes = get_gene_pos(chrom_file, htgs)
    chromosomes = sorted(genes['chrom'].unique())
    chr_sizes = get_chr_size(fai)
    n_chr = len(chromosomes)
    label = SPECIES_LABELS.get(species, species)
    
    # init plot
    n_cols = min(n_chr, 3)
    n_rows = int(np.ceil(n_chr / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 3.0, n_rows * 2),
        sharey=True, squeeze=False
    )
    large_genome = any(size > 1e7 for size in chr_sizes.values())
    tick_factor = 4e6 if large_genome else 2e6

    for idx, chrom in enumerate(chromosomes):
        row = idx // n_cols
        col = idx  % n_cols
        ax  = axes[row][col]

        chr_genes = genes[genes['chrom'] == chrom].reset_index(drop=True)
        n = len(chr_genes)
        midpoints = chr_genes['position'].values
        
        # --- get values to plot ---
        chr_size = chr_sizes[int(chrom[-1])]
        x_windows = [(i, min(chr_size, i + WINDOW_SIZE)) for i in range(0, chr_size + 1, WINDOW_SIZE + 1)]
        density_x, density_y = [(x1 + x2) / 2 for x1, x2 in x_windows], []
        
        for x1, x2 in x_windows:
            window_genes = chr_genes[(chr_genes['position'] >= x1) & (chr_genes['position'] < x2)]
            if len(window_genes) > 0:
                density_y.append(window_genes['is_hgt'].mean())
            else:
                density_y.append(0)
        
        # --- subtelomeric boundaries ---
        core_pos_start, core_pos_end = split_ratio * chr_size, (1 - split_ratio) * chr_size
        core_idx_start, core_idx_end = int(split_ratio * len(chr_genes)), int((1 - split_ratio) * len(chr_genes))
        core_idx_start_pos, core_idx_end_pos = midpoints[core_idx_start], midpoints[core_idx_end]
        
        # --- shade subtelomeric regions ---
        ax.axvspan(0, core_pos_start,  alpha=0.12, color=COLORS['subtelomere_pos'])
        ax.axvspan(core_pos_end, chr_size, alpha=0.12, color=COLORS['subtelomere_pos'])
        ax.axvspan(0, core_idx_start_pos,  alpha=0.12, color=COLORS['subtelomere_idx'])
        ax.axvspan(core_idx_end_pos, chr_size, alpha=0.12, color=COLORS['subtelomere_idx'])
        
        # --- boundary lines ---
        ax.axvline(core_pos_start,  color=COLORS['subtelomere_pos'], linewidth=0.9, linestyle='--', alpha=0.7)
        ax.axvline(core_pos_end, color=COLORS['subtelomere_pos'], linewidth=0.9, linestyle='--', alpha=0.7)
        ax.axvline(core_idx_start_pos,  color=COLORS['subtelomere_idx'], linewidth=0.9, linestyle='--', alpha=0.7)
        ax.axvline(core_idx_end_pos, color=COLORS['subtelomere_idx'], linewidth=0.9, linestyle='--', alpha=0.7)
        
        # --- density line ---
        if density_x:
            ax.plot(density_x, density_y,
                    color=COLORS['density'], linewidth=1.0, alpha=0.85)
            ax.fill_between(density_x, density_y,
                            color=COLORS['hgt'], alpha=0.25)

        ax.set_title(chrom, fontsize=8, fontstyle='italic')
        ax.set_xlim(0, chr_size)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_factor))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x / 1e6)}"))
        ax.set_xlabel('Position (Mb)', fontsize=7)
    
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=7)
        
        # y-label only on leftmost column
        if col == 0:
            ax.set_ylabel('#HTGs / #genes', fontsize=7)

    # hide unused axes in the last row
    for idx in range(n_chr, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    fig.suptitle(
        f'{label} - HTG density along chromosomes\n'
        f'(window size: {int(WINDOW_SIZE/1000)} kbp; subtelomeric regions split-ratio: {split_ratio})',
        fontsize=10, y=1.02
    )

    subtelo_pos_patch = mpatches.Patch(
        color=COLORS['subtelomere_pos'], alpha=0.3, label='Subtelomere: position-based'
    )
    subtelo_idx_patch = mpatches.Patch(
        color=COLORS['subtelomere_idx'], alpha=0.3, label='Subtelomere: gene index-based'
    )
    fig.legend(handles=[subtelo_pos_patch, subtelo_idx_patch], loc='upper right',
               fontsize=8, frameon=True, facecolor='#f2f2f2',
               edgecolor='none', bbox_to_anchor=(1.0, 1.0))

    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()
        
#-------------------------- RUN EXAMPLE ------------------------------#      
        
all_results = []        
results_output_file = './example/enrichment_results.tsv'
summary_output_file = './example/enrichment_summary.tsv'     
        
for sp in ['a_ricciae']:
    fai = f'./example/{sp}.assembly.fa.fai'
    chrom_file = f'./example/{sp}.chrom'
    hgt_results = f'./example/{sp}_hgt_results_nr.tsv'
    plot_file = f'./example/{sp}_hgt_density.svg'
    all_results.extend(compute_sp_hgt_enrichment(sp, chrom_file, fai, hgt_results, SPLIT_RATIOS))
    plot_hgt_density(sp, chrom_file, fai, get_htgs(hgt_results), 0.2, plot_file)

df = pd.DataFrame(all_results)
df.to_csv(results_output_file, sep='\t', index=False)
compute_enrichment_summary(df, summary_output_file)
