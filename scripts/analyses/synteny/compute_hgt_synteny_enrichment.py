# DISCLAIMER:
# This script was initially generated with assistance from Claude (Anthropic)
# and subsequently reviewed, modified, and adapted by the authors.
# The authors take full responsibility for the final code, its accuracy,
# implementation, and any resulting analyses.


"""
Compute HTGs enrichment in syntenic regions using a Fisher's exact test.

Outputs
-------
- TSV : per-species contingency tables and Fisher's exact test results
- Plot: grouped bar chart of % syntenic for HTGs vs self genes
        + odds ratio plot
"""

import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


plt.rcParams['svg.fonttype'] = 'none'


# ------------------------------------------------------------------ #
#                         CONSTANTS
# ------------------------------------------------------------------ #

SUBGENOME_1 = {'1', '2', '3'}
SUBGENOME_2 = {'4', '5', '6'}

BLOCK_RE = re.compile(r'^## Alignment\s+\d+:.*?N=(\d+)\s+(\S+)&(\S+)')
GENE_RE  = re.compile(r'^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)')


#------------------------------------------------------------------ #
#                              TOOLS
# ------------------------------------------------------------------ #

def chrom_number(chrom: str) -> str:
    return chrom[2:]


def same_subgenome(chrA: str, chrB: str) -> bool:
    numA, numB = chrom_number(chrA), chrom_number(chrB)
    return (numA in SUBGENOME_1) == (numB in SUBGENOME_1)


def load_htgs(species: dict, hgt_results_pattern: str, hgt_class_min=3) -> dict:
    htgs = {}
    for sp, shortname in species.items():
        df = pd.read_csv(hgt_results_pattern.format(shortname), sep='\t', header=0)
        htgs[sp] = set(df.loc[df['hgt_class'] >= hgt_class_min, 'prot_id'])
    return htgs


def species_from_chrom(chrom: str) -> str:
    return chrom[:2]


# ------------------------------------------------------------------ #
#                         PARSER
# ------------------------------------------------------------------ #

def build_syntenic_gene_sets(collinearity_file: str,
                              keep_species: set) -> dict:
    """
    For each species in keep_species, return the set of genes that
    appear in at least one inter-species collinear block with any
    other species (union across all pairwise comparisons).
    Cross-subgenome ohnologous blocks are excluded.

    Returns
    -------
    syntenic_genes : dict {species_code: set_of_gene_ids}
    """
    syntenic_genes = defaultdict(set)

    current_mode = None   # 'keep' | None
    current_spA  = None
    current_spB  = None

    with open(collinearity_file) as fh:
        for line in fh:
            m = BLOCK_RE.match(line)
            if m:
                chrA = m.group(2)
                chrB = m.group(3)
                spA  = species_from_chrom(chrA)
                spB  = species_from_chrom(chrB)

                # only inter-species blocks between kept species
                # and only orthologous (same) subgenomes
                if (spA != spB
                        and spA in keep_species
                        and spB in keep_species
                        and same_subgenome(chrA, chrB)):
                    current_mode = 'keep'
                    current_spA  = spA
                    current_spB  = spB
                else:
                    current_mode = None
                continue

            m = GENE_RE.match(line)
            if m and current_mode == 'keep':
                syntenic_genes[current_spA].add('-'.join(m.group(1).split('-')[1:]))
                syntenic_genes[current_spB].add('-'.join(m.group(2).split('-')[1:]))

    return dict(syntenic_genes)


# ------------------------------------------------------------------ #
#                    FISHER'S EXACT TEST
# ------------------------------------------------------------------ #

def run_enrichment(htgs: dict,
                   syntenic_genes: dict,
                   total_genes: dict) -> pd.DataFrame:
    """
    For each species, build the 2×2 table and run Fisher's exact test.

    Parameters
    ----------
    htgs          : {sp: set_of_hgt_gene_ids}
    syntenic_genes: {sp: set_of_syntenic_gene_ids}
    total_genes   : {sp: int}  total gene count

    Returns
    -------
    df : pd.DataFrame with columns:
         species, n_total, n_hgt, n_syntenic,
         hgt_syntenic, hgt_nonsyntenic,
         nonhgt_syntenic, nonhgt_nonsyntenic,
         pct_syntenic_hgt, pct_syntenic_nonhgt,
         odds_ratio, p_value, direction
    """
    rows = []
    for sp in sorted(htgs.keys()):
        hgt_set      = htgs[sp]
        syn_set      = syntenic_genes.get(sp, set())
        n_total      = total_genes[sp]
        n_hgt        = len(hgt_set)
        n_nonhgt     = n_total - n_hgt
        
        # 2×2 contingency table
        a = len(hgt_set    & syn_set)          # HGT     & syntenic
        b = len(hgt_set    - syn_set)          # HGT     & non-syntenic
        c = len(syn_set    - hgt_set)          # non-HGT & syntenic
        d = n_total - n_hgt - c                # non-HGT & non-syntenic

        table = [[a, b], [c, d]] 
        or_, p = fisher_exact(table, alternative='two-sided')

        rows.append({
            'species'             : sp,
            'n_total'             : n_total,
            'n_hgt'               : n_hgt,
            'n_syntenic'          : len(syn_set),
            'hgt_syntenic'        : a,
            'hgt_nonsyntenic'     : b,
            'nonhgt_syntenic'     : c,
            'nonhgt_nonsyntenic'  : d,
            'pct_syntenic_hgt'    : round(a / n_hgt    * 100, 2) if n_hgt    else np.nan,
            'pct_syntenic_nonhgt' : round(c / n_nonhgt * 100, 2) if n_nonhgt else np.nan,
            'odds_ratio'          : round(or_, 4),
            'p_value'             : round(p,   6),
            'direction'           : 'enriched' if or_ > 1 else 'depleted',
        })

    df = pd.DataFrame(rows)

    # Benjamini-Hochberg FDR correction across species
    reject, p_adj, _, _ = multipletests(df['p_value'], method='fdr_bh')
    df['p_adj_BH']    = p_adj.round(6)
    df['significant'] = reject

    return df


# ------------------------------------------------------------------ #
#                           PLOTS
# ------------------------------------------------------------------ #

SPECIES_NAMES = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}

def plot_results(df: pd.DataFrame, species_labels: dict, out_file: str):
    """
    Two-panel figure:
    Top   : grouped bar chart — % syntenic for HGT vs non-HGT genes
    Bottom: forest plot of odds ratios with 95% CI
    """
    import matplotlib.gridspec as gridspec
    from scipy.stats import norm

    order = [sp for sp in species_labels if sp in df['species'].values]
    df    = df.set_index('species').loc[order].reset_index()
    labels = [species_labels.get(sp, sp) for sp in order]
    x      = np.arange(len(order))
    width  = 0.32

    fig = plt.figure(figsize=(max(8, len(order) * 1.5), 9))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[2, 1.4], hspace=0.45)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # ---- top: % syntenic bar chart ---- #
    hgt_color    = '#c00000'
    nonhgt_color = '#2f5597'

    bars_hgt = ax1.bar(x - width / 2, df['pct_syntenic_hgt'],    width,
                       color=hgt_color,    alpha=0.85, edgecolor='white',
                       label='HTGs')
    bars_non = ax1.bar(x + width / 2, df['pct_syntenic_nonhgt'], width,
                       color=nonhgt_color, alpha=0.85, edgecolor='white',
                       label='Self genes')

    for bars in [bars_hgt, bars_non]:
        for bar in bars:
            h = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                     f'{h:.0f}%', ha='center', va='bottom',
                     fontsize=7, color='#333333')

    # significance stars above HGT bar
    for i, row in df.iterrows():
        if row['significant']:
            y_max = max(row['pct_syntenic_hgt'], row['pct_syntenic_nonhgt']) + 3
            star  = '***' if row['p_adj_BH'] < 0.001 else \
                    '**'  if row['p_adj_BH'] < 0.01  else '*'
            ax1.text(x[i], y_max, star,
                     ha='center', va='bottom', fontsize=9, color='#333333')

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=40, ha='right',
                         fontstyle='italic', fontsize=10)
    ax1.set_ylabel('% of genes in syntenic regions', fontsize=10)
    ax1.set_title('HGT enrichment in syntenic regions', fontsize=12)
    ax1.legend(frameon=True, facecolor='#f2f2f2', edgecolor='none', fontsize=10, 
               loc='upper right', bbox_to_anchor=(1.0, 1.03))
    ax1.spines[['top', 'right']].set_visible(False)

    # ---- bottom: odds ratio forest plot ---- #
    # approximate 95% CI from the 2×2 table using log-OR ± 1.96 * SE
    log_or = np.log(df['odds_ratio'].replace(0, np.nan))
    se_log = np.sqrt(
        1 / df['hgt_syntenic'].replace(0, np.nan) +
        1 / df['hgt_nonsyntenic'].replace(0, np.nan) +
        1 / df['nonhgt_syntenic'].replace(0, np.nan) +
        1 / df['nonhgt_nonsyntenic'].replace(0, np.nan)
    )
    ci_lo = np.exp(log_or - 1.96 * se_log)
    ci_hi = np.exp(log_or + 1.96 * se_log)

    colors = [hgt_color if sig else '#aaaaaa'
              for sig in df['significant']]

    ax2.axvline(1.0, color='#888888', linewidth=0.9,
                linestyle='--', zorder=0)

    for i, (_, or_, lo, hi, col) in enumerate(zip(
            order, df['odds_ratio'], ci_lo, ci_hi, colors)):
        ax2.plot([lo, hi], [i, i], color=col, linewidth=1.5, zorder=1)
        ax2.scatter(or_, i, color=col, s=55, zorder=2,
                    edgecolors='white', linewidths=0.4)

    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels(labels, fontstyle='italic', fontsize=10)
    ax2.set_xlabel('Odds ratio (syntenic enrichment of HTGs)', fontsize=10)
    ax2.set_title("Fisher's exact test odds ratio", fontsize=11)
    ax2.spines[['top', 'right']].set_visible(False)

    plt.savefig(out_file, bbox_inches='tight')
    print(f'Saved → {out_file}')


# ------------------------------------------------------------------ #
#                           EXECUTION
# ------------------------------------------------------------------ #

collinearity_file = 'rotifers_mcscanx.collinearity'
gff_file = 'mcscanx/default/rotifers_mcscanx.gff'
hgt_results = '{}_hgt_results_nr.tsv'
out_stats = 'hgt_synteny_enrichment.tsv'
out_plot = 'hgt_synteny_enrichment.svg'

species_labels = {
    'av': 'A. vaga',
    'as': 'Adineta sp.',
    'ar': 'A. ricciae',
    'rr': 'R. rotatoria',
    'mq': 'M. quadricornifera',
    'hs': 'Habrotrocha sp.',
    'pr': 'P. roseola',
}
species_shortnames = {
    'av': 'a_vaga',
    'as': 'a_sp_wild',
    'ar': 'a_ricciae',
    'rr': 'r_rotatoria',
    'mq': 'm_quadricornifera',
    'hs': 'h_sp_wild',
    'pr': 'p_roseola',
}
keep_species = set(species_labels.keys())

# gene totals from GFF
total_genes = defaultdict(int)
with open(gff_file) as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 3:
            continue
        sp = species_from_chrom(cols[0])
        if sp in keep_species:
            total_genes[sp] += 1

species = list(keep_species)
print('Loading HTGs...')
htgs = load_htgs(species_shortnames, hgt_results)

print('Building syntenic gene sets (union across all pairwise comparisons)...')
syntenic_genes = build_syntenic_gene_sets(collinearity_file, keep_species)
for sp in sorted(keep_species):
    print(f'  {sp}: {len(syntenic_genes.get(sp, set())):,} syntenic genes '
          f'/ {total_genes[sp]:,} total')

print('\nRunning Fisher\'s exact tests...')
results = run_enrichment(htgs, syntenic_genes, total_genes)
print(results[['species', 'pct_syntenic_hgt', 'pct_syntenic_nonhgt',
               'odds_ratio', 'p_value', 'p_adj_BH', 'significant']].to_string(index=False))

results.to_csv(out_stats, sep='\t', index=False)
print(f'\nSaved → {out_stats}')

plot_results(results, species_labels, out_plot)
