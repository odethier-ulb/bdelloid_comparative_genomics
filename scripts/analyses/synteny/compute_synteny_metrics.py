import os
import re
import numpy as np
import pandas as pd
from itertools import combinations
from collections import defaultdict

"""
Compute synteny metrics between species pairs based on colinear regions retrieved with MCScanX.
"""

# Constants
BLOCK_RE = re.compile(
    r'^## Alignment\s+\d+:.*?N=(\d+)\s+(\S+)&(\S+)\s+(plus|minus)'
)
GENE_RE = re.compile(r'^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)')

SPECIES = ['av', 'ar', 'as', 'rr', 'mq', 'pr', 'hs']
HOMOEOLOGOUS_PAIRS = set([(1, 4), (2, 5), (3, 6)])

MCSCANX_GFF = 'rotifers_mcscanx.gff'
MCSCANX_COLINEARITY = 'rotifers_mcscanx.collinearity'
OUT_DIR = './/bdelloids/synteny_metrics'
OUT_FILE = f'{OUT_DIR}/bdelloid_synteny_metrics.tsv'

# ---------------------------- PARSERS -----------------------------#

def load_gff(gff_file, sp1, sp2):
    """
    Load MCScanX GFF input and export dictionaries of gene pos and index
    per species/chromosome.
    """
    df = pd.read_csv(gff_file, sep='\t', header=None,
                     names=['chr', 'gene', 'start', 'end'])
    gene_order = {sp1: defaultdict(list), sp2: defaultdict(list)}
    for row in df.itertuples(index=False):
        sp = row.chr[:2]
        if sp in gene_order:
            gene_order[sp][row.chr].append((row.gene, row.start))
    for sp in gene_order:
        for chrom in gene_order[sp]:
            gene_order[sp][chrom].sort(key=lambda x: x[1])
    gene_index = {
        sp: {chrom: {g: i for i, (g, _) in enumerate(genes)}
             for chrom, genes in gene_order[sp].items()}
        for sp in gene_order
    }
    return gene_order, gene_index


def chr_num(chrom):
    """Return the chromosome number from its name."""
    return int(re.search(r'(\d+)$', chrom).group(1))


def parse_pair_collinearity(colinearity_file, sp1, sp2, homoeolog_pairs=None):
    """
    Parse the MCScanX colinearity output file to retrieve
    colinear blocks and their corresponding anchored genes
    for species sp1 and sp2.

    Skip results between homoelog pairs (if provided).
    """
    homoeolog_pairs = homoeolog_pairs or set()
    blocks = []
    current = None

    with open(colinearity_file) as fh:
        for line in fh:
            m = BLOCK_RE.match(line)
            if m:
                _, chrA, chrB, strand = m.groups()
                spA, spB = chrA[:2], chrB[:2]

                if spA == spB or {spA, spB} != {sp1, sp2}:
                    current = None
                    continue

                pair_nums = tuple(sorted((chr_num(chrA), chr_num(chrB))))
                if pair_nums in homoeolog_pairs:
                    current = None
                    continue

                swap = (spA != sp1)
                chrom_sp1 = chrB if swap else chrA
                chrom_sp2 = chrA if swap else chrB

                current = {
                    'chrom_sp1': chrom_sp1, 'chrom_sp2': chrom_sp2,
                    'strand': strand, 'swap': swap,
                    'genes_sp1': [], 'genes_sp2': [],
                }
                blocks.append(current)
            else:
                m = GENE_RE.match(line)
                if m and current is not None:
                    geneA, geneB = m.groups()
                    if current['swap']:
                        current['genes_sp1'].append(geneB)
                        current['genes_sp2'].append(geneA)
                    else:
                        current['genes_sp1'].append(geneA)
                        current['genes_sp2'].append(geneB)

    return blocks


# -------------------------- COMPUTE METRICS -----------------------#

def compute_synteny_metrics(gff_file, collinearity_file, sp1, sp2, total_genes,
                             homoeolog_pairs=None):
    """
    For a pair of species, compute:
    - Jaccard index: number genes in colinear blocks over all the genes of both species.
    - Breakpoints: number of gap between regions of genes in colinear blocks.
    - Transpositions: number of intra-chromosomal rearrangements.
    - Non-reciprocal translocations: number of transfers of colinear blocks between 
        non-homologous chromosomes. 
    - Inversions: number of colinear blocks in reverse order.

    Skip results between homoelog pairs (if provided).
    """
    homoeolog_pairs = homoeolog_pairs or set()
    sp_first, _ = sorted([sp1, sp2])
    gene_order, gene_index = load_gff(gff_file, sp1, sp2)

    # Init a structure to keep track of which genes are in 
    # syntenic blocks
    synteny_covered = {
        sp: {chrom: [False] * len(genes) for chrom, genes in gene_order[sp].items()}
        for sp in (sp1, sp2)
    }
    # Init a structure to keep track of syntenic blocks to check their
    # orders.
    synteny_blocks = {sp: defaultdict(list) for sp in (sp1, sp2)}
    # Init synteny metrics
    n_inversions = 0
    n_nr_translocations = 0 # nr = non-reciprocal
    blocks = parse_pair_collinearity(collinearity_file, sp1, sp2, homoeolog_pairs)

    # Order blocks
    def block_sp_first_positions(block):
        """Get block's index positions."""
        chrom = block['chrom_sp1'] if sp_first == sp1 else block['chrom_sp2']
        genes = block['genes_sp1'] if sp_first == sp1 else block['genes_sp2']
        return chrom, [gene_index[sp_first][chrom][g]
                       for g in genes if g in gene_index[sp_first].get(chrom, {})]
    keyed_blocks = []
    for block in blocks:
        chrom, idxs = block_sp_first_positions(block)
        keyed_blocks.append((chrom, min(idxs), block))
    # Sort by sp_first's chromosome and pos
    keyed_blocks.sort(key=lambda t: (t[0], t[1]))
    # Assign block index
    idx_counter = defaultdict(int)
    for chrom, _, block in keyed_blocks:
        idx_counter[chrom] += 1
        block['global_idx'] = idx_counter[chrom]
    blocks = [b for _, _, b in keyed_blocks]

    # Loop over blocks
    for block in blocks:
        chrom1, chrom2 = block['chrom_sp1'], block['chrom_sp2']
        orthologous = chr_num(chrom1) == chr_num(chrom2)

        # Update coverage for each gene covered by the pos
        # of the first and last genes in a block
        for sp, chrom, genes in [(sp1, chrom1, block['genes_sp1']),
                                 (sp2, chrom2, block['genes_sp2'])]:
            idxs = [gene_index[sp][chrom][g] for g in genes
                   if g in gene_index[sp].get(chrom, {})]
            if not idxs:
                continue
            lo, hi = min(idxs), max(idxs)
            for i in range(lo, hi + 1):
                synteny_covered[sp][chrom][i] = True

        # Count non-reciprocal translocations.
        # Homoeologous blocks never reach this point: they were already
        # excluded in parse_pair_collinearity, so any non-orthologous
        # block remaining here is a genuine translocation.
        if not orthologous:
            n_nr_translocations += 1
            continue   # excluded from inversion/transposition detection
        # Count inversions
        if block['strand'] == 'minus':
            n_inversions += 1
        # Save current block index
        for sp, chrom in [(sp1, chrom1), (sp2, chrom2)]:
            genes = block['genes_sp1'] if sp == sp1 else block['genes_sp2']
            idxs = [gene_index[sp][chrom][g] for g in genes
                   if g in gene_index[sp].get(chrom, {})]
            if idxs:
                synteny_blocks[sp][chrom].append((min(idxs), block['global_idx']))

    # Count the transpositions
    n_transpositions = 0
    for sp in (sp1, sp2):
        for chrom, entries in synteny_blocks[sp].items():
            entries.sort(key=lambda x: x[0])
            for i in range(1, len(entries)):
                if entries[i][1] < entries[i - 1][1]:
                    n_transpositions += 1

    # Count the breakpoints
    n_breakpoints = 0
    for sp in (sp1, sp2):
        for chrom, covered in synteny_covered[sp].items():
            covered_idxs = [i for i, c in enumerate(covered) if c]
            if len(covered_idxs) < 2:
                continue
            first, last = covered_idxs[0], covered_idxs[-1]
            in_gap = False
            for i in range(first, last + 1):
                if not covered[i]:
                    if not in_gap:
                        n_breakpoints += 1
                        in_gap = True
                else:
                    in_gap = False

    # Compute a Jaccard index
    n_covered_sp1 = sum(sum(c) for c in synteny_covered[sp1].values())
    n_covered_sp2 = sum(sum(c) for c in synteny_covered[sp2].values())
    jaccard = (n_covered_sp1 + n_covered_sp2) / (total_genes[sp1] + total_genes[sp2])

    return {
        'species_1': sp1,
        'species_2': sp2,
        'jaccard': jaccard,
        'n_breakpoints': n_breakpoints,
        'n_transpositions': n_transpositions,
        'n_inversions': n_inversions,
        'n_nr_translocations': n_nr_translocations,
    }


#----------------------------- EXPORT -------------------------------#

def build_matrix(results_df, metric, species):
    """ 
    Export the results as a symmetric matrix.
    """
    matrix = pd.DataFrame(np.nan, index=species, columns=species)
    for _, row in results_df.iterrows():
        sp1, sp2 = row['species_1'], row['species_2']
        matrix.loc[sp1, sp2] = row[metric]
        matrix.loc[sp2, sp1] = row[metric]
    return matrix


def export_all_matrices(results_tsv, species, out_dir, out_prefix):
    """
    Read a results TSV (as produced by the main loop below) and export
    one symmetric matrix TSV per metric, ready for heatmap plotting.
    """
    results_df = pd.read_csv(results_tsv, sep='\t')
    metrics = ['jaccard', 'n_breakpoints', 'n_transpositions',
              'n_inversions', 'n_nr_translocations']
    for metric in metrics:
        matrix = build_matrix(results_df, metric, species)
        out_file = os.path.join(out_dir, f'{out_prefix}_{metric}_matrix.tsv')
        matrix.to_csv(out_file, sep='\t')
        
        
#-------------------------------- RUN-------------------------------#

if __name__ == '__main__':
    # Retrieve total number of genes per species 
    total_genes = {}
    with open(MCSCANX_GFF) as f:
        for line in f:
            chrom = line.split('\t')[0]
            sp = chrom[:2]
            if sp in SPECIES:
                total_genes[sp] = total_genes.get(sp, 0) + 1
    # Compute results
    results = []
    for sp1, sp2 in combinations(sorted(SPECIES), 2):
        print(f'Processing {sp1}-{sp2}...')
        res = compute_synteny_metrics(MCSCANX_GFF, MCSCANX_COLINEARITY, sp1, sp2, total_genes, HOMOEOLOGOUS_PAIRS)
        results.append(res)

    # Export results
    df = pd.DataFrame(results)
    df.to_csv(OUT_FILE, sep='\t', index=False)
    export_all_matrices(OUT_FILE, SPECIES, OUT_DIR, 'bdelloid')
