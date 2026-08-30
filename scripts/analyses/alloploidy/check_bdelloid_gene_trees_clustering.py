import logging
from Bio import Phylo
from collections import defaultdict
from scipy.stats import fisher_exact, combine_pvalues

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

"""
Parse a list of MUL trees built using pairs of homoeologous genes from bdelloid rotifers and 
single copy genes for other (non-bdelloid) species, and analyse whether genes cluster by 
their chromosome of origin (i.e. by homoeolog) across species.
"""


BDELLOIDS = {'avaga', 'aspwild', 'aricciae', 'rrotatoria', 'hspwild', 'proseola'}
HOMOEOLOG_PAIRS = [(1, 4), (3, 6), (2, 5)]
 
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
 
 
def _build_pair_lookup(pairs):
    """
    Build a dict mapping each chromosome to (pair_label, role), where role is 'a' or 'b'
    depending on its position in the HOMOEOLOG_PAIRS tuple.
    """
    lookup = {}
    for chr_a, chr_b in pairs:
        label = f'{chr_a}-{chr_b}'
        lookup[chr_a] = (label, 'a')
        lookup[chr_b] = (label, 'b')
    return lookup
 
 
PAIR_LOOKUP = _build_pair_lookup(HOMOEOLOG_PAIRS)
 
 
#------------------------- TOOLS -------------------------#
 
 
def all_parents(tree):
    """
    From https://biopython.org/wiki/Phylo_cookbook,
    generate a dictionary mapping all nodes to their parents.
    """
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents
 
 
def is_bdelloid(clade):
    """
    Return true if the clade is leaf node of a bdelloid gene
    """
    return clade.is_terminal and clade.name.rsplit('_')[-1] in BDELLOIDS
 
 
def get_bdelloid_clades(tree):
    """
    Return a list of all bdelloid clades found in the tree
    """
    return [c for c in tree.get_terminals() if is_bdelloid(c)]
 
 
def contains_all(root, clades):
    """
    Return true if all clade in clades can be found starting from the root.
    """
    return {c.name for c in clades} <= {c.name for c in root.get_terminals()}
 
 
def find_species_root(clade, species, parents):
    """
    Return the first clade containing all the given species.
    """
    if contains_all(clade, species):
        return clade
    else:
        return find_species_root(parents[clade], species, parents)
 
 
def find_balanced_tree(root, min_clade):
    """
    Starting from the root, find a balanced tree having at least min_clade terminal
    clades. In this case, a tree is considered unbalanced if one of the children
    of the root has only 1 terminal clade.
    """
    left, right = root.clades
    left_t, right_t = len(left.get_terminals()), len(right.get_terminals())
 
    if left_t > 1 and right_t > 1 and (left_t + right_t) >= min_clade:
        return root
    elif left_t == 1 and right_t >= min_clade:
        return find_balanced_tree(right, min_clade)
    elif right_t == 1 and left_t >= min_clade:
        return find_balanced_tree(left, min_clade)
    else:
        return None
 
 
def count_homoeolog_genes(clades, line):
    """
    Given a list of terminal clades, return how many belong to each chromosome found
    among bdelloid genes, as a dict {chromosome: count}. Genes on chromosomes that
    aren't part of any known homoeolog pair are warned about and skipped.
    """
    counts = defaultdict(int)
    for clade in clades:
        sp = clade.name.rsplit('_')[-1]
        if sp not in BDELLOIDS:
            LOGGER.warning(f'non-bdelloid species found in bdelloid tree {line}: {sp}')
            continue
        chr_ = int(clade.name.rsplit('|')[-1].rsplit('_')[0])
        if chr_ not in PAIR_LOOKUP:
            LOGGER.warning(f'chromosome {chr_} at line {line} is not part of any '
                            f'defined homoeolog pair, skipping gene')
            continue
        counts[chr_] += 1
    return counts
 
 
def resolve_pair(chrs, line):
    """
    Given the set of bdelloid chromosomes found on both sides of a balanced split,
    check that they correspond to exactly one known homoeolog pair (one chromosome
    on each side). Return (pair_label, chr_a, chr_b) or None if the tree can't be
    resolved to a single clean pair.
    """
    if len(chrs) != 2:
        LOGGER.warning(f'expected 2 chromosomes forming one homoeolog pair at line '
                        f'{line}, found {sorted(chrs)}, skipping tree')
        return None
 
    labels = {PAIR_LOOKUP[c][0] for c in chrs}
    if len(labels) != 1:
        LOGGER.warning(f'chromosomes {sorted(chrs)} at line {line} do not belong to '
                        f'the same homoeolog pair, skipping tree')
        return None
 
    pair_label = labels.pop()
    chr_a = next(c for c in chrs if PAIR_LOOKUP[c][1] == 'a')
    chr_b = next(c for c in chrs if PAIR_LOOKUP[c][1] == 'b')
    return pair_label, chr_a, chr_b
 
 
#-------------------------- MAIN -------------------------#
 
 
def check_homoeolog_clustering(trees_file, out_file, summary_file):
    """
    Check, per homoeolog pair, whether bdelloid genes cluster by chromosome of origin
    in the given MUL trees and runs a per-tree Fisher's exact test on
    the 2x2 table [[chr_a_left, chr_a_right], [chr_b_left, chr_b_right]].
    """
    trees, tot_trees, valid_trees = Phylo.parse(trees_file, 'newick'), 0, 0
    rows = []
    pvalues_by_pair = defaultdict(list)
 
    for tree in trees:
        tot_trees += 1
        parents = all_parents(tree)
        bdelloids = get_bdelloid_clades(tree)
        bd_root = find_species_root(parents[bdelloids[0]], bdelloids, parents)
        balanced_bd_root = find_balanced_tree(bd_root, len(bdelloids) - 2)
 
        if balanced_bd_root is None:
            LOGGER.warning(f'Unbalanced tree at line: {tot_trees}')
            continue
 
        left, right = balanced_bd_root.clades
        counts_l = count_homoeolog_genes(left.get_terminals(), tot_trees)
        counts_r = count_homoeolog_genes(right.get_terminals(), tot_trees)
        chrs = set(counts_l) | set(counts_r)
 
        resolved = resolve_pair(chrs, tot_trees)
        if resolved is None:
            continue
        pair_label, chr_a, chr_b = resolved
 
        a_l, b_l = counts_l.get(chr_a, 0), counts_l.get(chr_b, 0)
        a_r, b_r = counts_r.get(chr_a, 0), counts_r.get(chr_b, 0)
 
        valid_trees += 1
        odds_ratio, p_value = fisher_exact([[a_l, a_r], [b_l, b_r]])
        pvalues_by_pair[pair_label].append(p_value)
 
        rows.append({
            'tree_index': tot_trees, 'pair': pair_label, 'chr_a': chr_a, 'chr_b': chr_b,
            'chr_a_left': a_l, 'chr_b_left': b_l, 'chr_a_right': a_r, 'chr_b_right': b_r,
            'odds_ratio': odds_ratio, 'p_value': p_value,
        })
 
    with open(out_file, 'w') as out:
        out.write('tree_index\tpair\tchr_a\tchr_b\tchr_a_left\tchr_b_left\t'
                   'chr_a_right\tchr_b_right\todds_ratio\tp_value\n')
        for r in rows:
            out.write(f"{r['tree_index']}\t{r['pair']}\t{r['chr_a']}\t{r['chr_b']}\t"
                       f"{r['chr_a_left']}\t{r['chr_b_left']}\t{r['chr_a_right']}\t"
                       f"{r['chr_b_right']}\t{r['odds_ratio']}\t{r['p_value']}\n")
 
    LOGGER.info(f'Total trees: {tot_trees}, valid (balanced, single clean pair): '
                f'{valid_trees}, skipped: {tot_trees - valid_trees}')
 
    summary_rows = []
    for pair_label in sorted(pvalues_by_pair, key=lambda p: HOMOEOLOG_PAIRS.index(
            tuple(int(x) for x in p.split('-')))):
        pvals = pvalues_by_pair[pair_label]
        n_sig = sum(1 for p in pvals if p < 0.05)
        combined_stat, combined_p = combine_pvalues(pvals, method='fisher')
        LOGGER.info(f"Pair {pair_label}: {len(pvals)} trees, {n_sig} individually "
                    f'significant ({100 * n_sig / len(pvals):.1f}%), '
                    f"Fisher's combined statistic={combined_stat:.2f}, "
                    f'p-value={combined_p:.3e}')
        summary_rows.append({
            'pair': pair_label, 'n_trees': len(pvals), 'n_significant': n_sig,
            'pct_significant': 100 * n_sig / len(pvals),
            'combined_statistic': combined_stat, 'combined_p_value': combined_p,
        })
 
    with open(summary_file, 'w') as out:
        out.write('pair\tn_trees\tn_significant\tpct_significant\t'
                   'combined_statistic\tcombined_p_value\n')
        for r in summary_rows:
            out.write(f"{r['pair']}\t{r['n_trees']}\t{r['n_significant']}\t"
                       f"{r['pct_significant']:.1f}\t{r['combined_statistic']:.2f}\t"
                       f"{r['combined_p_value']:.3e}\n")
 
    return rows, summary_rows


        
#-------------------------- RUN --------------------------#

trees_file = 'gene_trees_rooted.nwk'
out_file = 'bdelloid_mul_trees_clustering_analyses.tsv'
summary_file = 'bdelloid_homoeolog_pair_summary.tsv'
 
check_homoeolog_clustering(trees_file, out_file, summary_file)
