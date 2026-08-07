import logging
from Bio import Phylo
from scipy.stats import fisher_exact, combine_pvalues

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

"""
Parse a list of MUL trees built using pair of homoe/ohnolog genes from bdelloid rotifers and single
copy genes for other (non-bdelloid) species and analyse if bdelloid genes tend to cluster according
to their subgenome.
"""

BDELLOIDS = {'avaga', 'aspwild', 'aricciae', 'rrotatoria', 'hspwild', 'proseola'}
# Assign a subgenome to the chromomsomes based on the hypothesis that:
#   - chromosomes 1, 3, 5 belong to the subgenome 1
#   - chromosomes 2, 4, 6 belong to the subgenome 2
SUBGENOMES = {1: 'SG1', 2: 'SG2', 3: 'SG1', 4: 'SG2', 5: 'SG1', 6: 'SG2'}
# ...
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


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
    
    
def count_subgenome_genes(clades, line):
    """
    Given a list of terminal clades, return how many belong to the subgenomes 1 and 2
    """
    sg1, sg2, chrs = 0, 0, set()
    for clade in clades:
        sp = clade.name.rsplit('_')[-1]
        if not sp in BDELLOIDS:
            LOGGER.warning(f'non-bdelloid species found in bdelloid tree {line}: {sp}')
        else:
            chr = int(clade.name.rsplit('|')[-1].rsplit('_')[0])
            if SUBGENOMES[chr] == 'SG1':
                sg1 += 1
            else:
                sg2 += 1
            chrs.add(chr)
    return sg1, sg2, chrs 


        
#-------------------------- MAIN -------------------------#


def check_subgenome_clustering(trees_file, out_file):
    """
    Check if bdelloid genes cluster by subgenome in the given MUL trees.

    For each tree with a resolvable balanced split, runs a per-tree Fisher's
    exact test on the 2x2 table [[sg1_left, sg1_right], [sg2_left, sg2_right]].
    Results are written to out_file, and an aggregate p-value across all trees 
    is computed with Fisher's method for combining independent p-values.
    """
    trees, tot_trees, valid_trees = Phylo.parse(trees_file, 'newick'), 0, 0
    rows, p_values = [], []

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
        sg1_l, sg2_l, chrs_l = count_subgenome_genes(left.get_terminals(), tot_trees)
        sg1_r, sg2_r, chrs_r = count_subgenome_genes(right.get_terminals(), tot_trees)
        chrs = chrs_l.union(chrs_r)  
        
        valid_trees += 1

        odds_ratio, p_value = fisher_exact([[sg1_l, sg1_r], [sg2_l, sg2_r]])
        p_values.append(p_value)
        rows.append({
            'tree_index': tot_trees, 'chrs': ','.join(str(c) for c in sorted(chrs)),
            'sg1_left': sg1_l, 'sg2_left': sg2_l, 'sg1_right': sg1_r, 'sg2_right': sg2_r,
            'odds_ratio': odds_ratio, 'p_value': p_value,
        })

    with open(out_file, 'w') as out:
        out.write('tree_index\tchrs\tsg1_left\tsg2_left\tsg1_right\tsg2_right\todds_ratio\tp_value\n')
        for r in rows:
            out.write(f"{r['tree_index']}\t{r['chrs']}\t{r['sg1_left']}\t{r['sg2_left']}\t"
                       f"{r['sg1_right']}\t{r['sg2_right']}\t{r['odds_ratio']}\t{r['p_value']}\n")

    n_sig = sum(1 for p in p_values if p < 0.05)
    LOGGER.info(f'Total trees: {tot_trees}, valid (balanced): {valid_trees}, '
                f'unbalanced/skipped: {tot_trees - valid_trees}')
    LOGGER.info(f'Trees individually significant (p<0.05): {n_sig}/{valid_trees} '
                f'({100 * n_sig / valid_trees:.1f}%)')

    if p_values:
        combined_stat, combined_p = combine_pvalues(p_values, method='fisher')
        LOGGER.info(f"Fisher's combined test across all {len(p_values)} trees: "
                    f'statistic={combined_stat:.2f}, p-value={combined_p:.3e}')

    return rows, p_values

        
#-------------------------- RUN --------------------------#


trees_file = 'gene_trees_rooted.nwk'
out_file = 'bdelloid_mul_trees_clustering_analyses.tsv'

check_subgenome_clustering(trees_file, out_file)
