
import dendropy
import itertools
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import mantel

"""
Perform a Mantel test with 10,000 permutations to check if
there is a correlation between the Jaccard index and the
patristic distance between species pairs.
"""

#------------------- CONSTANTS -----------------#

SPECIES_LABELS = {
    'avaga': 'av',
    'aspwild': 'as',
    'aricciae': 'ar',
    'rrotatoria': 'rr',
    'mquadricornifera': 'mq',
    'hspwild': 'hs',
    'proseola': 'pr',
}
TREE_FILE = 'bdelloid_tree.nwk'
JACCARD_FILE = 'bdelloid_jaccard_indices.tsv'


#--------------------- TOOLS -------------------#

def get_jaccard_matrix(jaccard_file):
    """Return a DataFrame of inter-species Jaccard indices."""
    df = pd.read_csv(jaccard_file, sep='\t', header=0, index_col=0)
    df = df.fillna(1.0)
    jaccard_df = 1.0 - df
    return jaccard_df


def get_patristic_matrix(treefile, label_map, species_order):
    """Return a DataFrame of pairwise patristic distances for focal taxa."""
    tree = dendropy.Tree.get(path=treefile, schema='newick')
    pdm  = tree.phylogenetic_distance_matrix()
    taxa = {t.label: t for t in tree.taxon_namespace}
    code2label = {v: k for k, v in label_map.items()}
    matrix = pd.DataFrame(0.0, index=species_order, columns=species_order)

    for cA, cB in itertools.combinations(species_order, 2):
        tA = taxa[code2label[cA]]
        tB = taxa[code2label[cB]]
        d  = pdm(tA, tB)
        matrix.loc[cA, cB] = d
        matrix.loc[cB, cA] = d

    return matrix


def compute_mantel_test(jaccard_df, patristic_df):
    """Compute correlation between distance matrices using the Mantel test."""
    jaccard_dist = DistanceMatrix(jaccard_df)
    patristic_dist = DistanceMatrix(patristic_df)
    coeff, p_value, n = mantel(jaccard_dist, patristic_dist, permutations=10000)
    print(f'Mantel test coeff = {coeff}; p_value = {p_value}; n = {n}')
    

#--------------------- RUN -------------------#

jaccard_df = get_jaccard_matrix(JACCARD_FILE)
patristic_df = get_patristic_matrix(TREE_FILE, SPECIES_LABELS,jaccard_df.index.to_list())
compute_mantel_test(jaccard_df, patristic_df)
