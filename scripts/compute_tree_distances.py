# Usage: compute_tree_distances.py --pdist_files pdist1.txt,pdist2.txt,...,pdistN.txt

import argparse 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('pdist_files')
args = parser.parse_args()

pdist_files = args.pdist_files.split(',')
n_trees = len(pdist_files)

def tree_dist(pdist1, pdist2):
    n_cells = pdist1.shape[0]
    return np.linalg.norm(pdist1 - pdist2) / np.sqrt(n_cells*(n_cells-1)/2)

# Load cell-cell distance matrices
pcdists = []
for pdist_file in pdist_files:
    pcdists.append(np.loadtxt(pdist_file, delimiter=','))

# Compute pairwise tree distances between trees
ptdists = [] 
for i in range(n_trees):
    for j in range(i):
        ptdists.append(tree_dist(pcdists[i], pcdists[j]))

# Compute average 
print(np.mean(ptdists))
