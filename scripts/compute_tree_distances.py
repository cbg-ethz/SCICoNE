# Usage: compute_tree_distances.py pdist1.txt,pdist2.txt,...,pdistN.txt [--max_dist 1]

import argparse 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('pdist_files')
parser.add_argument('--max_dist', default=1)
args = parser.parse_args()

pdist_files = args.pdist_files.split(',')
n_trees = len(pdist_files)
max_dist = args.max_dist

def tree_dist(pdist1, pdist2, max_dist=1):
    n_cells = pdist1.shape[0]

    # Fraction of cells with distance bigger than max_dist
    dists = np.abs(pdist1-pdist2)
    return np.sum(dists[-1] > max_dist) / n_cells

# Load cell-cell distance matrices
pcdists = []
for pdist_file in pdist_files:
    pcdists.append(np.loadtxt(pdist_file, delimiter=','))

# Compute pairwise tree distances between trees
ptdists = [] 
for i in range(n_trees):
    for j in range(i):
        ptdists.append(tree_dist(pcdists[i], pcdists[j], max_dist))


# Compute average 
print(np.mean(ptdists))
