# Usage: compute_tree_distances.py pdist1.txt,pdist2.txt,...,pdistN.txt [--max_dist 1 --scores score1,score2,...scoreN]

import argparse 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('pdist_files')
parser.add_argument('--max_dist', default=1)
parser.add_argument('--scores', default="1,1", type=str)
args = parser.parse_args()

pdist_files = args.pdist_files.split(',')
max_dist = args.max_dist
scores = args.scores.split(',')

def tree_dist(pdist1, pdist2, max_dist=1):
    n_cells = pdist1.shape[0]

    # Fraction of cells with distance bigger than max_dist
    dists = np.abs(pdist1-pdist2)
    return np.sum(dists[-1] > max_dist) / n_cells

if len(scores) > 5:
    # Get top-5 scoring trees
    top_trees = np.argsort(np.array(scores).astype(float))[-5:]
else:
    top_trees = np.arange(len(pdist_files))

# Load cell-cell distance matrices
pcdists = []
for pdist_file in np.array(pdist_files)[top_trees]:
    pcdists.append(np.loadtxt(pdist_file, delimiter=','))

# Compute pairwise tree distances between trees
ptdists = [] 
for i in range(len(pcdists)):
    for j in range(i):
        ptdists.append(tree_dist(pcdists[i], pcdists[j], max_dist))


# Compute average 
print(np.mean(ptdists))
