# Usage: compute_pdist.py <region_sizes_file> <tree_file> <cnvs_file>

import numpy as np
from scicone import Tree
import sys

region_sizes_file = sys.argv[1]
tree_file = sys.argv[2]
cell_node_ids_file = sys.argv[3]

tree_name = tree_file.split('_')[0]

region_sizes = np.loadtxt(region_sizes_file, delimiter=',').ravel()
n_regions = region_sizes.shape[0]
n_bins = np.sum(region_sizes)

cell_node_ids = np.loadtxt(cell_node_ids_file, delimiter='\t')

tree = Tree('', '')
tree.outputs['region_sizes'] =  region_sizes
tree.outputs['region_neutral_states'] = np.ones((n_regions,)) * 2
tree.outputs['cell_node_ids'] = cell_node_ids
tree.read_from_file(tree_file)
#tree.create_cell_node_ids()
tree.count_nodes_bins()

#n_regions = np.sum(np.any(np.abs(np.diff(cnvs, axis=1)) > 0, axis=0))
#nodes = np.sort(np.array(list(tree.node_dict.keys())).astype(int))
#for id1 in nodes:
#    for id2 in nodes:
#        d = tree.get_distance(str(int(id1)), str(int(id2)), distance='events')
#        s = n_bins/n_regions
#        print(id1, id2, d, s, d/s)

dist = tree.get_pairwise_cell_distances(distance='events')
dist *= n_regions/n_bins

np.savetxt(tree_name + '_pdist.txt', dist, delimiter=',')
