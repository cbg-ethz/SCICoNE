import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
from pathlib import Path
import os

from scicone import Tree

def set_events(tree):
        tree['event'] = np.zeros((tree['cnv'].shape[0],))
        # Fill events
        def descend(root):
                for child in root['children']:
                        child['event'] = child['cnv'] - root['cnv']
                        descend(child)
                return
        descend(tree)
        return tree

def get_n_undoing(tree):
        n_undoing = 0. 
        def descend(root):
                for child in root['children']:
                        par = root['cnv']
                        chi = child['cnv']
                        dif = chi-par
                        is_event = np.logical_or(np.logical_and(par > 2, dif < 0), np.logical_and(par < 2, dif > 0))
                        evs = np.where(is_event)[0]
                        child['n_undoing'] += len(evs)
                        print(child['label'], dif[is_event][:10], par[is_event][:10])
                        descend(child)
                return
        descend(tree)
        node_dict = tree_to_dict(tree)
        for node in node_dict:
                n_undoing += node_dict[node]['n_undoing']
        return n_undoing
        

def create_cell_node_ids(self):
        nodes = list(self.node_dict.keys())
        n_cells = self.outputs['inferred_cnvs'].shape[0]
        self.outputs['cell_node_ids'] = np.zeros((n_cells, 2))
        self.outputs['cell_node_ids'][:,0] = np.arange(n_cells)
        for node in nodes:
              idx = np.where(np.all(self.outputs['inferred_cnvs'] == self.node_dict[node]['cnv'], axis=1))[0]
              print(node, idx)
              self.outputs['cell_node_ids'][idx,1] = node
        print(np.unique(self.outputs['cell_node_ids'][:,1], return_counts=True))

def path_to_node(tree, node_id):
        path = []
        path.append(node_id)
        parent_id = tree[node_id]['parent_id']
        while parent_id != '-1' and parent_id != 'NULL':
                path.append(parent_id)
                parent_id = tree[parent_id]['parent_id']
        return path[::-1][:]

def path_between_nodes(tree, nodeA, nodeB):
        pathA = np.array(path_to_node(tree, nodeA))
        pathB = np.array(path_to_node(tree, nodeB))
        path = []
        # Get MRCA
        i = -1
        for node in pathA:
                if node in pathB:
                        i += 1
                else:
                        break
        mrca = pathA[i]
        pathA = np.array(pathA[::-1])
        # Get path from A to MRCA
        path = path + list(pathA[:np.where(pathA == mrca)[0][0]])
        # Get path from MRCA to B
        path = path + list(pathB[np.where(pathB == mrca)[0][0]:])
        return path

def get_distance(tree, id1, id2, n_nodes=False, paths=None):
        if paths is None:
                path = path_between_nodes(tree, id1, id2)
        else:
                try:
                        path = paths[id1][id2]
                except KeyError:
                        path = paths[id2][id1]

        dist = 0
        if n_nodes:
                dist = len(path)
        else:
                for i, node in enumerate(path):
                        if i <= len(path) - 2:
                                if tree[node]['cnv'][0] != -1 and tree[path[i+1]]['cnv'][0] != -1:
                                        dist += np.sum(np.abs(tree[node]['cnv'] - tree[path[i+1]]['cnv']))

        return dist

def get_pairwise_cell_distances(tree, n_cells, attachments, paths=None, n_nodes=False):
        mat = np.zeros((n_cells, n_cells))

        for i in range(1, n_cells):
                id1 = attachments['node'].loc[i]
                for j in range(i):
                        id2 = attachments['node'].loc[j]
                        mat[i][j] = get_distance(tree, str(id1), str(id2), paths=paths, n_nodes=n_nodes)

        return mat


def tree_to_dict(tree):
        tree_dict = dict()

        def descend(node, par_id):
                label = node['label']
                tree_dict[label] = dict()
                tree_dict[label]["cnv"] = node["cnv"]
                tree_dict[label]["parent_id"] = par_id
                tree_dict[label]["n_undoing"] = node["n_undoing"]
                for child in node['children']:
                        descend(child, node['label'])

        descend(tree, '-1')

        return tree_dict

def dict_to_tree(tree_dict, root_node):
        # Make children
        for i in tree_dict:
                tree_dict[i]["children"] = []
        for i in tree_dict:
                for j in tree_dict:
                        if tree_dict[j]["parent_id"] == i:
                                tree_dict[i]["children"].append(j)

        # Make tree
        root = {
            "cnv": tree_dict[root_node]["cnv"],
            "children": [],
            "label": root_node,
            "parent_id": '-1',
            "n_undoing": 0,
        }

        # Recursively construct tree
        def descend(super_tree, label):
            for i, child in enumerate(tree_dict[label]["children"]):
                super_tree["children"].append(
                    {
                        "cnv": tree_dict[child]["cnv"],
                        "children": [],
                        "label": child,
                        "parent_id": label,
                        "n_undoing": 0,
                    }
                )

                descend(super_tree["children"][-1], child)

        descend(root, root_node)
        return root

def get_conet_tree_distances(file_conet_tree, file_attachments, file_cnvs):
        tree = pd.read_csv(file_conet_tree, sep="-", header=None, names=['from', 'to'])
        attachments = pd.read_csv(file_attachments, sep=";", header=None, names=['cellname', 'idx', 'lstart', 'lend'])
        cnvs = np.loadtxt(file_cnvs, delimiter=',').T
        
        tree['from'] = tree['from'].apply(lambda x: x.replace('(','').replace(')','') )
        tree['to'] = tree['to'].apply(lambda x: x.replace('(','').replace(')','') )

        attachments['node'] = attachments[['lstart', 'lend']].apply(lambda x: ','.join(x[x.notnull()]), axis = 1)

        # Build the tree dict from files
        print("Building tree...")
        tree_dict = dict()
        for i, row in tree.iterrows():
                parent_id = row['from']
                child_id = row['to']
                # Get the first cell attached to this node
                cell_idx = np.where(attachments['node'].values == child_id)[0]
                if len(cell_idx) > 0:
                        cell_idx = cell_idx[0]
                        node_cnv = cnvs[cell_idx].ravel()
                else:
                        node_cnv = -1 * np.ones((cnvs.shape[1],))
                tree_dict[str(child_id)] = dict(parent_id=str(parent_id), cnv=node_cnv)

        # Add root
        tree_dict['0,0'] = dict(parent_id='NULL', cnv=2*np.ones(cnvs.shape[1],))

        # Fill nodes
        def fill(node_dict, root_node):
                # Create recursive dict
                tree = dict_to_tree(node_dict, root_node)

                # Fill empty with mean of children
                def descend(root):
                        cnvs = []
                        for child in root['children']:
                                cr = descend(child)
                                cnvs.extend([cr])
                        if root['cnv'][0] == -1:
                                print(root['label'])
                                cnvs = np.vstack(cnvs)
                                root['cnv'] = np.round(np.mean(cnvs, axis=0))
                                print(f"Filling this node: {root['label']}")
                                print(root['cnv'][:10])
                                print(cnvs[:,:10])
                        return root['cnv']
                descend(tree)

                # Update node_dict
                tree_dict = tree_to_dict(tree)
                return tree_dict

        tree_dict = fill(tree_dict, '0,0')

        # Compute all paths between nodes
        nodes = list(tree_dict.keys())
        reps = []
        for i, row in attachments.iterrows():
                if row['lstart'] == row['lend']:
                        reps.append(row['node'])

        reps = np.unique(reps)
        if len(reps):
                attachments = attachments.replace(reps,['0,0'])

        print("Computing distances...")
        # Compute cell-cell distances
        mat = get_pairwise_cell_distances(tree_dict, n_cells, attachments, n_nodes=False)
        return mat, tree_dict

def get_true_conet_tree_distances(file_conet_tree, file_attachments, file_cnvs):
        tree = pd.read_csv(file_conet_tree, sep=" ", header=None, names=['from', 'to', 'none'])
        tree = tree.drop(columns=['none'])
        attachments = pd.read_csv(file_attachments, sep=",", header=None, names=['lstart', 'lend'])
        cnvs = np.loadtxt(file_cnvs, delimiter=',')

        tree['from'] = tree['from'].apply(lambda x: x.replace('(','').replace(')','') )
        tree['to'] = tree['to'].apply(lambda x: x.replace('(','').replace(')','') )

        attachments['lstart'] = attachments['lstart'].astype(int).astype(str)
        attachments['lend'] = attachments['lend'].astype(int).astype(str)
        attachments['node'] = attachments[['lstart', 'lend']].apply(lambda x: ','.join(x[x.notnull()]), axis = 1)

        # Build the tree dict from files
        print("Building tree...")
        tree_dict = dict()
        for i, row in tree.iterrows():
                parent_id = row['from']
                child_id = row['to']
                # Get the first cell attached to this node
                cell_idx = np.where(attachments['node'].values == child_id)[0]
                if len(cell_idx) > 0:
                        cell_idx = cell_idx[0]
                        node_cnv = cnvs[cell_idx].ravel()
                else:
                        node_cnv = -1 * np.ones((cnvs.shape[1],))
                tree_dict[str(child_id)] = dict(parent_id=str(parent_id), cnv=node_cnv)

        # Add root
        tree_dict['0,0'] = dict(parent_id='NULL', cnv=2*np.ones(cnvs.shape[1],))

        # Fill nodes
        def fill(node_dict, root_node):
                # Create recursive dict
                tree = dict_to_tree(node_dict, root_node)

                # Fill empty with mean of children
                def descend(root):
                        cnvs = []
                        for child in root['children']:
                                cr = descend(child)
                                cnvs.extend([cr])
                        if root['cnv'][0] == -1:
                                print(root['label'])
                                cnvs = np.vstack(cnvs)
                                root['cnv'] = np.round(np.mean(cnvs, axis=0))
                                print(f"Filling this node: {root['label']}")
                                print(root['cnv'][:10])
                                print(cnvs[:,:10])
                        return root['cnv']
                descend(tree)
                
                # Update node_dict
                tree_dict = tree_to_dict(tree)
                return tree_dict

        tree_dict = fill(tree_dict, '0,0')
        print(tree_dict)

        # Compute all paths between nodes
        nodes = list(tree_dict.keys())
        reps = []
        for i, row in attachments.iterrows():
                if row['lstart'] == row['lend']:
                        reps.append(row['node'])

        reps = np.unique(reps)
        if len(reps):
                attachments = attachments.replace(reps,['0,0'])

        print("Computing true distances...")
        # Compute cell-cell distances
        mat = get_pairwise_cell_distances(tree_dict, n_cells, attachments, n_nodes=False)
        return mat, tree_dict

errs = 0
compute_errors = True
simsdir = "/cluster/work/bewi/members/pedrof/sc-dna/sims2020_20_nodes"
scicone_path = "/cluster/work/bewi/members/pedrof/tupro_code/SCICoNE/build3"
output_path = "/cluster/work/bewi/members/pedrof/"

all_folders = glob.glob(simsdir + "/results2023_boxplot3/*")
all_folders = np.sort([f for f in all_folders if "nu_on" not in f and "inference_full" not in f])

#methods = ['cluster_tree', 'full_tree', 'cluster_tree_sum', 'full_tree_sum', 'conet']
methods = [ 
'cluster_tree_unfair', 'full_tree_unfair', 'cluster_tree_sum_unfair', 'full_tree_sum_unfair',
'cluster_tree_fair', 'full_tree_fair', 'cluster_tree_sum_fair', 'full_tree_sum_fair',
'cluster_tree_unfair_adapted', 'full_tree_unfair_adapted', 'cluster_tree_sum_unfair_adapted', 'full_tree_sum_unfair_adapted',
'cluster_tree_fair_adapted', 'full_tree_fair_adapted', 'cluster_tree_sum_fair_adapted', 'full_tree_sum_fair_adapted', 'conet']
#methods = ['conet']

#methods = ['cluster_tree_fair', 'full_tree_fair', 'conet']

n_cells = 200
n_reps = 40
n_nodes = 20
n_bins = 1500
sim_folder = [folder for folder in all_folders if "simulation" in folder][0]

files_missing = []
n_files_per_method = n_reps
n_files_missing = {key: 0 for key in methods}

row_list = []
for rep_id in tqdm(range(n_reps)):
        print(f'{rep_id}:starting...')
        true_tree_path = sim_folder + f"/{n_nodes}nodes/" + str(rep_id) + "_tree.txt"
        true_cnvs_path = sim_folder + f"/{n_nodes}nodes/" + str(rep_id) + "_ground_truth.txt"
        true_attachments_path = sim_folder + f"/{n_nodes}nodes/" + str(rep_id) + "_attachments.txt"
        if compute_errors:
                true_cnvs = np.loadtxt(true_cnvs_path, delimiter=',')
                n_regions = np.count_nonzero(np.sum((np.diff(true_cnvs,axis=1) != 0) != 0, axis=0))
                
                # Load in true CONET tree
                real_distances, true_tree = get_true_conet_tree_distances(true_tree_path, true_attachments_path, true_cnvs_path)
                real_distances *= n_regions/n_bins                
                print("Real distances")
                print(real_distances[:10,:10])
                print(f'Norm: {np.linalg.norm(real_distances)}')
                nodes_tree = len(true_tree.keys())
                print(f"Number of nodes in real tree: {nodes_tree}")
                tree = dict_to_tree(true_tree, '0,0')
                n_undoing = get_n_undoing(tree)
                print("Undoing")
                print(f'Real tree total number of undoing events: {n_undoing}')

        for method in methods:
                method_ = str(method)
                if 'fair' in method_:
                        adapted = ''
                        if 'adapted' in method_:
                                adapted = '_adapted'
                                fair = method_.split('_')[-2]
                                method_ = method_.rsplit('_',2)[0]
                        else:
                                fair = method_.split('_')[-1]
                                method_ = method_.rsplit('_',1)[0]
                print(adapted, fair, method_)
                #if 'ginkgo' in method:
                #        method_ = 'medalt_tree_distances'
                method_folder = [folder for folder in all_folders if method_ in folder][0]
                if 'medalt' not in method and 'conet' not in method:
                        if "sum" in method_:
                                method_ = method_.replace('_sum', '')
                        inferred_tree_path = method_folder + f"/{fair}/{n_nodes}nodes/" + str(rep_id) + '_' + method_ + adapted + ".txt"
                        inferred_cnvs_path = method_folder + f"/{fair}/{n_nodes}nodes/" + str(rep_id) + '_' + method_ + f"_cnvs{adapted}.csv"
                        print(f"Looking for {inferred_cnvs_path}")
                        my_file = Path(inferred_cnvs_path)
                        if not my_file.is_file():
                                print(f'{inferred_cnvs_path} not found!')
                                n_files_missing[method] = n_files_missing[method] + 1
                                files_missing.append(inferred_cnvs_path)
                                continue
                elif 'medalt' in method:
                        distances_path = method_folder + f"/{n_nodes}nodes/" + str(rep_id) + '_' + method[:-1] + ".txt"
                        my_file = Path(distances_path)
                        if not my_file.is_file():
                                print(f'{distances_path} not found!')
                                n_files_missing[method] = n_files_missing[method] + 1
                                files_missing.append(inferred_cnvs_path)
                                continue
                elif 'conet' in method:
                        file_cnvs = method_folder + f"/{n_nodes}nodes/{rep_id}_conet_cbs_inferred.txt"
                        my_file = Path(file_cnvs)
                        if not my_file.is_file():
                                print(f'{file_cnvs} not found!')
                                n_files_missing[method] = n_files_missing[method] + 1
                                files_missing.append(file_cnvs)
                                continue
                
                if compute_errors:
                        if 'tree' in method:
                                # Load in tree
                                inf_regions_path = f"/cluster/work/bewi/members/pedrof/sc-dna/sims2020_20_nodes/results2023_boxplot3/bp_detection__trees_{fair}/" + f"{n_nodes}nodes/" + str(rep_id) + "_segmented_region_sizes.txt"
                                inf_regions = np.loadtxt(inf_regions_path, delimiter=',').ravel()

                                tree_inf = Tree(scicone_path, output_path, n_bins=n_bins)
                                tree_inf.outputs['inferred_cnvs'] = np.loadtxt(inferred_cnvs_path, delimiter=',')
                                tree_inf.outputs['region_sizes'] = inf_regions
                                tree_inf.outputs['region_neutral_states'] = np.ones((tree_inf.outputs['region_sizes'].shape[0],)) * 2
                                print('Inf: reading tree...')
                                tree_inf.read_from_file(inferred_tree_path)
                                print('Inf: creating cell_node_ids...')
                                create_cell_node_ids(tree_inf)
                                print('Inf: getting pairwise distances...')
                                tree_inf.count_nodes_bins()
                                inf_distances = tree_inf.get_pairwise_cell_distances(distance='events')
                                inf_distances *= n_regions/n_bins                                              
                                nodes_tree = len(tree_inf.node_dict.keys())
                                print(f"Number of nodes in SCICoNE tree: {nodes_tree}")
                                tree = dict_to_tree(tree_inf.node_dict, '0')
                                n_undoing = get_n_undoing(tree)
                                print("Undoing")
                                print(f'SCICoNE total number of undoing events: {n_undoing}')
                                #print(tree_inf.tree_str)
                                #print(f"n_regions: {len(tree_inf.outputs['region_sizes'])}")
                                #for node in tree_inf.node_dict:
                                #       print(tree_inf.node_dict[node]['region_event_dict'])
                                #       print(f"n_bins: {tree_inf.node_dict[node]['n_bins']}")

                                #node1 = list(tree_inf.node_dict.keys())[1]
                                #node2 = list(tree_inf.node_dict.keys())[-1]
                                #print(f"{node1}->{node2}: {tree_inf.path_between_nodes(node1, node2)}")
                        elif 'conet' in method:
                                 file_conet_tree = method_folder + f"/{n_nodes}nodes/{rep_id}_conet_cbs_tree.txt"
                                 file_attachments = method_folder + f"/{n_nodes}nodes/{rep_id}_conet_cbsattachments.txt"
                                 file_cnvs = method_folder + f"/{n_nodes}nodes/{rep_id}_conet_cbs_inferred.txt"
                                 try:
                                         inf_distances, conet_tree = get_conet_tree_distances(file_conet_tree, file_attachments, file_cnvs)
                                         inf_distances *= n_regions/n_bins 
                                         nodes_tree = len(conet_tree.keys())
                                         print(f"Number of nodes in CONET tree: {nodes_tree}") 
                                         tree = dict_to_tree(conet_tree, '0,0')
                                         n_undoing = get_n_undoing(tree)
                                         print("Undoing")
                                         print(f'CONET total number of undoing events: {n_undoing}')
                                 except Exception as e:
                                         print(e)
                                         print("CONET FAILED!") 
                                         errs += 1
                                         print(errs)
                                         continue                                 
                        else:
                                distances_path = method_folder + f"/{n_nodes}nodes/" + str(rep_id) + '_' + method[:-1] + ".txt"
                                inf_distances = np.loadtxt(distances_path, delimiter=',')
                                inf_distances *= n_regions/n_bins
                        

                        print(f"\n**{method}**")
                        print("Inferred distances")
                        print(inf_distances)
                        print(f'Norm: {np.linalg.norm(inf_distances)}') 
                        delta = np.linalg.norm(inf_distances - real_distances) / np.sqrt(n_cells*(n_cells-1)/2)
                        print(f'Tree delta: {delta}')
                        cnvs_delta = np.sqrt(np.mean((true_cnvs-tree_inf.outputs['inferred_cnvs'])**2))
                        print(f'CNVs delta: {cnvs_delta}')
                        tup = (rep_id, n_nodes, method, delta)
                        row_list.append(tup)
                        print(f'{rep_id}:{method} done!')

if len(files_missing) >= 0:
        print('No files missing.')
        if compute_errors:
                df = pd.DataFrame(data=row_list, columns=['rep_id', 'n_nodes', 'method', 'delta'])
                results_file = simsdir + f'/results2023_boxplot3/{n_nodes}_nodes_tau.csv'
                df.to_csv(results_file)
                print(df)
                print(f'Results file is available at {results_file}.')
else:
        print(f'FAILED: There are {len(files_missing)} files missing in total.')
        print(f'Number of missing CNV output files per method (out of {n_files_per_method}):')
        for method in n_files_missing:
                print(f'{method}: {n_files_missing[method]}')

