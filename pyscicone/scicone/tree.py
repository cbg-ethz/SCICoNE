import os
import string
import subprocess, re
import numpy as np
import pandas as pd
from graphviz import Source
from scicone.constants import *

class Tree(object):
    def __init__(self, binary_path, output_path, postfix='PYSCICONETREETEMP', persistence=False, ploidy=2, copy_number_limit=6):
        self.binary_path = binary_path

        self.output_path = output_path
        self.persistence = persistence
        self.postfix = postfix if postfix != "" else "PYSCICONETREETEMP"

        self.traces = None
        self.best_scores = None
        self.graphviz_str = ""

        self.ploidy = ploidy
        self.copy_number_limit = copy_number_limit

        self.num_nodes = 0

        self.tree_str = ""
        self.posterior_score = 0.
        self.tree_prior_score = 0.
        self.event_prior_score = 0.
        self.log_likelihood = 0.
        self.root_score = 0.
        self.score = 0.
        self.nu = 1.
        self.node_dict = dict() #node_dict = dict(key=p_id, val=event_vec)

        self.outputs = dict()#tree_inferred=None, inferred_cnvs=None,
                            # rel_markov_chain=None, cell_node_ids=None,
                            # cell_region_cnvs=None, acceptance_ratio=None,
                            # gamma_values=None, attachment_scores=None, markov_chain=None)

        self.cell_node_labels = []

    def get_n_children_per_node(self):
        n_children = dict()
        for node_id in self.node_dict:
            l = len([self.node_dict[i]['parent_id'] for i in self.node_dict if self.node_dict[i]['parent_id'] == node_id])
            n_children[node_id] = l

        return n_children

    def get_avg_n_children_per_node(self, all=False):
        n_children = self.get_n_children_per_node()
        avg = np.mean([n_children[n] for n in n_children if (n_children[n] > 0 and not all) or all])
        return avg

    def get_node_depths(self):
        # Go to each node and follow up to the root

        def get_n_levels(node_id, n_levels=-1):
            try:
                n_levels = get_n_levels(self.node_dict[node_id]['parent_id'], n_levels+1)
            except:
                return n_levels
            return n_levels

        depths = dict()
        for node_id in self.node_dict:
            depths[node_id] = get_n_levels(node_id)

        return depths

    def get_max_depth(self):
        node_depths = self.get_node_depths()
        return max([node_depths[node] for node in node_depths])

    def get_node_malignancies(self):
        # Compute distance to root of each clone
        dists = dict()
        root = self.node_dict['0']['cnv']
        for node_id in self.node_dict:
            dists[node_id] = np.linalg.norm(self.node_dict[node_id]['cnv'] - root, ord=1)
        return dists

    def set_node_cnvs(self):
        # Set root state
        self.node_dict['0']['cnv'] = np.ones(self.outputs['inferred_cnvs'].shape[1],)
        bin_start = 0
        bin_end = 0
        for region, state in enumerate(self.outputs['region_neutral_states']):
            region_size = self.outputs['region_sizes'][region].astype(int)
            bin_end = bin_start + region_size
            self.node_dict['0']['cnv'][bin_start:bin_end] = state
            bin_start = bin_end

        for node_id in self.node_dict:
            if node_id != '0':
                parent_id = self.node_dict[node_id]['parent_id']
                bins_event = np.zeros(self.node_dict['0']['cnv'].size)
                for region in self.node_dict[node_id]['region_event_dict']:
                    region = int(region)
                    bin_start = np.sum(self.outputs['region_sizes'][:region].astype(int))
                    region_size = self.outputs['region_sizes'][region].astype(int)
                    bin_end = bin_start + region_size
                    bins_event[bin_start:bin_end] = int(self.node_dict[node_id]['region_event_dict'][str(region)])
                self.node_dict[node_id]['cnv'] = self.node_dict[parent_id]['cnv'] + bins_event

    def update_tree_str(self):
        tree_strs = []
        tree_strs.append(f"Tree posterior: {self.posterior_score}")
        tree_strs.append(f"Tree prior: {self.tree_prior_score}")
        tree_strs.append(f"Event prior: {self.event_prior_score}")
        tree_strs.append(f"Log likelihood: {self.log_likelihood}")
        tree_strs.append(f"Root score: {self.root_score}")
        tree_strs.append(f"Tree score: {self.score}")
        tree_strs.append(f"Nu: {self.nu}")

        for node in self.node_dict:
            event_str = ','.join([f'{key}:{val}' for key, val in self.node_dict[node]["region_event_dict"].items()])
            event_str = f'[{event_str}]'
            node_str = f'node {node}: p_id:{self.node_dict[node]["parent_id"]},{event_str}'
            tree_strs.append(node_str)

        tree_str = '\n'.join(tree_strs) + '\n'

        self.tree_str = tree_str

        return tree_str

    def read_tree_str(self, tree_str, num_labels=False):
        self.tree_str = tree_str

        list_tree_file = tree_str.split('\n')

        for line in list_tree_file:
            if line.startswith("Tree posterior:"):
                self.posterior_score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Tree prior:"):
                self.tree_prior_score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Event prior:"):
                self.event_prior_score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Log likelihood:"):
                self.log_likelihood = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Root score:"):
                self.root_score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Tree score:"):
                self.score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Nu:"):
                self.nu = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("node"):
                colon_splits = line.split(':', 2)
                node_id = colon_splits[0].split(" ")[1]
                parent_id, event_str = colon_splits[2].split(',', 1)
                event_str = event_str[1:-1] # Remove '[' and ']'
                if parent_id != 'NULL':
                    event_dict = dict([s.split(':') for s in event_str.split(',')])
                else:
                    event_dict = dict()
                self.node_dict[node_id] = dict(parent_id=parent_id, region_event_dict=event_dict)

        # Set each node's CNV profile
        self.set_node_cnvs()

        # # Sort the nodes by depth and distance to the root
        # node_ids = []
        # node_depths = self.get_node_depths()
        # node_mals = self.get_node_malignancies()
        # node_ids = pd.DataFrame({'depth':node_depths, 'mal':node_mals}).sort_values(
        #                         ['depth', 'mal'], ascending=[True, True]).index.tolist()


        # Label the non-empty nodes
        node_ids = np.array(list(self.node_dict.keys())).astype(int)
        node_ids = list(node_ids.astype(str))
        nodes, counts = np.unique(self.outputs['cell_node_ids'][:,-1], return_counts=True)
        node_sizes = dict(zip(nodes.astype(int).astype(str), counts))
        i = 0
        for node in node_ids:
            self.node_dict[node]['label'] = ""
            self.node_dict[node]['size'] = 0
            try: # only add label if node has cells attached
                if node_sizes[node] > 0:
                    self.node_dict[node]['size'] = int(node_sizes[node])
                    if num_labels:
                        self.node_dict[node]['label'] = str(i)
                    else:
                        self.node_dict[node]['label'] = list(string.ascii_uppercase)[i]
                    i += 1
            except KeyError:
                pass

        # Indicate cell node assignments via labels
        self.cell_node_labels = [self.node_dict[str(int(node))]['label'] for node in self.outputs['cell_node_ids'][:,-1]]

        self.num_nodes = len(list(self.node_dict.keys()))

    def read_from_file(self, file, num_labels=False):
        """
            reads the file containing a tree and converts it to graphviz format
            :param file: path to the tree file.
        """
        with open(file) as f:
            list_tree_file = list(f)

        self.tree_str = ''.join(list_tree_file)
        self.read_tree_str(self.tree_str, num_labels=num_labels)

    def learn_tree(self, segmented_data, segmented_region_sizes, n_iters=1000, move_probs=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.01],
                    n_nodes=3,  seed=42, postfix="", initial_tree=None, nu=1.0, cluster_sizes=None, region_neutral_states=None, alpha=0., max_scoring=True, copy_number_limit=2,
                    c_penalise=10.0, lambda_r=0.2, lambda_c=0.1, ploidy=2, verbosity=1, verbose=False, num_labels=False):
        if postfix == "":
            postfix = self.postfix

        n_cells, n_regions = segmented_data.shape
        move_probs_str = ",".join(str(p) for p in move_probs)

        # Save temporary files
        temp_segmented_data_file = f"{postfix}_tree_temp_segmented_data.txt"
        np.savetxt(temp_segmented_data_file, segmented_data, delimiter=',')

        temp_segmented_region_sizes_file = f"{postfix}_tree_temp_segmented_region_sizes.txt"
        np.savetxt(temp_segmented_region_sizes_file, segmented_region_sizes, delimiter=',')

        temp_cluster_sizes_file = ""
        if cluster_sizes is not None:
            temp_cluster_sizes_file = f"{postfix}_tree_temp_cluster_sizes.txt"
            np.savetxt(temp_cluster_sizes_file, cluster_sizes, delimiter=',')

        temp_region_neutral_states_file = ""
        if region_neutral_states is not None:
            temp_region_neutral_states_file = f"{postfix}_tree_temp_region_neutral_states_file.txt"
            np.savetxt(temp_region_neutral_states_file, region_neutral_states, delimiter=',')

        if initial_tree is not None:
            print("Using initial tree.")
            temp_tree_file = f"{postfix}_temp_tree.txt"
            f = open(temp_tree_file, "w")
            f.write(initial_tree.tree_str)
            f.close()
            nu = initial_tree.nu

            try:
                cmd = [self.binary_path, f"--d_matrix_file={temp_segmented_data_file}", f"--n_regions={n_regions}",\
                    f"--n_cells={n_cells}", f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={postfix}",\
                    f"--copy_number_limit={copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--tree_file={temp_tree_file}", f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}",\
                    f"--max_scoring={max_scoring}", f"--c_penalise={c_penalise}", f"--lambda_r={lambda_r}",
                    f"--lambda_c={lambda_c}", f"--region_neutral_states_file={temp_region_neutral_states_file}"]
                if verbose:
                    print(' '.join(cmd))
                cmd_output = subprocess.run(cmd)
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

            os.remove(temp_tree_file)
        else:
            try:
                cmd = [self.binary_path, f"--d_matrix_file={temp_segmented_data_file}", f"--n_regions={n_regions}",\
                    f"--n_cells={n_cells}", f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={postfix}",\
                    f"--copy_number_limit={copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}", f"--max_scoring={max_scoring}",\
                    f"--c_penalise={c_penalise}", f"--lambda_r={lambda_r}", f"--lambda_c={lambda_c}",
                    f"--region_neutral_states_file={temp_region_neutral_states_file}"]
                if verbose:
                    print(' '.join(cmd))
                cmd_output = subprocess.run(cmd)
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        self.outputs['region_sizes'] = segmented_region_sizes
        if region_neutral_states is not None:
            self.outputs['region_neutral_states'] = region_neutral_states
        else:
            self.outputs['region_neutral_states'] = np.ones((n_regions,)) * ploidy

        try:
            os.remove(temp_segmented_data_file)
        except Exception as e:
            print(f'Could not delete {temp_segmented_data_file}: {e}')
        try:
            os.remove(temp_segmented_region_sizes_file)
        except Exception as e:
            print(f'Could not delete {temp_segmented_data_file}: {e}')

        if cluster_sizes is not None:
            os.remove(temp_cluster_sizes_file)
        if region_neutral_states is not None:
            os.remove(temp_region_neutral_states_file)

        # Read in all the outputs
        try:
            cwd = os.getcwd()
            for fn in os.listdir(cwd):
                if postfix in fn and (".csv" in fn or ".tsv" in fn): # Read in all data
                    key = fn.split(postfix)[1].split('_', 1)[-1].split('.')[0]
                    delim = ','
                    if ".tsv" in fn:
                        delim = '\t'
                    self.outputs[key] = np.loadtxt(fn, delimiter=delim)
                    if len(self.outputs[key].shape) == 1: # protect against 1 cluster
                        self.outputs[key] = self.outputs[key].reshape(1, -1)

                    if not self.persistence:
                        os.remove(fn)
            # Read tree after everything else is read
            for fn in os.listdir(cwd):
                if postfix in fn and "tree_inferred.txt" in fn: # Parse tree structure, score and nu
                    self.read_from_file(fn, num_labels=num_labels)
                    if not self.persistence:
                        os.remove(fn)
        except Exception as e:
            print(f'Could not load {fn}: {e}')
            # print("OSError: ", e.output, e.stdout, e.stderr)

        return cmd_output

    def set_gene_event_dicts(self, region_gene_map):
        for node in self.node_dict:
            self.node_dict[node]['gene_event_dict'] = dict()
            if node != '0':
                for region in self.node_dict[node]['region_event_dict']:
                    for gene in region_gene_map[int(region)]:
                        if gene != "":
                            self.node_dict[node]['gene_event_dict'][gene] = self.node_dict[node]['region_event_dict'][region]

    def set_graphviz_str(self, root_label='Neutral', node_sizes=True, node_labels=True, color="#E6E6FA", event_fontsize=14, nodesize_fontsize=14,
                                nodelabel_fontsize=14, gene_labels=False, gene_list=None, tumor_type=None):
        if gene_labels is True and gene_list is None:
            # Load COSMIC gene list
            bpath = os.path.join(os.path.dirname(__file__), 'data')
            df = pd.read_csv(os.path.join(bpath, 'cancer_gene_census.csv'))
            gene_list = df['Gene Symbol'].tolist()
            if tumor_type is not None:
                gene_list = df.loc[np.where(np.array([tumor_type in val for val in df['Tumour Types(Somatic)'].values.astype(str)]))[0], 'Gene Symbol'].tolist()

        if node_sizes:
            nodes, counts = np.unique(self.outputs['cell_node_ids'][:,-1], return_counts=True)
            node_sizes = dict()
            for i, node in enumerate(nodes):
                node_sizes[str(int(node))] = counts[i]
        else:
            node_sizes = None

        graphviz_header = [
                    "digraph {",
                    f'node [style=filled,color="{color}",fontsize={event_fontsize},margin=0,shape=oval]'
                    f'edge [arrowhead=none, color="{color}"]',
                ]

        graphviz_labels = []
        graphviz_links = []

        for key in self.node_dict:
            node_id = key
            p_id = self.node_dict[key]['parent_id']
            if node_id == '0':
                str_merged_labels = root_label
            elif node_id == '-100':
                str_merged_labels = 'Whole-genome duplication'
            else:
                merged_labels = []

                event_dict = self.node_dict[key]['region_event_dict']
                first_region = list(event_dict.keys())[0]
                previous_event = event_dict[first_region]
                last_region = first_region
                for i, region in enumerate(event_dict):
                    if i > 0:
                        event = event_dict[region]
                        if int(region) == int(last_region) + 1 and event == previous_event:
                            last_region = region  # update the end
                        else:
                            if first_region == last_region:
                                merged_labels.append(f"{int(previous_event):+}R{first_region}")
                            else:
                                merged_labels.append(f"{int(previous_event):+}R{first_region}:{last_region}")
                            first_region = last_region = region
                        previous_event = event
                if first_region == last_region:
                    merged_labels.append(f"{int(previous_event):+}R{first_region}")
                else:
                    merged_labels.append(f"{int(previous_event):+}R{first_region}:{last_region}")

                # Add line breaks
                region_str_merged_labels = " ".join(
                        f"{x}<br/>" if i % 5 == 0 and i > 0 else str(x)
                        for i, x in enumerate(merged_labels)
                    )
                # Remove trailing line break
                if ''.join(list(region_str_merged_labels)[-len('<br/>'):]) == '<br/>':
                    region_str_merged_labels = ''.join(list(region_str_merged_labels)[:-len('<br/>')])

                if gene_labels:
                    try:
                        merged_labels = []
                        event_dict = self.node_dict[key]['gene_event_dict']

                        # Get all events
                        unique_events = np.unique(np.array(list(self.node_dict[key]['region_event_dict'].values())))

                        # Sort by descending
                        sorted_unique_events = np.sort(unique_events.astype(int))[::-1]

                        for event in sorted_unique_events:
                            if event > 0:
                                color = 'red'
                            else:
                                color = 'blue'

                            # Get all genes with this event
                            genes = []
                            for gene in event_dict:
                                if int(event_dict[gene]) == event:
                                    if gene_list is not None:
                                        if np.any(gene == np.array(gene_list)):
                                            genes.append(gene)
                                    else:
                                        genes.append(gene)

                            # Sort alphabetically
                            genes.sort()

                            # Add line breaks
                            merged_genes = ", ".join(
                                    f"{x}<br/>" if i % 5 == 0 and i > 0 else str(x)
                                    for i, x in enumerate(genes)
                                )
                            merged_genes = merged_genes.replace('<br/>,',',<br/>')

                            # Remove trailing line break
                            if ''.join(list(merged_genes)[-len('<br/>'):]) == '<br/>':
                                merged_genes = ''.join(list(merged_genes)[:-len('<br/>')])

                            event_str = ''
                            if len(genes) != 0:
                                event_str = f'<font point-size="{event_fontsize}" color="{color}">{event:+}</font>: ' + merged_genes + '<br/><br/>'
                                merged_labels.append(event_str)

                            str_merged_labels = ''.join(merged_labels)

                            # Remove trailing line breaks
                            if ''.join(list(str_merged_labels)[-len('<br/><br/>'):]) == '<br/><br/>':
                                str_merged_labels = ''.join(list(str_merged_labels)[:-len('<br/><br/>')])

                            # Count amplified and deleted regions
                            number_of_amps = np.sum(['+' in s for s in list(region_str_merged_labels)])
                            number_of_dels = np.sum(['-' in s for s in list(region_str_merged_labels)])
                            if len(merged_labels) == 0:
                                str_merged_labels = f'<font point-size="{event_fontsize}">({number_of_amps}+, {number_of_dels}-)</font>'
                            else:
                                str_merged_labels = f'<font point-size="{event_fontsize}">({number_of_amps}+, {number_of_dels}-)</font>' + '<br/><br/>' + str_merged_labels
                    except:
                        str_merged_labels = region_str_merged_labels
                else:
                    str_merged_labels = region_str_merged_labels

            # Add node label at the top
            if node_labels:
                if self.node_dict[node_id]['label'] != "":
                    try:
                        labcolor = LABEL_COLORS_DICT[self.node_dict[node_id]["label"]]
                    except KeyError:
                        labcolor = LABEL_COLORS_DICT_NUM[self.node_dict[node_id]["label"]]
                    node_label = (
                        f'<font point-size="{nodelabel_fontsize}" color="{labcolor}"><b>'
                        + str(self.node_dict[node_id]['label'])
                        + "</b></font>"
                    )
                    str_merged_labels = node_label + "<br/><br/>" + str_merged_labels

            # Add node size
            if node_sizes is not None:
                try:
                    node_size = node_sizes[node_id]
                except KeyError:
                    node_size = 0
                str_merged_labels = str_merged_labels + "<br/><br/>"
                str_merged_labels = (
                    str_merged_labels
                    + f'<font point-size="{nodesize_fontsize}">'
                    + str(int(node_size))
                    + " cell"
                )
                if int(node_size) > 1 or int(node_size) == 0:
                    str_merged_labels = str_merged_labels + "s"
                str_merged_labels = str_merged_labels + "</font>"


            graphviz_labels.append(
                f"{node_id}[label=<{str_merged_labels}>]"
            )  # use < > to allow HTML
            if p_id != 'NULL':
                graphviz_links.append(f"{p_id} -> {node_id}")

        self.graphviz_str = '\n'.join(graphviz_header + graphviz_labels + graphviz_links + ["}"])

    def plot_tree(self, root_label='Neutral', node_sizes=True, node_labels=True, color="#E6E6FA", event_fontsize=14, nodesize_fontsize=14, nodelabel_fontsize=14,
                    gene_labels=False, gene_list=None, tumor_type=None):

        self.set_graphviz_str(root_label=root_label, node_sizes=node_sizes, node_labels=node_labels, color=color,
                                event_fontsize=event_fontsize, nodesize_fontsize=nodesize_fontsize,
                                nodelabel_fontsize=nodelabel_fontsize,
                                gene_labels=gene_labels, gene_list=gene_list, tumor_type=tumor_type)

        s = Source(self.graphviz_str)
        return s


    def adjust_to_wgd(self, threshold=0.98):
        # if data is None and self.data is not None:
        #     data = self.data['filtered_counts']

        turned_diploid = False
        # Get the cells assigned to tetraploid root
        cells = np.where(self.outputs['cell_node_ids'] == '0')
        #
        # smaller_lib_size = True if np.median(np.sum(data[cells], axis=1)) < np.median(np.sum(data[~cells], axis=1)) else False
        # smaller_variance = True if np.median(np.var(data[cells], axis=1)) < np.median(np.var(data[~cells], axis=1)) else False
        tetraploid_nodes = np.unique(self.outputs['cell_node_ids'][cells]).astype(int).astype(str)
        # if smaller_lib_size and smaller_variance: # They are diploid
        if len(cells) > 0:
            self.outputs['inferred_cnvs'][cells] = (self.outputs['inferred_cnvs'][cells]/2).astype(int)
            # Re-assign to root
            self.outputs['cell_node_ids'][cells] = 0
            self.node_dict['0']['cnv'] = (self.node_dict['0']['cnv']/2).astype(int)
            turned_diploid = True

        # Add WGD
        if turned_diploid:
            # Remove the diploid nodes
            for c in tetraploid_nodes:
                if str(int(c)) != '0':
                    parent = self.node_dict[str(int(c))]['parent_id']

                    # Get children
                    for node_id in self.node_dict:
                        if self.node_dict[node_id]['parent_id'] == str(int(c)):
                            # Update parent
                            self.node_dict[node_id]['parent_id'] = parent

                    # Remove node
                    del self.node_dict[str(int(c))]

            # Add a WGD node below the root
            self.node_dict['-100'] = dict()
            self.node_dict['-100']['parent_id'] = '0'
            self.node_dict['-100']['label'] = ""
            self.node_dict['-100']['region_event_dict'] = dict()
            self.node_dict['-100']['size'] = 0
            for node_id in self.node_dict:
                if node_id != '-100' and node_id != '0':
                    if self.node_dict[node_id]['parent_id'] == '0':
                        self.node_dict[node_id]['parent_id'] = '-100'
