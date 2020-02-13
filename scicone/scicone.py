import os, shutil
import subprocess, re
from snakemake.workflow import Workflow, Rules
import snakemake.workflow
from snakemake import shell
from snakemake.logging import setup_logger
import numpy as np
import pandas as pd
import graphviz
from graphviz import Source

import multiprocessing
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist

import phenograph
from collections import Counter

import matplotlib.pyplot as plt
import seaborn as sns

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

        self.tree_str = ""
        self.posterior_score = 0.
        self.tree_prior_score = 0.
        self.event_prior_score = 0.
        self.log_likelihood = 0.
        self.root_score = 0.
        self.score = 0.
        self.nu = 1.

        self.outputs = dict()#tree_inferred=None, inferred_cnvs=None,
                            # rel_markov_chain=None, cell_node_ids=None,
                            # cell_region_cnvs=None, acceptance_ratio=None,
                            # gamma_values=None, attachment_scores=None, markov_chain=None)

    def read_from_file(self, file, output_path=None):
        """
            reads the file containing a tree and converts it to graphviz format
            :param file: path to the tree file.
            :param output_path: (optional) path to the output file.
        """
        with open(file) as f:
            list_tree_file = list(f)

        self.tree_str = ''.join(list_tree_file)

        graphviz_header = [
            "digraph {",
            'node [style=filled,color="#D4C0D6",fontsize=20,margin=0,shape=oval]'
            'edge [arrowhead=none, color="#602A86"]',
        ]

        graphviz_labels = []
        graphviz_links = []

        graphviz_labels.append("0[label=<<font point-size='30'> Neutral </font>>]")  # root

        for line in list_tree_file:
            if line.startswith("Tree posterior:"):
                self.posterior_score = float(line.split(' ')[-1].split('\n')[0])
            #
            # elif line.startswith("Tree prior:"):
            #     self.tree_prior_score = float(line.split(' ')[-1].split('\n')[0])
            #
            # elif line.startswith("Event prior:"):
            #     self.event_prior_score = float(line.split(' ')[-1].split('\n')[0])
            #
            # elif line.startswith("Log likelihood:"):
            #     self.log_likelihood = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Root score:"):
                self.root_score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Tree score:"):
                self.score = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("Nu:"):
                self.nu = float(line.split(' ')[-1].split('\n')[0])

            elif line.startswith("node 0:"):
                continue

            elif line.startswith("node"):
                comma_splits = line.split(",")

                comma_first = re.split(" |:", comma_splits[0])
                node_id = comma_first[1]
                p_id = comma_first[4]
                comma_rest = comma_splits[1:]
                comma_rest[0] = comma_rest[0].lstrip("[")
                comma_rest[-1] = comma_rest[-1].rstrip("]\n")
                merged_labels = []
                [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
                k_end = k_begin
                for term in comma_rest[1:]:  # events vector
                    [k, v] = (int(x) for x in term.split(":"))
                    if k == k_end + 1 and v == previous_v:
                        k_end = k  # update the end
                    else:
                        if k_begin == k_end:
                            merged_labels.append(f"{previous_v:+}R{k_begin}")
                        else:
                            merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
                        k_begin = k_end = k
                    previous_v = v
                # print the last one
                if k_begin == k_end:
                    merged_labels.append(f"{previous_v:+}R{k_begin}")
                else:
                    merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")

                str_merged_labels = (
                    "<i> </i>"
                    + " ".join(
                        f"{x}<i> </i><br/>" if i % 10 == 0 and i > 0 else str(x)
                        for i, x in enumerate(merged_labels)
                    )
                    + "<i> </i>"
                )

                newline = "<br/>"
                endline = "<i> </i>"

                # If there are newlines followed by endlines, switch
                str_merged_labels = str_merged_labels.replace(
                    newline + endline, endline + newline
                )

                # Remove whatever comes after the last endline position
                new_end_pos = (
                    [m.start() for m in re.finditer(endline, str_merged_labels)][-1]
                    + len(endline)
                    - 1
                )
                if len(str_merged_labels) > new_end_pos + 1:
                    str_merged_labels = str_merged_labels[: new_end_pos + 1]

                # Replace multiple endlines with one
                while endline * 2 in str_merged_labels:
                    str_merged_labels = str_merged_labels.replace(endline * 2, endline)

                graphviz_labels.append(
                    f"{node_id}[label=<{str_merged_labels}>]"
                )  # use < > to allow HTML
                graphviz_links.append(f"{p_id} -> {node_id}")

        self.graphviz_str = graphviz_header + graphviz_labels + graphviz_links + ["}"]

        if output_path is not None:
            with open(output_path, "w") as file:
                for line in self.graphviz_str:
                    file.write(f"{line}\n")

        self.graphviz_str = '\n'.join(self.graphviz_str)

    def learn_tree(self, segmented_data, segmented_region_sizes, n_iters=1000, move_probs=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.01],
                    n_nodes=3,  seed=42, postfix="", initial_tree=None, nu=1.0, cluster_sizes=None, region_neutral_states=None, alpha=0., max_scoring=True,copy_number_limit=2, verbosity=1):
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
            temp_tree_file = f"{postfix}_temp_tree.txt"
            f = open(temp_tree_file, "w")
            f.write(initial_tree.tree_str)

            try:
                cmd_output = subprocess.run([self.binary_path, f"--d_matrix_file={temp_segmented_data_file}", f"--n_regions={n_regions}",\
                    f"--n_cells={n_cells}", f"--ploidy={self.ploidy}", f"--verbosity={verbosity}", f"--postfix={postfix}",\
                    f"--copy_number_limit={copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--tree_file={temp_tree_file}", f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}",\
                    f"--max_scoring={max_scoring}", f"--region_neutral_states_file={temp_region_neutral_states_file}"])
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

            os.remove(temp_tree_file)
        else:
            try:
                cmd_output = subprocess.run([self.binary_path, f"--d_matrix_file={temp_segmented_data_file}", f"--n_regions={n_regions}",\
                    f"--n_cells={n_cells}", f"--ploidy={self.ploidy}", f"--verbosity={verbosity}", f"--postfix={postfix}",\
                    f"--copy_number_limit={copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}", f"--max_scoring={max_scoring}",\
                    f"--region_neutral_states_file={temp_region_neutral_states_file}"])
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.remove(temp_segmented_data_file)
        os.remove(temp_segmented_region_sizes_file)
        os.remove(temp_cluster_sizes_file)
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
                    if not self.persistence:
                        os.remove(fn)
                elif postfix in fn and "tree_inferred.txt" in fn: # Parse tree structure, score and nu
                    self.read_from_file(fn)
                    if not self.persistence:
                        os.remove(fn)
        except Exception:
            print('Could not load {}'.format(fn))
            # print("OSError: ", e.output, e.stdout, e.stderr)

        return cmd_output

    def plot_tree(self, show_genes=True):
        if self.graphviz_str != "":
            s = Source(self.graphviz_str)
            return s
        else:
            raise Exception("No graphviz string available. Please read from file or learn a new tree from data.")

class SCICoNE(object):
    """
    This class  provides an interface to interact with the outputs from the C++
    program that infers a copy number phylogeny from single-cell WGS data.

    It exposes functions to visualize breakpoint detection results, annotate
    the bins, and plot the resulting trees. Its attributes can be further
    processed using scgenpy, a generic package for pre, post-processing and
    visualization of single-cell copy number data.
    """
    def __init__(self, binary_path, output_path, persistence=False, postfix="PYSCICONETEMP", n_cells=0, n_bins=0):
        """Create a SCICoNE object.

        binary_path : type
            Path to SCICoNE binaries.
        output_path : type
            Path to SCICoNE output files.
        persistence : boolean
            Wether to delete output files from C++ after loading them into the class.
        """
        self.binary_path = binary_path
        self.simulation_binary = os.path.join(self.binary_path, 'simulation')
        self.bp_binary = os.path.join(self.binary_path, 'breakpoint_detection')
        self.inference_binary = os.path.join(self.binary_path, 'inference')
        self.score_binary = os.path.join(self.binary_path, 'score')

        self.output_path = output_path
        self.persistence = persistence
        self.postfix = postfix

        self.n_cells = n_cells
        self.n_bins = n_bins
        self.best_cluster_tree = None
        self.cluster_tree_robustness_score = 0.
        self.best_full_tree = None
        self.full_tree_robustness_score = 0.
        self.tree_list = []

    def simulate_data(self, n_cells=200, n_nodes=5, n_bins=1000, n_regions=40, n_reads=10000, nu=1.0, min_reg_size=10, max_regions_per_node=1, ploidy=2, region_neutral_states=None, verbosity=0):
        output_path = os.path.join(self.output_path, f"{self.postfix}_simulation")

        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        cwd = os.getcwd()

        region_neutral_states_file = ""
        if region_neutral_states is not None:
            # Use provided region_neutral_states instead of assuming the same state for all regions
            region_neutral_states_file = f"{self.postfix}_pre_region_neutral_states_file.txt"
            np.savetxt(region_neutral_states_file, region_neutral_states, delimiter=',')


        done = False
        while not done:
            try:
                cmd_output = subprocess.run([self.simulation_binary, f"--n_cells={n_cells}", f"--n_nodes={n_nodes}",\
                    f"--n_regions={n_regions}", f"--n_bins={n_bins}", f"--n_reads={n_reads}", f"--nu={nu}",\
                    f"--min_reg_size={min_reg_size}", f"--max_regions_per_node={max_regions_per_node}",\
                    f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={self.postfix}",\
                    f"--region_neutral_states_file={region_neutral_states_file}"])
            except subprocess.SubprocessError as e:
                print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            # The binary generates 4 files: *_d_mat.csv, *_ground_truth.csv, *_region_sizes.txt and *_tree.txt.
            # We read each one of these files into np.array, np.array, np.array and Tree structures, respectively,
            # and output a dictionary with keys "d_mat", "ground_truth", "region_sizes" and "tree".
            d_mat_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv"
            ground_truth_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv"
            region_sizes_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt"
            tree_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt"

            if region_neutral_states is not None:
                os.remove(region_neutral_states_file)

            try:
                # Move to output directory
                cwd = os.getcwd()
                os.rename(os.path.join(cwd, d_mat_file), os.path.join(output_path, d_mat_file))
                os.rename(os.path.join(cwd, ground_truth_file), os.path.join(output_path, ground_truth_file))
                os.rename(os.path.join(cwd, region_sizes_file), os.path.join(output_path, region_sizes_file))
                os.rename(os.path.join(cwd, tree_file), os.path.join(output_path, tree_file))
                done = True
            except OSError as e:
                pass

        d_mat_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv")
        ground_truth_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv")
        region_sizes_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt")
        tree_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt")

        try:
            output = dict()
            output['d_mat'] = np.loadtxt(d_mat_file, delimiter=',')
            output['ground_truth'] = np.loadtxt(ground_truth_file, delimiter=',')
            output['region_sizes'] = np.loadtxt(region_sizes_file, delimiter=',')

            tree = Tree(self.inference_binary, self.output_path)
            tree.read_from_file(tree_file)

            output['tree'] = tree

        except OSError as e:
            pass
            # print("OSError: ", e.output, e.stdout, e.stderr)

        # Now that we've read all outputs into memory, we delete the temporary files if persistence==False
        if not self.persistence:
            shutil.rmtree(output_path)

        return output

    def detect_breakpoints(self, data, window_size=30, threshold=3.0, bp_limit=300, lr=None, sp=None, evaluate_peaks=True, compute_lr=True, compute_sp=True, input_breakpoints=None):
        n_cells = data.shape[0]
        n_bins = data.shape[1]
        verbosity = 1 # > 0 to generate all files

        lr_file = ""
        if lr is not None:
            compute_lr = False
            # Use provided lr values instead of computing new ones
            if lr.shape[0] < lr.shape[1]:
                lr = lr.T # make sure it's bins by cells
            lr_file = f"{self.postfix}_pre_lr_vec.txt"
            np.savetxt(lr_file, lr, delimiter=',')

        sp_file = ""
        if sp is not None:
            compute_sp = False
            sp_file = f"{self.postfix}_pre_sp_vec.txt"
            np.savetxt(sp_file, sp, delimiter=',')

        input_breakpoints_file = ""
        if input_breakpoints is not None:
            input_breakpoints_file = f"{self.postfix}_pre_input_breakpoints_file.txt"
            np.savetxt(input_breakpoints_file, input_breakpoints, delimiter=',')

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{self.postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')

            cmd_output = subprocess.run([self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}",\
                f"--bp_limit={bp_limit}", f"--compute_lr={compute_lr}", f"--lr_file={lr_file}",\
                f"--compute_sp={compute_sp}", f"--sp_file={sp_file}", f"--verbosity={verbosity}",\
                f"--evaluate_peaks={evaluate_peaks}", f"--postfix={self.postfix}",\
                f"--input_breakpoints_file={input_breakpoints_file}"])

            # Delete the data file
            os.remove(data_file)

            if lr is not None:
                os.remove(lr_file)

            if sp is not None:
                os.remove(sp_file)
        except subprocess.SubprocessError as e:
            print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

        output = dict()

        try:
            cwd = os.getcwd()
            for fn in os.listdir(cwd):
                if self.postfix in fn:
                    key = fn.split(self.postfix)[1].split('_', 1)[1].split('.')[0]
                    output[key] = np.loadtxt(fn, delimiter=',')
                    if not self.persistence:
                        os.remove(fn)
        except OSError as e:
            print("OSError: ", e.output, e.stdout, e.stderr)

        return output

    def condense_regions(self, data, segmented_region_sizes):
        n_cells = data.shape[0]
        n_regions = len(segmented_region_sizes)
        sum_region_sizes = np.sum(segmented_region_sizes)
        condensed_mat = np.zeros((n_cells, n_regions))

        for i in range(n_cells):
            region_id = 0
            region_count = 0
            # import ipdb; ipdb.set_trace() # debugging starts here
            for j in range(data.shape[1]):
                to_add = data[i][j]
                condensed_mat[i][region_id] += to_add
                region_count += 1
                if region_count == segmented_region_sizes[region_id]:
                    region_id += 1
                    region_count = 0

        if not np.allclose(condensed_mat.sum(axis=1), data.sum(axis=1)):
            raise AssertionError(
                "not all values of the sums before & after "
                "segmentation are close")

        return condensed_mat

    def condense_segmented_clusters(self, segmented_data):
        # Cluster the segmented counts
        n_cells = segmented_data.shape[0]
        n_regions = segmented_data.shape[1]
        n_neighbours = int(n_cells / 10)
        communities, graph, Q = phenograph.cluster(data=segmented_data, k=n_neighbours, n_jobs=1, jaccard=True)
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]
        community_dict = dict((Counter(communities)))
        community_ids = sorted(list(community_dict))

        # Compute average counts of each cluster
        avg_segmented_counts = np.empty(segmented_data.shape)
        condensed_avg_segmented_counts = np.empty((len(community_ids), n_regions))
        cluster_sizes = np.zeros((len(community_ids),))
        for id in community_ids:
            avg_segmented_counts[np.where(communities==id)[0]] = np.mean(segmented_data[np.where(communities==id)[0], :], axis=0)
            condensed_avg_segmented_counts[id] = avg_segmented_counts[np.where(communities==id)[0][0],:]
            cluster_sizes[id] = np.where(communities==id)[0].shape[0]

        return condensed_avg_segmented_counts, cluster_sizes, communities

    def learn_single_tree(self, segmented_data, segmented_region_sizes, **tree_kwargs):
        tree = Tree(self.inference_binary, self.output_path)
        tree.learn_tree(segmented_data, segmented_region_sizes, **tree_kwargs)
        return tree

    def learn_tree_parallel(self, segmented_data, segmented_region_sizes, n_reps=10, **tree_kwargs):
        pool = Pool(n_reps)
        results = []
        for i in range(n_reps):
            kwargs = dict(seed=i+42, postfix=f"rep{i}_PYSCICONETREETEMP_{self.postfix}")
            kwargs.update(tree_kwargs)
            results.append(pool.apply_async(self.learn_single_tree, (segmented_data, segmented_region_sizes), kwargs))

        # Close the pool and wait for each running task to complete
        pool.close()
        pool.join()
        scores = []
        trees = []
        for result in results:
            tree = result.get()
            trees.append(tree)
            scores.append(tree.score)

        # Get best tree
        scores = np.array(scores)
        scores[scores==0.] = -np.inf # just in case
        best_tree = trees[np.argmax(scores)]

        # Assess convergence
        robustness_score = np.count_nonzero(np.isclose(np.array(scores), best_tree.score, atol=5)) / n_reps

        return best_tree, robustness_score, trees

    def learn_tree(self, segmented_data, segmented_region_sizes, n_reps=10, cluster=True, full=True, cluster_tree_n_iters=4000, full_tree_n_iters=4000, max_tries=2, robustness_thr=0.5, **kwargs):
        if "region_neutral_states" in kwargs:
            region_neutral_states = np.array(kwargs['region_neutral_states'])

            if np.any(region_neutral_states) < 0:
                raise Exception("Neutral states can not be negative!")

            # If there is a region with neutral state = 0, remove it to facilitate tree inference
            zero_neutral_regions = np.where(region_neutral_states==0)[0]
            if len(zero_neutral_regions) > 0:
                full_segmented_region_sizes = segmented_region_sizes.astype(int)
                full_segmented_data = segmented_data
                full_region_neutral_states = region_neutral_states
                segmented_data = np.delete(segmented_data, zero_neutral_regions, axis=1)
                segmented_region_sizes = np.delete(segmented_region_sizes, zero_neutral_regions)
                region_neutral_states = np.delete(region_neutral_states, zero_neutral_regions)

                kwargs['region_neutral_states'] = region_neutral_states

        if cluster:
            # Get the average read counts
            clustered_segmented_data, cluster_sizes, cluster_assignments = self.condense_segmented_clusters(segmented_data)

            cnt = 0
            robustness_score = 0.
            tree = None
            while robustness_score < robustness_thr:
                if cnt >= max_tries:
                    break
                nu = tree.nu if tree is not None else 1.0
                tree, robustness_score, trees = self.learn_tree_parallel(clustered_segmented_data, segmented_region_sizes, n_reps=n_reps, nu=nu, cluster_sizes=cluster_sizes, initial_tree=tree, n_iters=cluster_tree_n_iters, **kwargs)
                cnt += 1

            print(f"Cluster tree finished with a robustness score of {robustness_score} after {cnt} tries")

            # Expand clusters back into cells to get the cell per bin genotype
            cluster_bin_genotypes = tree.outputs['inferred_cnvs']
            cell_bin_genotypes = np.empty((segmented_data.shape[0], cluster_bin_genotypes.shape[1]))
            for id in np.unique(cluster_assignments):
                cell_bin_genotypes[np.where(cluster_assignments==id)[0]] = cluster_bin_genotypes[id]

            tree.outputs['inferred_cnvs'] = cell_bin_genotypes

            # Add back regions with neutral state = 0
            if "region_neutral_states" in kwargs:
                if len(zero_neutral_regions) > 0:
                    rec_cell_bin_genotypes = np.zeros((cell_bin_genotypes.shape[0], np.sum(full_segmented_region_sizes)))
                    new_start = 0
                    start = 0
                    for region in range(len(full_region_neutral_states)):
                        new_region_size = full_segmented_region_sizes[region]

                        new_end = new_start+new_region_size
                        if full_region_neutral_states[region] == 0:
                            rec_cell_bin_genotypes[:, new_start:new_end] = 0
                            new_start = new_end
                        else:
                            end = start+new_region_size
                            rec_cell_bin_genotypes[:, new_start:new_end] = cell_bin_genotypes[:, start:end]
                            start = end

                        new_start = new_end

                    tree.outputs['inferred_cnvs'] = rec_cell_bin_genotypes

            self.best_cluster_tree = tree
            self.cluster_tree_robustness_score = robustness_score

        if full:
            if self.best_cluster_tree is not None:
                tree = self.best_cluster_tree

                print('Initializing nu for full tree.')
                # Update the nu on the full data (the nu on the clustered data is very different) with this tree
                nu = tree.nu
                tree = self.learn_single_tree(segmented_data, segmented_region_sizes, nu=nu, initial_tree=tree, n_iters=5000, move_probs=[0,0,0,0,0,0,0,0,0,0,0,1,0], postfix=f"nu_tree_{self.postfix}", **kwargs)
                print('Done. Will start from nu={}'.format(tree.nu))

                cnt = 0
                robustness_score = 0.
                while robustness_score < robustness_thr:
                    if cnt >= max_tries:
                        break
                    nu = tree.nu
                    tree, robustness_score, trees = self.learn_tree_parallel(segmented_data, segmented_region_sizes, n_reps=n_reps, nu=nu, initial_tree=tree, n_iters=full_tree_n_iters)
                    cnt += 1

                cell_bin_genotypes = tree.outputs['inferred_cnvs']
                # Add back regions with neutral state = 0
                if "region_neutral_states" in kwargs:
                    if len(zero_neutral_regions) > 0:
                        rec_cell_bin_genotypes = np.zeros((cell_bin_genotypes.shape[0], np.sum(full_segmented_region_sizes)))
                        new_start = 0
                        start = 0
                        for region in range(len(full_region_neutral_states)):
                            new_region_size = full_segmented_region_sizes[region]

                            new_end = new_start+new_region_size
                            if full_region_neutral_states[region] == 0:
                                rec_cell_bin_genotypes[:, new_start:new_end] = 0
                                new_start = new_end
                            else:
                                end = start+new_region_size
                                rec_cell_bin_genotypes[:, new_start:new_end] = cell_bin_genotypes[:, start:end]
                                start = end

                            new_start = new_end

                        tree.outputs['inferred_cnvs'] = cell_bin_genotypes

                self.best_full_tree = tree
                self.full_tree_robustness_score = robustness_score
                self.tree_list = trees

                print(f"Full tree finished with a robustness score of {robustness_score} after {cnt} tries")
            else:
                raise Exception("Full trees require a cluster tree to start from. Please re-run with cluster=True.")

    def plot_bps(self, cluster_cells=True):
        pass

def filter_lr(lr_matrix, H=None):
    freq = np.fft.fftfreq(lr_matrix.shape[-1], 1) # 1 Hz sampling rate

    filtered_lr = np.empty(lr_matrix.shape)
    for c in range(lr_matrix.shape[0]):
        X = np.fft.fft(lr_matrix[c])
        Y = X * H
        y = np.fft.ifft(Y)
        filtered_lr[c] = y

    return filtered_lr


sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build_fix_root',
                  '/cluster/work/bewi/members/pedrof/', postfix="XYZA", persistence=False)

n_regions = 10
region_neutral_states = np.ones((n_regions,)) * 2
region_neutral_states[:2] = 1
region_neutral_states
data = sci.simulate_data(n_cells=200, n_nodes=1, n_bins=1000, n_regions=n_regions,
            n_reads=2000, nu=10, max_regions_per_node=10, min_reg_size=20, region_neutral_states=region_neutral_states)
true_tree = data['tree']
data.keys()
true_tree.plot_tree()


# Get true breakpoints
cell_genotypes = pd.DataFrame(data['ground_truth'])
cell_bps = cell_genotypes.diff(periods=1, axis=1)
cell_bps = cell_bps.fillna(value=0.0)
cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
grouped_cell_bps = cell_bps.sum(axis=0)
ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
ground_truth = ground_truth.index.tolist()

bps_indicator = np.zeros(grouped_cell_bps.shape[0])
bps_indicator[grouped_cell_bps>=1] = 1

Z = ward(pdist(cell_genotypes))
hclust_index = leaves_list(Z)
data = data['d_mat'][hclust_index]
n_bins = data.shape[1]
cell_genotypes = cell_genotypes.to_numpy()
cell_genotypes = cell_genotypes[hclust_index]
cell_bps = cell_bps.to_numpy()[hclust_index]

fig = plt.figure()
normalized_counts = data / np.sum(data, axis=1)[:, np.newaxis] * n_bins
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(normalized_counts, cmap=cmap, vmax=2)
plt.colorbar()
ax = plt.gca()
ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--')
plt.show()

plt.pcolor(cell_genotypes, cmap=cmap)
plt.colorbar()
ax = plt.gca()
ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--')
plt.show()
np.where(bps_indicator)

window_size = 20
freq = np.fft.fftfreq(n_bins, 1)
df = 0.015
gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
g = gpl + gmn

# Generate list of known breakpoints
input_breakpoints = np.array([200, 500])
bps = sci.detect_breakpoints(data, window_size=window_size, threshold=3.0, evaluate_peaks=True, compute_sp=True,  input_breakpoints=input_breakpoints)
bps['segmented_regions']

filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
fbps = sci.detect_breakpoints(data, window_size=window_size, threshold=4.0, lr=filtered_lr, evaluate_peaks=True, compute_sp=True, input_breakpoints=input_breakpoints)
fbps['segmented_regions']

plt.pcolor(normalized_counts, cmap=cmap, vmax=2)
plt.colorbar()
ax = plt.gca()
ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
plt.show()

segmented_regions = fbps['segmented_regions']
segmented_region_sizes = fbps['segmented_region_sizes']

condensed_data = sci.condense_regions(data, segmented_region_sizes)
condensed_normalised_data = sci.condense_regions(normalized_counts, segmented_region_sizes)

region_neutral_states = np.ones((condensed_data.shape[1]),) * 2
region_neutral_states[:3] = 1
region_neutral_states
sci.learn_tree(condensed_data, segmented_region_sizes, n_reps=5, full=False, cluster=True,
                cluster_tree_n_iters=10000, full_tree_n_iters=10000, max_tries=1, robustness_thr=0.5, max_scoring=True, verbosity=1, region_neutral_states=region_neutral_states)

sci.best_cluster_tree.plot_tree()

sci.best_full_tree.plot_tree()

plt.pcolor(sci.best_cluster_tree.outputs['inferred_cnvs'], cmap=cmap)
ax = plt.gca()
ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
plt.show()

plt.pcolor(sci.best_full_tree.outputs['inferred_cnvs'], cmap=cmap)
ax = plt.gca()
ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
plt.show()

sci.learn_tree(condensed_data, segmented_region_sizes, n_reps=1, full=True, cluster=True,
                cluster_tree_n_iters=10000, full_tree_n_iters=10000, max_tries=1, robustness_thr=0.5, max_scoring=False, verbosity=1)

sci.best_cluster_tree.plot_tree()

sci.best_full_tree.plot_tree()

plt.pcolor(sci.best_cluster_tree.outputs['inferred_cnvs'], cmap=cmap)
ax = plt.gca()
ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
plt.show()

plt.pcolor(sci.best_full_tree.outputs['inferred_cnvs'], cmap=cmap)
ax = plt.gca()
ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
plt.show()

threshold_coeffs[2::2]
for thres in threshold_coeffs[2::2]:

#
# sci.best_cluster_tree.plot_tree()
#
#
# sci.best_full_tree.plot_tree()
#
# plt.pcolor(sci.best_cluster_tree.outputs['inferred_cnvs'], cmap=cmap)
# ax = plt.gca()
# ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
# plt.show()
#
# plt.pcolor(sci.best_full_tree.outputs['inferred_cnvs'], cmap=cmap)
# ax = plt.gca()
# ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
# plt.show()
#
#
#
# nbps = sci.detect_breakpoints(data, window_size=window_size, threshold=3, sp=bps['sp_vec'])
# len(nbps['segmented_regions'])
#
# plt.plot(np.log(bps['sp_vec']))
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
# # ax.vlines(nbps['segmented_regions'], *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
#
# plt.pcolor(bps['lr_vec'].T, cmap=cmap, vmin=0, vmax=2)
# plt.colorbar()
#
# # Filter LR
# window_size = 20
# freq = np.fft.fftfreq(n_bins, 1)
# df = 0.015
# gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
# gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
# g = gpl + gmn
#
# filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
# plt.pcolor(filtered_lr, cmap=cmap, vmin=0, vmax=2)
# plt.colorbar()
#
# flbps = sci.detect_breakpoints(data, window_size=window_size, threshold=4, lr=filtered_lr)
# len(flbps['segmented_regions'])
#
# fig = plt.figure(figsize=(16,4))
# plt.plot(np.log(bps['sp_vec']))
# plt.plot(np.log(flbps['sp_vec']) - 15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
# ax.vlines(flbps['segmented_regions'], *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
#
#
#
#
# # ROC curve
# from sklearn.metrics import auc
# def generate_bp_roc(inferred_bps_dict, true_bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=True):
#     bps = pd.DataFrame(inferred_bps_dict['all_bps_comparison'])
#     bps.columns = ['idx','log_sp','range']
#     bps.sort_values('idx')['idx'].tolist()
#     bps.index = bps['idx']
#
#     bps['ranking'] = bps['log_sp'] / bps['range']
#     bps = bps.dropna()
#     min_ranking = np.min(bps['ranking'])
#     # threshold_coeffs = sorted(bps['ranking'].values)
#
#     # correcting for the bps 1-2 bins nearby
#     if correct_close:
#         for index, row in bps.iterrows():
#             idx_val = bps.loc[index, 'idx']
#             for gt in np.where(true_bps_indicator)[0]:
#                 if (abs(idx_val - gt) <= np.ceil(window_size/2) and idx_val != gt):
#                     # print('correcting ' + str(idx_val) + '->' + str(gt))
#                     bps.loc[index,'idx'] = gt
#
#     # Add remaining bins to bps to make sure all bins leads to TPR and FPR == 1
#     if add_dummy_bps:
#         new_bps_indicator = np.ones(true_bps_indicator.shape[0])
#         new_bps_indicator[bps['idx'].astype(int)] += 1
#         dummy_bps = np.where(new_bps_indicator==1)[0]
#
#         dummy_bps_df = pd.DataFrame({'idx': dummy_bps, 'log_sp': np.zeros(dummy_bps.shape[0]), 'range': np.zeros(dummy_bps.shape[0]), 'ranking': np.zeros(dummy_bps.shape[0])})
#         dummy_bps_df.index = dummy_bps_df['idx']
#         dummy_bps_df['ranking'] = min_ranking - 1.0
#         threshold_coeffs = np.append(np.min(dummy_bps_df['ranking']), threshold_coeffs)
#
#         bps = pd.concat([bps, dummy_bps_df])
#
#     tpr_values = []
#     fpr_values = []
#     fdr_values = []
#     roc_values = []
#     bps_indicators = []
#     for thr in threshold_coeffs:
#         inferred_bps = []
#         for index, row in bps.iterrows():
#             if row['ranking'] >= thr:
#                 inferred_bps.append(row['idx'])
#             else:
#                 break
#
#         inferred_bps_indicator = np.zeros(true_bps_indicator.shape[0])
#         inferred_bps_indicator[np.array(inferred_bps).astype(int)] = 1.
#         inferred_bps_indicator = np.array(inferred_bps_indicator)
#
#         tp = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==1)[0]]==1)
#         fp = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==0)[0]]==1)
#
#         tpr = tp / np.count_nonzero(true_bps_indicator==1)
#         fpr = fp / np.count_nonzero(true_bps_indicator==0)
#
#         # nothing was detected, maybe because the threshold was too high
#         # set fdr to zero
#         if (fp + tp != 0):
#             fdr = fp / (fp + tp)
#         else:
#             fdr = 0
#
#         tpr_values.append(tpr)
#         fpr_values.append(fpr)
#         fdr_values.append(fdr)
#         bps_indicators.append(inferred_bps_indicator)
#
#     roc_curve = dict()
#     roc_curve['tpr'] = tpr_values
#     roc_curve['fpr'] = fpr_values
#     roc_curve['fdr'] = fdr_values
#     roc_curve['bps'] = bps_indicators
#     roc_curve['auc'] = auc(fpr_values, tpr_values)
#
#     return roc_curve
#
# n_reps = 10
# nu = 4.
# n_cells = 200
# n_nodes = 10
# n_bins = 10000
# n_regions = 10
# n_reads_list = [2000, 4000, 8000]
# window_size = 20
# threshold_coeffs = np.linspace(0.1, 20.0, 100)
#
#
# def simulate_data(n_reads, postfix):
#     sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build',
#                       '/cluster/work/bewi/members/pedrof/sc-dna/scripts/', persistence=False, postfix=postfix)
#
#     # Generate data
#     data = sci.simulate_data(n_cells=n_cells, n_nodes=n_nodes, n_bins=n_bins, n_regions=n_regions,
#                             n_reads=n_reads, nu=nu, max_regions_per_node=10, min_reg_size=2*window_size)
#
#     # Get true breakpoints
#     cell_genotypes = pd.DataFrame(data['ground_truth'])
#     cell_bps = cell_genotypes.diff(periods=1, axis=1)
#     cell_bps = cell_bps.fillna(value=0.0)
#     cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
#     grouped_cell_bps = cell_bps.sum(axis=0)
#     ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
#     ground_truth = ground_truth.index.tolist()
#
#     bps_indicator = np.zeros(grouped_cell_bps.shape[0])
#     bps_indicator[grouped_cell_bps>=1] = 1
#
#     return data['d_mat'], bps_indicator
#
# def get_roc_curves(data, bps_indicator, postfix):
#     sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build',
#                       '/cluster/work/bewi/members/pedrof/sc-dna/scripts/', persistence=False, postfix=postfix)
#
#     # Infer breakpoints with new and old methods
#     bps = sci.detect_breakpoints(data, window_size=window_size, threshold=0.1, bp_limit=data.shape[1])
#     basic_roc_curve = generate_bp_roc(bps, bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=False)
#
#     filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
#     bps = sci.detect_breakpoints(data, window_size=window_size, threshold=0.1, lr=filtered_lr)
#     filtered_roc_curve = generate_bp_roc(bps, bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=False)
#
#     print("Done!")
#
#     return basic_roc_curve, filtered_roc_curve
#
# nu = 5
# n_bins = 10000
# data, bps_indicator = simulate_data(40000, "")
# sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build',
#                   '/cluster/work/bewi/members/pedrof/sc-dna/scripts/', persistence=False, postfix="heyhey")
# # Infer breakpoints with new and old methods
# bps = sci.detect_breakpoints(data, window_size=window_size, threshold=3, bp_limit=data.shape[1], compute_sp=False, evaluate_peaks=False)
#
# window_size = 20
# freq = np.fft.fftfreq(n_bins, 1)
# df = 0.015
# gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
# gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
# g = gpl + gmn
# filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
#
# flbps = sci.detect_breakpoints(data, window_size=window_size, threshold=3, lr=filtered_lr)
# len(flbps['segmented_regions'])
# logsp = flbps['all_bps_comparison'][:, 1]
# logsprange = flbps['all_bps_comparison'][:, 2]
# ranking = np.argsort(logsp/logsprange)[::-1]
#
# roc_curve = generate_bp_roc(flbps, bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=False)
# flbps['segmented_regions']
# np.where(bps_indicator)[0]
#
# np.count_nonzero(logsp/logsprange > 20)
# plt.plot(logsp/logsprange, marker='.')
# # plt.ylim([0, 100])
# # plt.xlim([0, 20])
# newbps = flbps['all_bps_comparison'][:,0][ranking][:10]
#
# plt.figure(figsize=(16,4))
# plt.plot(np.log(flbps['sp_vec']))
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6, color='black')
# ax.vlines(newbps, *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
# newbps[3]
# fig = plt.figure(figsize=(16,4))
# plt.pcolor(bps['lr_vec'].T[:, 7000:8000], cmap=cmap, vmin=0)
# plt.colorbar()
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator[7000:8000])[0], *ax.get_ylim(), linestyle='--', alpha=0.6, color='yellow')
# # ax.vlines(newbps[7:9]-7000, *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
#
# newbps[3]
# fig = plt.figure(figsize=(16,4))
# plt.pcolor(bps['lr_vec'].T[:, 2500:4000], cmap=cmap, vmin=0)
# plt.colorbar()
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0][2500:4000], *ax.get_ylim(), linestyle='--', alpha=0.6, color='black')
# ax.vlines(newbps[3]-2500, *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
#
# for n_reads in n_reads_list:
#     print('{} reads...'.format(n_reads))
#
#     basic_roc_curves = []
#     filtered_roc_curves = []
#     for i in range(n_reps):
#         print(f"Rep. number {i+1}")
#         data, bps_indicator = simulate_data(n_reads, f"PYSCICONETEMP_{i}_{n_reads}")
#         basic_roc, filtered_roc = get_roc_curves(data_list[i], bps_indicator_list[i], f"PYSCICONETEMP_{i}_{n_reads}")
#         basic_roc_curves.append(basic_roc)
#         filtered_roc_curves.append(filtered_roc)
#
#     for result in results:
#         basic, filtered = result.get()
#         basic_roc_curves.append(basic)
#         filtered_roc_curves.append(filtered)
#
#     basic_fprs = []
#     basic_tprs = []
#     for roc_curve in basic_roc_curves:
#         basic_fprs.append(roc_curve['fpr'])
#         basic_tprs.append(roc_curve['tpr'])
#     filtered_fprs = []
#     filtered_tprs = []
#     for roc_curve in filtered_roc_curves:
#         filtered_fprs.append(roc_curve['fpr'])
#         filtered_tprs.append(roc_curve['tpr'])
#     basic_mean_fprs = np.mean(basic_fprs, axis=0)
#     basic_mean_tprs = np.mean(basic_tprs, axis=0)
#     filtered_mean_fprs = np.mean(filtered_fprs, axis=0)
#     filtered_mean_tprs = np.mean(filtered_tprs, axis=0)
#
#     mean_fprs = filtered_mean_fprs
#     mean_tprs = filtered_mean_tprs
#     plt.figure(figsize=(8,8))
#     plt.plot(mean_fprs, mean_tprs, color="black", alpha=1)
#     cmap = sns.color_palette("Dark2", n_colors=8)
#     plt.scatter(mean_fprs[5], mean_tprs[5], color=cmap[0], alpha=1, s=100, label='1')
#     plt.scatter(mean_fprs[7], mean_tprs[7], color=cmap[1], alpha=1, s=100, label='1.5')
#     plt.scatter(mean_fprs[10], mean_tprs[10], color=cmap[2], alpha=1, s=100, label='2')
#     plt.scatter(mean_fprs[12], mean_tprs[12], color=cmap[3], alpha=1, s=100, label='2.5')
#     plt.scatter(mean_fprs[15], mean_tprs[15], color=cmap[4], alpha=1, s=100, label='3')
#     plt.scatter(mean_fprs[17], mean_tprs[17], color=cmap[5], alpha=1, s=100, label='3.5')
#     plt.scatter(mean_fprs[20], mean_tprs[20], color=cmap[6], alpha=1, s=100, label='4')
#     plt.scatter(mean_fprs[22], mean_tprs[22], color=cmap[7], alpha=1, s=100, label='4.5')
#     plt.grid(True)
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.title(f"Mean ROC curves ({n_reads/n_bins}x)")
#     plt.ylim(0, 1.01)
#     plt.legend()
#     plt.show()
#
# plt.figure(figsize=(8,8))
# plt.plot(basic_mean_fprs, basic_mean_tprs, alpha=1, label='basic')
# plt.plot(filtered_mean_fprs, filtered_mean_tprs, alpha=1, label='filtered')
# plt.grid(True)
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title(f"Mean ROC curves")
# plt.ylim(0, 1.01)
# plt.legend()
# plt.show()
#
#
# mean_fprs = filtered_mean_fprs
# mean_tprs = filtered_mean_tprs
# plt.figure(figsize=(8,8))
# plt.plot(mean_fprs, mean_tprs, color="black", alpha=1)
# cmap = sns.color_palette("Dark2", n_colors=8)
# plt.scatter(mean_fprs[5], mean_tprs[5], color=cmap[0], alpha=1, s=100, label='1')
# plt.scatter(mean_fprs[7], mean_tprs[7], color=cmap[1], alpha=1, s=100, label='1.5')
# plt.scatter(mean_fprs[10], mean_tprs[10], color=cmap[2], alpha=1, s=100, label='2')
# plt.scatter(mean_fprs[12], mean_tprs[12], color=cmap[3], alpha=1, s=100, label='2.5')
# plt.scatter(mean_fprs[15], mean_tprs[15], color=cmap[4], alpha=1, s=100, label='3')
# plt.scatter(mean_fprs[17], mean_tprs[17], color=cmap[5], alpha=1, s=100, label='3.5')
# plt.scatter(mean_fprs[20], mean_tprs[20], color=cmap[6], alpha=1, s=100, label='4')
# plt.scatter(mean_fprs[22], mean_tprs[22], color=cmap[7], alpha=1, s=100, label='4.5')
# plt.grid(True)
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title(f"Mean ROC curves ({n_reads/n_bins}x)")
# plt.ylim(0, 1.01)
# plt.legend()
# plt.show()
#
#
#
#
#
#
# plt.pcolor(bps['lr_vec'].T[:,300:400], cmap=cmap, vmin=0, vmax=2)
# plt.colorbar()
#
# # Filter Sp
# import pywt
# sig = np.log(flbps['sp_vec'])
# noise = np.percentile(sig, 75) - np.percentile(sig, 25)
# # Create wavelet object and define parameters
# w = pywt.Wavelet('db4')
# maxlev = pywt.dwt_max_level(len(sig), w.dec_len)
# # maxlev = 2 # Override if desired
# print("maximum level is " + str(maxlev))
# threshold = 0.2 * noise # Threshold for filtering
# # Decompose into wavelet components, to the level selected:
# coeffs = pywt.wavedec(sig, 'db4', level=maxlev)
# # plt.figure()
# for i in range(1, len(coeffs)):
#     # plt.subplot(maxlev, 1, i)
#     # plt.plot(coeffs[i])
#     coeffs[i] = pywt.threshold(coeffs[i], threshold*max(coeffs[i]))
#     # plt.plot(coeffs[i])
# datarec = pywt.waverec(coeffs, 'db4')
#
# fig = plt.figure(figsize=(16,4))
# plt.plot(sig)
# plt.plot(datarec.real-15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
#
# normedseg = sig[2049+20:3142-20]/np.median(sig[2049+20:3142-20])
# r = np.percentile(normedseg, 75) - np.median(normedseg)
# plt.plot(normedseg)
# ax = plt.gca()
# ax.hlines(1 + 3*r, xmin=0, xmax=normedseg.shape[0])
#
# normedseg = datarec.real[2049+20:3142-20]/np.median(datarec.real[2049+20:3142-20])
# r = np.percentile(normedseg, 75) - np.percentile(normedseg, 25)
# plt.plot(normedseg)
# ax = plt.gca()
# ax.hlines(1 + 3*r, xmin=0, xmax=normedseg.shape[0])
#
#
# np.where(bps_indicator)
# fbps = sci.detect_breakpoints(data, window_size=window_size, threshold=4, sp=np.exp(datarec.real))
# len(fbps['segmented_regions'])
#
# logsp = flbps['all_bps_comparison'][:, 1]
# logsprange = flbps['all_bps_comparison'][:, 2]
# ranking = np.argsort(logsp/logsprange)[::-1]
#
# plt.plot((logsp/logsprange)[ranking][:20], marker='.')
#
#
# fig = plt.figure(figsize=(16,8))
# plt.plot(sig)
# plt.plot(datarec.real-15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6, color='black')
# ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--', alpha=0.6, color='blue')
#
# sig = datarec.real
# freq = np.fft.fftfreq(n_bins, 1) # 1 Hz sampling rate
# X = np.fft.fft(sig - np.mean(sig))
# plt.plot(freq, np.abs(X))
# plt.xlim([0, 0.2])
# np.where(np.isclose(freq, 0.125))
# H = np.zeros(X.shape)
# H[:200], H[-200:] = 1, 1
# plt.plot(freq, np.abs(X))
# plt.plot(freq, H*np.max(np.abs(X)))
# Y = X * H
# y = np.fft.ifft(Y)
#
# fig = plt.figure(figsize=(16,4))
# plt.plot(sig)
# plt.plot(y.real-15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
# normedseg = y[2049+20:3142-20]/np.median(y[2049+20:3142-20])
# r = np.percentile(normedseg, 75) - np.median(normedseg)
# plt.plot(normedseg)
# ax = plt.gca()
# ax.hlines(1 + 3*r, xmin=0, xmax=normedseg.shape[0])
#
#
# plt.plot(np.exp(np.abs(y.real)))
#
# fbps = sci.detect_breakpoints(data, window_size=window_size, threshold=3, sp=flbps['sp_vec'])
# len(fbps['segmented_regions'])
#
# plt.pcolor(normalized_counts, cmap=cmap, vmax=2)
# plt.colorbar()
# ax = plt.gca()
# ax.vlines(fbps['segmented_regions'], *ax.get_ylim(), linestyle='--')
# plt.show()
#
# segmented_regions = filtered_bps['segmented_regions']
# segmented_region_sizes = filtered_bps['segmented_region_sizes']
#
# condensed_data = sci.condense_regions(data, segmented_region_sizes)
#
# sci.learn_tree(condensed_data, segmented_region_sizes, n_reps=10, full=True, cluster=True, cluster_tree_n_iters=4000, full_tree_n_iters=0000, max_tries=3, robustness_thr=0.5)
#
# sci.best_cluster_tree.plot_tree()
#
# sci.best_full_tree.plot_tree()
#
# plt.pcolor(sci.best_cluster_tree.outputs['inferred_cnvs'], cmap=cmap)
# ax = plt.gca()
# ax.vlines(filtered_bps['segmented_regions'], *ax.get_ylim(), linestyle='--')
# plt.show()
#
#
# plt.pcolor(sci.best_full_tree.outputs['inferred_cnvs'], cmap=cmap)
# ax = plt.gca()
# ax.vlines(filtered_bps['segmented_regions'], *ax.get_ylim(), linestyle='--')
# plt.show()
#
# #####
#
# mat = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/simulation_/10nodes_10regions_20000reads/0_d_mat.txt', delimiter=',')
# lr_vec = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/bp_detection__trees/10nodes_10regions_20000reads/0_lr_vec.csv', delimiter=',')
# gt = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/simulation_/10nodes_10regions_20000reads/0_ground_truth.txt', delimiter=',')
# n_bins = mat.shape[1]
#
# cell_genotypes = pd.DataFrame(gt)
# cell_bps = cell_genotypes.diff(periods=1, axis=1)
# cell_bps = cell_bps.fillna(value=0.0)
# cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
# grouped_cell_bps = cell_bps.sum(axis=0)
# ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
# ground_truth = ground_truth.index.tolist()
#
# bps_indicator = np.zeros(grouped_cell_bps.shape[0])
# bps_indicator[grouped_cell_bps>=1] = 1
#
# def filter_lr(lr_matrix, H=None):
#     freq = np.fft.fftfreq(lr_matrix.shape[-1], 1) # 1 Hz sampling rate
#
#     filtered_lr = np.empty(lr_matrix.shape)
#     for c in range(lr_matrix.shape[0]):
#         X = np.fft.fft(lr_matrix[c])
#         Y = X * H
#         y = np.fft.ifft(Y)
#         filtered_lr[c] = y
#
#     return filtered_lr
#
# window_size = 20
# freq = np.fft.fftfreq(n_bins, 1)
# df = 0.02
# gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
# gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
# g = gpl + gmn
#
#
# fig = plt.figure(figsize=(26,4))
# plt.pcolor(lr_vec.T[:, :5000], cmap=cmap, vmin=0, vmax=40)
# plt.colorbar()
#
# segmented_regions = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/bp_detection__trees/10nodes_10regions_20000reads/0_segmented_regions.txt', delimiter=',')
# sp_vec = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/bp_detection__trees/10nodes_10regions_20000reads/0_sp_vec.csv', delimiter=',')
# log_sp_vec = np.log(sp_vec)
# plt.plot(sp_vec)
# fig = plt.figure(figsize=(26,4))
# plt.plot(log_sp_vec)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--')
# # ax.vlines(segmented_regions, *ax.get_ylim(), linestyle='--', color='blue', alpha=0.6)
#
#
# all_bps_comparison = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/bp_detection__trees/10nodes_10regions_20000reads/0_all_bps_comparison.csv', delimiter=',')
# bps = pd.DataFrame(all_bps_comparison)
# bps.columns = ['idx','log_sp','range']
# bps.sort_values('idx')['idx'].tolist()
# bps.index = bps['idx']
#
# bps['ranking'] = bps['log_sp'] / bps['range']
# bps = bps.dropna()
# min_ranking = np.min(bps['ranking'])
# # threshold_coeffs = sorted(bps['ranking'].values)
#
# thr = 10
# inferred_bps = []
# for index, row in bps.iterrows():
#     if row['ranking'] >= thr:
#         inferred_bps.append(row['idx'])
#     else:
#         break
#
# inferred_bps_indicator = np.zeros(n_bins)
# inferred_bps_indicator[np.array(inferred_bps).astype(int)] = 1.
# inferred_bps_indicator = np.array(inferred_bps_indicator)
#
# np.count_nonzero(inferred_bps_indicator)
#
#
# freq = np.fft.fftfreq(n_bins, 1) # 1 Hz sampling rate
# X = np.fft.fft(log_sp_vec - np.mean(log_sp_vec))
#
# plt.plot(freq, np.abs(X))
# plt.xlim([0, 0.2])
#
# np.where(np.isclose(freq, 0.125))
#
# H = np.zeros(X.shape)
# H[:300], H[-300:] = 1, 1
# plt.plot(freq, np.abs(X))
# plt.plot(freq, H*np.max(np.abs(X)))
#
# Y = X * H
# y = np.fft.ifft(Y)
#
# plt.plot(freq, np.abs(Y))
# plt.xlim([0, 0.2])
#
# plt.figure(figsize=(26,8))
# plt.plot(log_sp_vec + 15)
# plt.plot(y + np.mean(log_sp_vec))
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
#
# plt.figure(figsize=(26,8))
# plt.plot(log_sp_vec[5000:7000] + 15)
# plt.plot(y[5000:7000] + np.mean(log_sp_vec))
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator[5000:7000])[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
# # The cost of denoising based on frequencies is the removal of edges. Let's try a wavelet instead
#
# import pywt
#
# noise = np.percentile(log_sp_vec, 75) - np.percentile(log_sp_vec, 50)
# noise
#
# # Create wavelet object and define parameters
# w = pywt.Wavelet('db4')
# maxlev = pywt.dwt_max_level(len(log_sp_vec), w.dec_len)
# # maxlev = 2 # Override if desired
# print("maximum level is " + str(maxlev))
# threshold = 0.5*noise # Threshold for filtering
#
# # Decompose into wavelet components, to the level selected:
# coeffs = pywt.wavedec(log_sp_vec, 'db4', level=maxlev)
#
# #cA = pywt.threshold(cA, threshold*max(cA))
# plt.figure()
# for i in range(1, len(coeffs)):
#     plt.subplot(maxlev, 1, i)
#     plt.plot(coeffs[i])
#     coeffs[i] = pywt.threshold(coeffs[i], threshold*max(coeffs[i]))
#     plt.plot(coeffs[i])
#
#
# datarec = pywt.waverec(coeffs, 'db4')
#
# plt.figure(figsize=(26,8))
# plt.plot(log_sp_vec[5000:7000] + 15)
# plt.plot(datarec[5000:7000])
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator[5000:7000])[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
# plt.xticks(np.arange(0, 2000, 20))
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(26,8))
# plt.plot(log_sp_vec + 15)
# plt.plot(y + np.mean(log_sp_vec))
# plt.plot(datarec- 15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
# # The wavelet still has this really high frequency oscilations, which I can remove with the fft
# freq = np.fft.fftfreq(n_bins, 1) # 1 Hz sampling rate
# X = np.fft.fft(datarec - np.mean(datarec))
#
# plt.plot(freq, np.abs(X))
# plt.xlim([0, 0.2])
#
# H = np.zeros(X.shape)
# H[:600], H[-600:] = 1, 1
# plt.plot(freq, np.abs(X))
# plt.plot(freq, H*np.max(np.abs(X)))
#
# Y = X * H
# y = np.fft.ifft(Y)
#
# plt.figure(figsize=(26, 8))
# plt.plot(datarec)
# plt.plot(y- 15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
#
# plt.figure(figsize=(26,8))
# plt.plot(log_sp_vec[5000:7000]+15)
# plt.plot(datarec[5000:7000])
# plt.plot(y[5000:7000] - 15)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator[5000:7000])[0], *ax.get_ylim(), linestyle='--', alpha=0.6)
# plt.xticks(np.arange(0, 2000, 20))
# plt.grid()
# plt.show()
#
#
# len(coeffs)
# filtered_coeffs = pywt.threshold(coeffs, 2*noise, 'soft')
# filtered_log_sp_vec = pywt.waverec(filtered_coeffs, 'coif5')
#
# fig = plt.figure(figsize=(26,4))
# plt.plot(filtered_log_sp_vec)
#
# plt.hist(np.log(sp_vec))
#
#
# rm_filtered_lr = filter_lr(lr_vec.T, H=g)
# rm_filtered_bps = sci.detect_breakpoints(mat, window_size=window_size, threshold=5, bp_limit=n_bins, lr=lr_vec.T)
# rm_filtered_bps['segmented_region_sizes'].shape
#
#
# # Get true breakpoints
# gt = np.loadtxt('/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/simulation_/10nodes_10regions_20000reads/0_ground_truth.txt', delimiter=',')
# cell_genotypes = pd.DataFrame(gt)
# cell_bps = cell_genotypes.diff(periods=1, axis=1)
# cell_bps = cell_bps.fillna(value=0.0)
# cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
# grouped_cell_bps = cell_bps.sum(axis=0)
# ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
# ground_truth = ground_truth.index.tolist()
#
# bps_indicator = np.zeros(grouped_cell_bps.shape[0])
# bps_indicator[grouped_cell_bps>=1] = 1
#
# Z = ward(pdist(cell_genotypes))
# hclust_index = leaves_list(Z)
# data = mat[hclust_index]
# n_bins = data.shape[1]
# cell_genotypes = cell_genotypes.to_numpy()
# cell_genotypes = cell_genotypes[hclust_index]
# cell_bps = cell_bps.to_numpy()[hclust_index]
#
# normalized_counts = data / np.sum(data, axis=1)[:, np.newaxis] * n_bins
# cmap = sns.diverging_palette(220, 10, as_cmap=True)
# plt.pcolor(normalized_counts, cmap=cmap, vmax=2)
# plt.colorbar()
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--')
# plt.show()
#
# plt.pcolor(cell_genotypes, cmap=cmap)
# ax = plt.gca()
# ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), linestyle='--')
# plt.show()
#
# #####
#
# def sweep_thresholds(inferred_bps_dict, threshold_coeffs, window_size, n_bins):
#     bps = pd.DataFrame(inferred_bps_dict['all_bps_comparison'])
#     bps.columns = ['idx','log_sp','range']
#     bps.sort_values('idx')['idx'].tolist()
#     bps.index = bps['idx']
#
#     bps['ranking'] = bps['log_sp'] / bps['range']
#     bps = bps.dropna()
#     min_ranking = np.min(bps['ranking'])
#     # threshold_coeffs = sorted(bps['ranking'].values)
#
#     bps_indicators = []
#     for thr in threshold_coeffs:
#         inferred_bps = []
#         for index, row in bps.iterrows():
#             if row['ranking'] >= thr:
#                 inferred_bps.append(row['idx'])
#             else:
#                 break
#
#         inferred_bps_indicator = np.zeros(n_bins)
#         inferred_bps_indicator[np.array(inferred_bps).astype(int)] = 1.
#         inferred_bps_indicator = np.array(inferred_bps_indicator)
#
#         bps_indicators.append(inferred_bps_indicator)
#
#     roc_curve = dict()
#     roc_curve['bps'] = bps_indicators
#
#     return roc_curve
#
# mat = np.loadtxt('/cluster/work/bewi/members/pedrof/data/SA501/cnv/cell_corrected_bin_filtered.txt', delimiter=',')
# lr_vec = np.loadtxt('/cluster/work/bewi/members/pedrof/vancouver_analysis/snake_analysis_files/vancouver_unfiltered_lr_vec.csv', delimiter=',')
# n_bins = mat.shape[1]
#
# def filter_lr(lr_matrix, H=None):
#     freq = np.fft.fftfreq(lr_matrix.shape[-1], 1) # 1 Hz sampling rate
#
#     filtered_lr = np.empty(lr_matrix.shape)
#     for c in range(lr_matrix.shape[0]):
#         X = np.fft.fft(lr_matrix[c])
#         Y = X * H
#         y = np.fft.ifft(Y)
#         filtered_lr[c] = y
#
#     return filtered_lr
#
# window_size = 20
# freq = np.fft.fftfreq(n_bins, 1)
# df = 0.02
# gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
# gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
# g = gpl + gmn
#
#
# rm_filtered_lr = filter_lr(lr_vec.T, H=g)
# rm_filtered_bps = sci.detect_breakpoints(mat, window_size=window_size, threshold=3, bp_limit=n_bins, lr=rm_filtered_lr)
# rm_filtered_bps['segmented_region_sizes'].shape
#
# thres_coeffs = [1, 2, 3, 4]
# rm_filtered_bps_sweep = sweep_thresholds(rm_filtered_bps, thres_coeffs, window_size, n_bins)
#
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# normalized_counts = mat / np.sum(mat, axis=1)[:, np.newaxis] * n_bins
#
# [(np.count_nonzero(rm_filtered_bps_sweep['bps'][i]), thres) for i, thres in enumerate(thres_coeffs)]
#
#
#
# i = 2
# thres = thres_coeffs[i]
# cmap = sns.diverging_palette(220, 10, as_cmap=True)
# fig = plt.figure(figsize=(16,4))
# plt.pcolor(normalized_counts, cmap=cmap, vmax=2)
# ax = plt.gca()
# ax.vlines(np.where(rm_filtered_bps_sweep['bps'][i])[0], *ax.get_ylim(), linestyle='--')
# plt.title(f'Threshold={thres}')
# plt.savefig('/cluster/work/bewi/members/pedrof/{}vancouver_lr_vec.png'.format(i), bbox_inches='tight')
