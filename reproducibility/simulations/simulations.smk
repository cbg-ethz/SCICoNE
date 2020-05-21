import os
import itertools
import pandas as pd
import numpy as np
from collections import Counter
from sklearn.metrics import roc_curve, auc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from tqdm import tqdm as tqdm
import re
import phenograph
from scgenpy import *
from scgenpy.preprocessing.utils import *
sns.set(rc={'figure.figsize':(15.7,8.27)})

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

class Tree(object):
    def __init__(self, binary_path, output_path, postfix='PYSCICONETREETEMP', persistence=False, ploidy=2, copy_number_limit=2):
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

        self.outputs = dict(tree_inferred=None, inferred_cnvs=None,
                            rel_markov_chain=None, cell_node_ids=None,
                            cell_region_cnvs=None, acceptance_ratio=None,
                            gamma_values=None)

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
                    n_nodes=3,  seed=42, postfix="", initial_tree=None, nu=1.0, cluster_sizes=None, alpha=0., verbosity=1):
        if postfix == "":
            postfix = self.postfix

        n_cells, n_regions = segmented_data.shape
        move_probs_str = ",".join(str(p) for p in move_probs)

        # Save temporary files
        temp_segmented_data_file = f"{postfix}_temp_segmented_data.txt"
        np.savetxt(temp_segmented_data_file, segmented_data, delimiter=',')

        temp_segmented_region_sizes_file = f"{postfix}_temp_segmented_region_sizes.txt"
        np.savetxt(temp_segmented_region_sizes_file, segmented_region_sizes, delimiter=',')

        if cluster_sizes is None:
            cluster_sizes = np.ones((n_cells,))
        temp_cluster_sizes_file = f"{postfix}_temp_cluster_sizes.txt"
        np.savetxt(temp_cluster_sizes_file, cluster_sizes, delimiter=',')

        if initial_tree is not None:
            temp_tree_file = f"{postfix}_temp_tree.txt"
            f = open(temp_tree_file, "w")
            f.write(initial_tree.tree_str)

            try:
                cmd_output = subprocess.run([self.binary_path, f"--d_matrix_file={temp_segmented_data_file}", f"--n_regions={n_regions}",\
                    f"--n_cells={n_cells}", f"--ploidy={self.ploidy}", f"--verbosity={verbosity}", f"--postfix={postfix}",\
                    f"--copy_number_limit={self.copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--tree_file={temp_tree_file}", f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}",\
                    "--max_scoring=true"])
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
                    f"--copy_number_limit={self.copy_number_limit}", f"--n_iters={n_iters}", f"--n_nodes={n_nodes}",\
                    f"--move_probs={move_probs_str}", f"--seed={seed}", f"--region_sizes_file={temp_segmented_region_sizes_file}",\
                    f"--nu={nu}", f"--cluster_sizes_file={temp_cluster_sizes_file}", f"--alpha={alpha}", "--max_scoring=true"])
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.remove(temp_segmented_data_file)
        os.remove(temp_segmented_region_sizes_file)
        os.remove(temp_cluster_sizes_file)

        # Read in all the outputs
        try:
            cwd = os.getcwd()
            for fn in os.listdir(cwd):
                if postfix in fn and (".csv" in fn or ".tsv" in fn): # Read in all data
                    for key in self.outputs.keys():
                        if key in fn:
                            delim = ','
                            if ".tsv" in fn:
                                delim = '\t'
                            self.outputs[key] = np.loadtxt(fn, delimiter=delim)
                            os.remove(fn)
                elif postfix in fn and "tree_inferred.txt" in fn: # Parse tree structure, score and nu
                    self.read_from_file(fn)
                    os.remove(fn)
        except OSError as e:
            pass
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

    def simulate_data(self, n_cells=200, n_nodes=5, n_bins=1000, n_regions=40, n_reads=10000, nu=1.0, min_reg_size=10, max_regions_per_node=1, ploidy=2, verbosity=0):
        output_path = os.path.join(self.output_path, f"{self.postfix}_simulation")

        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        cwd = os.getcwd()

        done = False
        while not done:
            try:
                cmd_output = subprocess.run([self.simulation_binary, f"--n_cells={n_cells}", f"--n_nodes={n_nodes}",\
                    f"--n_regions={n_regions}", f"--n_bins={n_bins}", f"--n_reads={n_reads}", f"--nu={nu}",\
                    f"--min_reg_size={min_reg_size}", f"--max_regions_per_node={max_regions_per_node}",\
                    f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={self.postfix}"])
            except subprocess.SubprocessError as e:
                print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            # The binary generates 4 files: *_d_mat.csv, *_ground_truth.csv, *_region_sizes.txt and *_tree.txt.
            # We read each one of these files into np.array, np.array, np.array and Tree structures, respectively,
            # and output a dictionary with keys "d_mat", "ground_truth", "region_sizes" and "tree".
            d_mat_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv"
            ground_truth_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv"
            region_sizes_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt"
            tree_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt"

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

    def _detect_breakpoints(self, data, window_size=30, threshold=3.0, bp_limit=300, lr=None):
        n_cells = data.shape[0]
        n_bins = data.shape[1]
        verbosity = 1 # > 0 to generate all files

        compute_lr = True
        lr_file = f"{self.postfix}_pre_lr_vec.txt"
        if lr is not None:
            compute_lr = False
            # Use provided lr values instead of computing new ones
            if lr.shape[0] < lr.shape[1]:
                lr = lr.T # make sure it's bins by cells
            np.savetxt(lr_file, lr, delimiter=',')

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{self.postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')

            cmd_output = subprocess.run([self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}",\
                f"--bp_limit={bp_limit}",\
                f"--verbosity={verbosity}", f"--postfix={self.postfix}"])

            # Delete the data file
            os.remove(data_file)

            if not compute_lr:
                os.remove(lr_file)
        except subprocess.SubprocessError as e:
            print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

        output = dict()

        try:
            cwd = os.getcwd()
            for fn in os.listdir(cwd):
                if self.postfix in fn:
                    key = fn.split('_', 1)[1].split('.')[0]
                    output[key] = np.loadtxt(fn, delimiter=',')
                    if not self.persistence:
                        os.remove(fn)
        except OSError as e:
            print("OSError: ", e.output, e.stdout, e.stderr)

        return output

    def detect_breakpoints(self, data, window_size=30, threshold=3.0, bp_limit=300, lr=None, sp=None, evaluate_peaks=True, compute_lr=True, compute_sp=True):
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

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{self.postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')

            cmd_output = subprocess.run([self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}",\
                f"--bp_limit={bp_limit}", f"--compute_lr={compute_lr}", f"--lr_file={lr_file}",\
                f"--compute_sp={compute_sp}", f"--sp_file={sp_file}", f"--verbosity={verbosity}",\
                f"--evaluate_peaks={evaluate_peaks}", f"--postfix={self.postfix}"])

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
        pool = ThreadPool(n_reps)
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

'''
parameters
'''

'''
Trees of size n=10, 20 and 40
n_regions = n, 2n and 4n
10,000 bins and 400 cells
20,000, 40,000 and 80,000 reads per cell
'''

all_n_tps = config["simulate"]["n_reps"]
n_inference_reps = config["inference"]["full_trees"]["n_reps"]
tree_rep = config["inference"]["cluster_trees"]["n_reps"]

n_nodes = config["simulate"]["n_nodes"] # values: [10,20,30]
n_regions = [n_nodes, 2*n_nodes, 4*n_nodes] # values: [n, 2n,4n]
n_bins = config["simulate"]["n_bins"]
n_reads = [2*n_bins, 4*n_bins, 8*n_bins]
nu = config["simulate"]["nu"]
window_size = config["bp_detection"]["window_size"]

n_cells = config["simulate"]["n_cells"]
n_iters = config["inference"]["full_trees"]["n_iters"]  # int(1000000*n_nodes/10)

sim_output_file_exts = ['d_mat.txt','ground_truth.txt','region_sizes.txt', 'tree.txt'] #, 'inferred_cnvs.txt', 'tree_inferred.txt', 'HMMcopy_inferred.txt','inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']
bp_output_file_exts = ['segmented_regions.txt', 'segmented_region_sizes.txt', 'sp_vec.csv', 'log_posterior_vec.csv', 'expected_k_vec.csv', 'all_bps_comparison.csv', 'lr_vec.csv']
# trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt'] # , 'inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

output_path = config["output"]
SIM_OUTPUT= os.path.join(config["output"], "simulation")
BP_OUTPUT = os.path.join(config["output"], "bp_detection")
PHENO_OUTPUT = os.path.join(config["output"], "phenograph")
TREES_OUTPUT = os.path.join(config["output"], "inference")
HMMCOPY_OUTPUT = os.path.join(config["output"], "hmm_copy")
HCLUST_OUTPUT = os.path.join(config["output"], "hclust")
sim_prefix=config["simulate"]["prefix"]

output_temp = config["inference"]["output_temp"]
try:
    os.mkdir(output_temp)
except:
    pass

binaries_path = config['binaries_path']

# the default case for cluster_fraction variable
try:
    cf = config["inference"]["full_trees"]["cluster_fraction"]
except KeyError:
    cf = 1.0

freq = np.fft.fftfreq(n_bins, 1)
df = 0.015
gpl = np.exp(- ((freq-1/(2*window_size))/(2*df))**2)  # pos. frequencies
gmn = np.exp(- ((freq+1/(2*window_size))/(2*df))**2)
g = gpl + gmn

rule all:
    input:
        simulations = expand(f'{SIM_OUTPUT}_{sim_prefix}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_d_mat.txt',
                                        regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        best_full_tree = expand(f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
                        regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        best_cluster_tree = expand(f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
                            regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        hmmcopy_inferred_cnvs = expand(HMMCOPY_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_hmm_copy_inferred.txt',
                                regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        hclust_inferred_cnvs = expand(HCLUST_OUTPUT+'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_hclust_inferred.txt',
                                regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        phenograph_assignments = expand(PHENO_OUTPUT+'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_phenograph_inferred.txt',
                                 regions=[x for x in n_regions],reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)])
    run:
        print("rule all")

rule run_sim:
    params:
        sim_bin = config["simulate"]["bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        all_n_tps = all_n_tps,
        n_iters = n_iters,
        nu = nu,
        min_reg_size = window_size,
        max_regions_per_node = 10
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ground_truth.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        tree = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree.txt'
    run:
        done = False
        while not done:
            try:
                cmd_output = subprocess.run([params.sim_bin, f"--n_cells={params.n_cells}", f"--n_nodes={params.n_nodes}",\
                    f"--n_regions={wildcards.regions}", f"--n_bins={params.n_bins}", f"--n_reads={wildcards.reads}", f"--nu={params.nu}",\
                    f"--min_reg_size={params.min_reg_size}", f"--max_regions_per_node={params.max_regions_per_node}",\
                    f"--ploidy=2", f"--verbosity=0", f"--postfix={wildcards.rep_id}", f"--seed={wildcards.rep_id}"])
            except subprocess.SubprocessError as e:
                print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            try:
                # Move to output directory
                os.rename(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.csv", output.d_mat)
                os.rename(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_ground_truth.csv", output.ground_truth)
                os.rename(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_region_sizes.txt", output.region_sizes)
                os.rename(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_tree.txt", output.tree)
                done = True
            except OSError as e:
                print(e)
                print('Retrying...')
                pass

rule detect_breakpoints:
    params:
        binary = config["bp_detection"]["bin"],
        window_size = config["bp_detection"]["window_size"],
        verbosity = config["bp_detection"]["verbosity"],
        threshold = config["bp_detection"]["threshold"],
        bp_limit = config["bp_detection"]["bp_limit"],
        postfix = sim_prefix
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        segmented_regions = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_regions.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt',

        # sp_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_sp_vec.csv',
        # log_posterior_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_log_posterior_vec.csv',
        # expected_k_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_expected_k_vec.csv',
        # all_bps_comparison = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_all_bps_comparison.csv',
        # lr_vec = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_lr_vec.csv'
    run:
        try:
            os.makedirs(params.binary)
        except FileExistsError:
            print("breakpoint detection directory already exists.")
        data = np.loadtxt(input.d_matrix_file, delimiter=',')
        n_cells = data.shape[0]
        n_bins = data.shape[1]

        freq = np.fft.fftfreq(n_bins, 1)
        df = 0.015
        gpl = np.exp(- ((freq-1/(2*params.window_size))/(2*df))**2)  # pos. frequencies
        gmn = np.exp(- ((freq+1/(2*params.window_size))/(2*df))**2)
        g = gpl + gmn

        sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build',
                          '/cluster/work/bewi/members/pedrof/sc-dna/scripts/', persistence=True, postfix=f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{params.postfix}")

        bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=0.1, bp_limit=data.shape[1], compute_sp=False, evaluate_peaks=False)

        filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
        bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=params.threshold, lr=filtered_lr)

        os.rename(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{params.postfix}" + "_segmented_regions.txt", output.segmented_regions)
        os.rename(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{params.postfix}" + "_segmented_region_sizes.txt", output.segmented_region_sizes)

        # # n_cells and n_bins are global
        # try:
        #     cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.d_matrix_file}", f"--n_bins={n_bins}",\
        #         f"--n_cells={n_cells}", f"--window_size={params.window_size}", f"--postfix={n_nodes}nodes_{wildcards.regions}regions_"
        #         f"{wildcards.reads}reads_{wildcards.rep_id}",\
        #         f"--verbosity=1", f"--threshold=-1", f"--bp_limit={params.bp_limit}",\
        #         f"--compute_lr=true"])
        # except subprocess.SubprocessError as e:
        #     print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        # else:
        #     print(f"subprocess out: {cmd_output}")
        #     print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")
        #
        # # Threshold was -1, so output is only lr_vec
        # lr = np.loadtxt(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_lr_vec.csv", delimiter=',')
        #
        # # Save unfiltered LR
        # os.path.join(f"{BP_OUTPUT}_{sim_prefix}_trees", str(n_nodes) + 'nodes_' + wildcards.regions +'regions_'+ wildcards.reads + 'reads', wildcards.rep_id +'_' + "unfiltered_lr_vec.csv")
        #
        # filtered_lr = filter_lr(lr.T, H=g)
        # if filtered_lr.shape[0] < filtered_lr.shape[1]:
        #     filtered_lr = filtered_lr.T # make sure it's bins by cells
        # temp_lr_file = f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_temp_lr.txt"
        # np.savetxt(temp_lr_file, filtered_lr, delimiter=',')
        #
        # try:
        #     cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.d_matrix_file}", f"--n_bins={n_bins}",\
        #         f"--n_cells={n_cells}", f"--window_size={params.window_size}", f"--postfix={n_nodes}nodes_{wildcards.regions}regions_"
        #         f"{wildcards.reads}reads_{wildcards.rep_id}",\
        #         f"--verbosity=1", f"--threshold={params.threshold}", f"--bp_limit={params.bp_limit}",\
        #         f"--compute_lr=false", f"lr_file={temp_lr_file}"])
        # except subprocess.SubprocessError as e:
        #     print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        # else:
        #     print(f"subprocess out: {cmd_output}")
        #     print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")
        #
        # os.remove(temp_lr_file)
        #
        # for filename in bp_output_file_exts:
        #     os.rename(f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{filename}", os.path.join(f"{BP_OUTPUT}_{sim_prefix}_trees", str(n_nodes) + 'nodes_' + wildcards.regions +'regions_'+ wildcards.reads + 'reads', wildcards.rep_id +'_' + filename))

rule segment_regions:
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' +\
         '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt"
    run:
        filtered_counts = np.loadtxt(input.d_matrix_file, delimiter=',')
        n_cells = filtered_counts.shape[0]
        region_sizes = np.loadtxt(input.segmented_region_sizes)
        n_regions = len(region_sizes)
        sum_region_sizes = np.sum(region_sizes)
        condensed_mat = np.zeros((n_cells, n_regions))

        print("segmenting the bins...")
        for i in tqdm(range(n_cells)):
            region_id = 0
            region_count = 0
            # import ipdb; ipdb.set_trace() # debugging starts here
            for j in range(filtered_counts.shape[1]):
                to_add = filtered_counts[i][j]
                condensed_mat[i][region_id] += to_add
                region_count += 1
                if region_count == region_sizes[region_id]:
                    region_id += 1
                    region_count = 0

        if not np.allclose(condensed_mat.sum(axis=1), filtered_counts.sum(axis=1)):
            raise AssertionError(
                "not all values of the sums before & after "
                "segmentation are close")

        print("saving the segmented regions...")
        np.savetxt(output.segmented_counts, condensed_mat, delimiter=",")
        np.savetxt(output.segmented_counts_shape, condensed_mat.shape)

# rule learn_trees:
#     params:
#         cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
#         cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
#         full_tree_n_iters = config["inference"]["full_trees"]["n_iters"],
#         posfix = f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}"
#     input:
#         segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
#          + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
#         segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
#          + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
#     output:
#         cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
#         cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree_cnvs.csv',
#         full_tree = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
#         full_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree_cnvs.csv'
#     threads: 10
#     run:
#         segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
#         segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
#         alpha = 1./segmented_counts.shape[1]
#
#         sci = SCICoNE(binaries_path, output_temp, persistence=False, n_cells=n_cells, n_bins=n_bins, postfix=params.posfix)
#
#         # Run cluster trees
#         sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=10, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters, max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha)
#
#         # Store best cluster tree
#         with open(output.cluster_tree, "w") as file:
#             for line in sci.best_cluster_tree.tree_str.splitlines():
#                 file.write(f"{line}\n")
#         np.savetxt(output.cluster_tree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')
#
#         # Run full trees starting from cluster tree
#         sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=10, cluster=False, full=True, full_tree_n_iters=params.full_tree_n_iters, max_tries=1, robustness_thr=0.5, alpha=alpha)
#
#         # Store best full tree
#         with open(output.full_tree, "w") as file:
#             for line in sci.best_full_tree.tree_str.splitlines():
#                 file.write(f"{line}\n")
#         np.savetxt(output.full_tree_inferred_cnvs, sci.best_full_tree.outputs['inferred_cnvs'], delimiter=',')

rule learn_cluster_trees:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        posfix = "ctree" + f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree_cnvs.csv'
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        alpha = 1./segmented_counts.shape[1]

        sci = SCICoNE(binaries_path, output_temp, persistence=False, n_cells=n_cells, n_bins=n_bins, postfix=params.posfix)

        # Run cluster trees
        sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=10, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters, max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha)

        # Store best cluster tree
        with open(output.cluster_tree, "w") as file:
            for line in sci.best_cluster_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.cluster_tree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')


rule learn_full_trees:
    params:
        full_tree_n_iters = config["inference"]["full_trees"]["n_iters"],
        posfix = "ftree" + f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        full_tree = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
        full_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree_cnvs.csv'
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        alpha = 1./segmented_counts.shape[1]

        tree = Tree(binaries_path + '/inference', "", postfix=params.posfix)
        tree.read_from_file(input.cluster_tree)
        sci.best_cluster_tree = tree

        sci = SCICoNE(binaries_path, output_temp, persistence=False, n_cells=n_cells, n_bins=n_bins, postfix=params.posfix)

        # Run full trees starting from cluster tree
        sci.learn_tree(segmented_counts, segmented_region_sizes, n_reps=10, cluster=False, full=True, full_tree_n_iters=params.full_tree_n_iters,
                        max_tries=1, robustness_thr=0.5, alpha=alpha)

        # Store best full tree
        with open(output.full_tree, "w") as file:
            for line in sci.best_full_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.full_tree_inferred_cnvs, sci.best_full_tree.outputs['inferred_cnvs'], delimiter=',')

rule plot_breakpoints_thresholds:
    params:
        n_bins = n_bins,
        n_nodes = n_nodes
    input:
        gt_files = expand(os.path.join(f"{SIM_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{{regions}}'+'regions_'+ '{{reads}}'+'reads', '{rep_id}' +'_' + 'ground_truth.txt')\
        , rep_id=[x for x in range(0,all_n_tps)]),
        bps_files = expand(os.path.join(f"{BP_OUTPUT}_{sim_prefix}", str(n_nodes) + 'nodes_' + '{{regions}}'+'regions_'+ '{{reads}}'+'reads', '{rep_id}' +'_' + 'all_bps_comparison.csv')\
        , rep_id=[x for x in range(0,all_n_tps)])
    output:
        mean_tprs = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'mean_tpr_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        mean_fprs = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'mean_fpr_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        sum_tps = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'sum_tps_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png'),
        sum_fps = os.path.join(f"{BP_OUTPUT}_{sim_prefix}", "plots", 'sum_fps_'+ f"{n_nodes}nodes_" + '{regions}regions_{reads}reads.png')
    run:
        all_bins = range(0,params.n_bins)
        all_tpr, all_fpr, all_n_tps, all_n_fps, all_bps_tables = ([] for i in range(5))

        threshold_coeffs = np.linspace(1.0,16.0, 50) # 50 thresholds btw 1 and 16

        for gt_file, bps_file in tqdm(zip(input.gt_files, input.bps_files)):
            bps = pd.read_csv(bps_file, header=None)
            bps.columns = ['idx','log_sp','stdev']
            bps['ranking'] = bps['log_sp'] / bps['stdev']
            # bps = bps.sort_values('ranking',ascending=False)
            bps = bps.dropna()

            all_bps_tables.append(bps)

            # get the ground truth
            cell_genotypes = pd.read_csv(gt_file, sep=',' ,header=None)
            cell_genotypes = cell_genotypes[cell_genotypes.columns[:-1]] # remove the last (only NaN) column
            cell_bps = cell_genotypes.diff(periods=1, axis=1) # apply diff to detect breakpoints
            cell_bps = cell_bps.fillna(value=0.0) # diff makes the 1st row NaN, make it zero
            cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
            grouped_cell_bps = cell_bps.sum(axis=0) # count the non-zeroes
            ground_truth = grouped_cell_bps[grouped_cell_bps > 0] # if they occur in at least 1 cell
            ground_truth = ground_truth.index.tolist()
            # end of ground truth

            # correcting for the bps 1-2 bins nearby
            for index, row in bps.iterrows():
                idx_val = bps.loc[index, 'idx']
                for gt in ground_truth:
                    if (abs(idx_val - gt) <=2 and idx_val != gt):
                        print('correcting ' + str(idx_val) + '->' + str(gt))
                        bps.loc[index,'idx'] = gt

            tpr_values, fpr_values, n_tps, n_fps = [], [], [], []
            for thr in threshold_coeffs:
                predicted_positives = []
                predicted_negatives = []
                for index, row in bps.iterrows():
                    if row['ranking'] > thr:
                        predicted_positives.append(row['idx'])
                    else:
                        break

                #import ipdb; ipdb.set_trace()
                predicted_negatives = [i for i in all_bins if i not in predicted_positives]

                true_positives = [i for i in predicted_positives if i in ground_truth]
                false_positives = [i for i in predicted_positives if i not in ground_truth]

                true_negatives = [i for i in predicted_negatives if i not in ground_truth]
                false_negatives = [i for i in predicted_negatives if i in ground_truth]

                # import ipdb; ipdb.set_trace()
                assert(len(ground_truth) == (len(true_positives) + len(false_negatives)))
                tpr = len(true_positives) / len(ground_truth) # len(ground_truth)
                fpr = len(false_positives) / (params.n_bins - len(ground_truth)) # (len(false_positives) + len(true_negatives))
                tpr_values.append(tpr)
                fpr_values.append(fpr)
                n_tps.append(len(true_positives))
                n_fps.append(len(false_positives))


            all_tpr.append(tpr_values)
            all_fpr.append(fpr_values)
            all_n_tps.append(n_tps)
            all_n_fps.append(n_fps)


        # Average TRP, FPR, number of predicted positives for each threshold value over all runs
        print("Plotting...")
        tpr_df = pd.DataFrame(all_tpr)
        tpr_df.columns = threshold_coeffs
        mean_tpr = tpr_df.mean()
        ax = sns.lineplot(x=mean_tpr.index,y=mean_tpr.values).set_title('Average TPR values')
        plt.savefig(output.mean_tprs)
        plt.clf()

        fpr_df = pd.DataFrame(all_fpr)
        fpr_df.columns = threshold_coeffs
        mean_fpr = fpr_df.mean()
        ax = sns.lineplot(x=mean_fpr.index,y=mean_fpr.values).set_title('Average FPR values')
        plt.savefig(output.mean_fprs)
        plt.clf()

        n_tps_df = pd.DataFrame(all_n_tps)
        n_tps_df.columns = threshold_coeffs
        sum_tps = n_tps_df.sum()
        ax = sns.lineplot(x=sum_tps.index,y=sum_tps.values).set_title('Sum TP values')
        plt.savefig(output.sum_tps)
        plt.clf()

        n_fps_df = pd.DataFrame(all_n_fps)
        n_fps_df.columns = threshold_coeffs
        sum_fps = n_fps_df.sum()
        ax = sns.lineplot(x=sum_fps.index,y=sum_fps.values).set_title('Sum FP values')
        plt.savefig(output.sum_fps)
        plt.clf()

rule hmm_copy_inference:
    params:
        script = config["hmm_copy"]["script"],
        n_nodes = n_nodes,
        scratch = config["hmm_copy"]["scratch"],
        mem = config["hmm_copy"]["mem"],
        time = config["hmm_copy"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_HMMcopy_inferred.txt
        hmmcopy_inferred_cnvs = HMMCOPY_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_hmm_copy_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat} {output.hmmcopy_inferred_cnvs}"

rule hierarchical_copy_inference:
    params:
        script = config["hclust"]["script"],
        n_nodes = n_nodes,
        scratch = config["hclust"]["scratch"],
        mem = config["hclust"]["mem"],
        time = config["hclust"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        hclust_inferred_cnvs = HCLUST_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_hclust_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat} {output.hclust_inferred_cnvs}"

rule phenograph_clustering:
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        clusters_phenograph_assignment = f'{PHENO_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv',
        phenograph_inferred_cnvs = PHENO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_phenograph_inferred.txt'
    threads:
        config["phenograph"]["threads"]
    run:
        data = np.loadtxt(input.d_mat, delimiter=",")
        n_cells = data.shape[0]
        print(f"n_cells: {str(n_cells)}")
        n_neighbours = int(n_cells / 10)
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        communities, graph, Q = phenograph.cluster(data=data, k=n_neighbours, n_jobs=threads, jaccard=True)

        print(f"Communities: {communities}")
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]

        norm_data = data / np.sum(data, axis=1)[:, np.newaxis] # normalize by lib size

        # Set each cell to the average read counts of that state
        inferred_states = np.empty(data.shape)
        community_ids = np.unique(communities)
        for id in community_ids:
            selected_cells = np.where(communities==id)[0]
            cluster_counts = np.sum(norm_data[selected_cells], axis=0) # total number of norm counts of each bin across cluster
            cluster_profile = 2*cluster_counts / np.median(cluster_counts)
            inferred_states[selected_cells] = cluster_profile[np.newaxis,:]

        inferred_states = np.round(inferred_states)

        np.savetxt(output.phenograph_inferred_cnvs, inferred_states, delimiter=",")

        communities_df.to_csv(output.clusters_phenograph_assignment, sep="\t", index=False)
