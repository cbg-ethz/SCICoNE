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
                    n_nodes=3,  seed=42, postfix="", initial_tree=None, nu=1.0, cluster_sizes=None, region_neutral_states=None, alpha=0., max_scoring=True, copy_number_limit=2,
                    c_penalise=10.0, verbosity=1):
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
                    f"--max_scoring={max_scoring}", f"--c_penalise={c_penalise}", f"--region_neutral_states_file={temp_region_neutral_states_file}"])
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
                    f"--c_penalise={c_penalise}", f"--region_neutral_states_file={temp_region_neutral_states_file}"])
            except subprocess.SubprocessError as e:
                print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
            else:
                pass
                # print(f"subprocess out: {cmd_output}")
                # print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.remove(temp_segmented_data_file)
        os.remove(temp_segmented_region_sizes_file)
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
        self.clustering_score = 0.
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

        postfix = self.postfix + "bp"

        lr_file = ""
        if lr is not None:
            compute_lr = False
            # Use provided lr values instead of computing new ones
            if lr.shape[0] < lr.shape[1]:
                lr = lr.T # make sure it's bins by cells
            lr_file = f"{postfix}_pre_lr_vec.txt"
            np.savetxt(lr_file, lr, delimiter=',')

        sp_file = ""
        if sp is not None:
            compute_sp = False
            sp_file = f"{postfix}_pre_sp_vec.txt"
            np.savetxt(sp_file, sp, delimiter=',')

        input_breakpoints_file = ""
        if input_breakpoints is not None:
            input_breakpoints_file = f"{postfix}_pre_input_breakpoints_file.txt"
            np.savetxt(input_breakpoints_file, input_breakpoints, delimiter=',')

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')

            cmd_output = subprocess.run([self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}",\
                f"--bp_limit={bp_limit}", f"--compute_lr={compute_lr}", f"--lr_file={lr_file}",\
                f"--compute_sp={compute_sp}", f"--sp_file={sp_file}", f"--verbosity={verbosity}",\
                f"--evaluate_peaks={evaluate_peaks}", f"--postfix={postfix}",\
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
                if postfix in fn:
                    key = fn.split(postfix)[1].split('_', 1)[1].split('.')[0]
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
        n_neighbours = max(1, int(n_cells / 10))
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

        return condensed_avg_segmented_counts, cluster_sizes, communities, Q

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
            clustered_segmented_data, cluster_sizes, cluster_assignments, Q = self.condense_segmented_clusters(segmented_data)

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
            cluster_node_ids = tree.outputs['cell_node_ids']

            cell_bin_genotypes = np.empty((segmented_data.shape[0], cluster_bin_genotypes.shape[1]))
            cell_node_ids = np.empty((segmented_data.shape[0], 2))
            cell_node_ids[:, 0] = range(segmented_data.shape[0])
            for id in np.unique(cluster_assignments):
                cell_bin_genotypes[np.where(cluster_assignments==id)[0]] = cluster_bin_genotypes[id]
                cell_node_ids[np.where(cluster_assignments==id)[0],1] = cluster_node_ids[id,1]

            tree.outputs['inferred_cnvs'] = cell_bin_genotypes
            tree.outputs['cell_node_ids'] = cell_node_ids

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
            self.clustering_score = Q

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
