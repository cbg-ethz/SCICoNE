from scicone.tree import Tree
from scicone import utils_10x, utils

import sys, os, shutil, subprocess
import multiprocessing
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
import phenograph
from collections import Counter

class SCICoNE(object):
    """
    This class  provides an interface to interact with the outputs from the C++
    program that infers a copy number phylogeny from single-cell WGS data.

    It exposes functions to visualize breakpoint detection results, annotate
    the bins, and plot the resulting trees. Its attributes can be further
    processed using scgenpy, a generic package for pre, post-processing and
    visualization of single-cell copy number data.
    """
    def __init__(self, binary_path=None, output_path='', persistence=False, postfix="PYSCICONETEMP", verbose=False):
        """Create a SCICoNE object.
        binary_path : type
            Path to SCICoNE binaries.
        output_path : type
            Path to SCICoNE output files.
        persistence : boolean
            Wether to delete output files from C++ after loading them into the class.
        """
        if binary_path is None:
            bpath = os.path.join(os.path.dirname(__file__), 'bin')
            try:
                assert os.path.isdir(bpath)
            except AssertionError:
                print("Could not find binaries, tried: {}".format(bpath), flush=True)

            # Determine if we're using Linux or Mac
            if sys.platform.startswith("linux"):
                simulation_binary = "linux-simulation"
                bp_binary = "linux-breakpoint_detection"
                inference_binary = "linux-inference"
                score_binary = "linux-score"
                tests_binary = "linux-tests"
            elif sys.platform == "darwin":
                simulation_binary = "simulation"
                bp_binary = "breakpoint_detection"
                inference_binary = "inference"
                score_binary = "score"
                tests_binary = "tests"
            else:
                raise RuntimeError("Operating system could not be determined or is not supported. "
                                   "sys.platform == {}".format(sys.platform), flush=True)
            # Prepend appropriate path separator
            self.simulation_binary = os.path.join(bpath, simulation_binary)
            self.bp_binary = os.path.join(bpath, bp_binary)
            self.inference_binary = os.path.join(bpath, inference_binary)
            self.score_binary = os.path.join(bpath, score_binary)
            self.tests_binary = os.path.join(bpath, tests_binary)
        else:
            self.binary_path = binary_path
            self.simulation_binary = os.path.join(self.binary_path, 'simulation')
            self.bp_binary = os.path.join(self.binary_path, 'breakpoint_detection')
            self.inference_binary = os.path.join(self.binary_path, 'inference')
            self.score_binary = os.path.join(self.binary_path, 'score')
            self.tests_binary = os.path.join(self.binary_path, 'tests')

        self.data = dict()
        self.bps = dict()
        self.region_gene_map = None

        self.output_path = output_path
        self.persistence = persistence
        self.postfix = postfix

        self.clustering_score = 0.
        self.best_cluster_tree = None
        self.cluster_tree_robustness_score = 0.
        self.cluster_tree_list = []
        self.best_full_tree = None
        self.full_tree_robustness_score = 0.
        self.tree_list = []

        self.verbose = verbose

    def run_tests(self):
        try:
            cmd_output = subprocess.run([self.tests_binary])
        except subprocess.SubprocessError as e:
            print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

    def read_10x(self, h5f_path, bins_to_exclude=None, downsampling_factor=1):
        self.data = utils_10x.read_hdf5(h5f_path, bins_to_exclude=bins_to_exclude, downsampling_factor=downsampling_factor)

    def simulate_data(self, n_cells=200, n_nodes=5, n_bins=1000, n_regions=40, n_reads=10000, nu=1.0, min_reg_size=10, max_regions_per_node=1, ploidy=2, region_neutral_states=None, verbosity=1, seed=42):
        if verbosity < 1:
            verbosity = 1

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
                cmd = [self.simulation_binary, f"--n_cells={n_cells}", f"--n_nodes={n_nodes}",\
                    f"--n_regions={n_regions}", f"--n_bins={n_bins}", f"--n_reads={n_reads}", f"--nu={nu}",\
                    f"--min_reg_size={min_reg_size}", f"--max_regions_per_node={max_regions_per_node}",\
                    f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={self.postfix}",\
                    f"--region_neutral_states_file={region_neutral_states_file}", f"--seed={seed}"]
                if self.verbose:
                    print(' '.join(cmd))
                cmd_output = subprocess.run(cmd)
            except subprocess.SubprocessError as e:
                print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            # The binary generates 4 files: *_d_mat.csv, *_ground_truth.csv, *_region_sizes.txt and *_tree.txt.
            # We read each one of these files into np.array, np.array, np.array and Tree structures, respectively,
            # and output a dictionary with keys "d_mat", "ground_truth", "region_sizes" and "tree".
            d_mat_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv"
            ground_truth_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv"
            region_sizes_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt"
            tree_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt"
            cell_node_ids_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_cell_node_ids.tsv"

            try:
                # Move to output directory
                cwd = os.getcwd()
                os.rename(os.path.join(cwd, d_mat_file), os.path.join(output_path, d_mat_file))
                os.rename(os.path.join(cwd, ground_truth_file), os.path.join(output_path, ground_truth_file))
                os.rename(os.path.join(cwd, region_sizes_file), os.path.join(output_path, region_sizes_file))
                os.rename(os.path.join(cwd, tree_file), os.path.join(output_path, tree_file))
                os.rename(os.path.join(cwd, cell_node_ids_file), os.path.join(output_path, cell_node_ids_file))
                done = True
            except OSError as e:
                pass

        d_mat_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv")
        ground_truth_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv")
        region_sizes_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt")
        tree_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt")
        cell_node_ids_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_cell_node_ids.tsv")
        try:
            output = dict()
            output['d_mat'] = np.loadtxt(d_mat_file, delimiter=',')
            output['ground_truth'] = np.loadtxt(ground_truth_file, delimiter=',')
            output['region_sizes'] = np.loadtxt(region_sizes_file, delimiter=',')

            tree = Tree(self.inference_binary, self.output_path)
            tree.outputs['cell_node_ids'] = np.loadtxt(cell_node_ids_file, delimiter='\t')
            tree.outputs['inferred_cnvs'] = output['ground_truth']
            tree.outputs['region_sizes'] = output['region_sizes']
            if region_neutral_states is not None:
                tree.outputs['region_neutral_states'] = region_neutral_states
                os.remove(region_neutral_states_file)
            else:
                tree.outputs['region_neutral_states'] = np.ones((n_regions,)) * ploidy
            tree.read_from_file(tree_file)
            output['tree'] = tree

        except OSError as e:
            pass
            # print("OSError: ", e.output, e.stdout, e.stderr)

        # Now that we've read all outputs into memory, we delete the temporary files if persistence==False
        if not self.persistence:
            shutil.rmtree(output_path)

        return output

    def detect_breakpoints(self, data=None, window_size=30, threshold=3.0, bp_limit=300, lr=None, sp=None,
                            evaluate_peaks=True, compute_lr=True, compute_sp=True, input_breakpoints=None, verbosity=1):
        if data is None:
            data = self.data['filtered_counts']

        n_cells = data.shape[0]
        n_bins = data.shape[1]
        verbosity = 1 if verbosity < 1 else verbosity # > 0 to generate all files

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
            input_breakpoints_new = []
            if input_breakpoints[0] != 0:
                input_breakpoints_new.append(0)
            for i in range(0, len(input_breakpoints)):
                input_breakpoints_new.append(input_breakpoints[i])
            if input_breakpoints_new[-1] != n_bins-1:
                input_breakpoints_new.append(n_bins-1)
            input_breakpoints_new = np.array(input_breakpoints_new)
            input_breakpoints_file = f"{postfix}_pre_input_breakpoints_file.txt"
            np.savetxt(input_breakpoints_file, input_breakpoints_new, delimiter=',')

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')
            # os.environ["OMP_NUM_THREADS"] = "4"
            cmd = [self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}",\
                f"--bp_limit={bp_limit}", f"--compute_lr={compute_lr}", f"--lr_file={lr_file}",\
                f"--compute_sp={compute_sp}", f"--sp_file={sp_file}", f"--verbosity={verbosity}",\
                f"--evaluate_peaks={evaluate_peaks}", f"--postfix={postfix}",\
                f"--input_breakpoints_file={input_breakpoints_file}"]
            if self.verbose:
                print(' '.join(cmd))
            cmd_output = subprocess.run(cmd)

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

        self.bps = output
        if 'unfiltered_chromosome_stops' in self.data.keys():
            print('Mapping to genes...')
            self.region_gene_map = utils.get_region_gene_map(self.data['bin_size'], self.data['unfiltered_chromosome_stops'],
                                    self.bps['segmented_regions'], self.data['excluded_bins'])
            print('Done.')

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
        n_neighbours = max(int(n_cells / 10), 2) # avoid errors
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        # Cluster the normalised segmented data
        normalised_segmented_data = segmented_data/np.sum(segmented_data, axis=1)[:,np.newaxis]
        communities, graph, Q = phenograph.cluster(data=normalised_segmented_data, k=n_neighbours, n_jobs=1, jaccard=True)
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]
        community_dict = dict((Counter(communities)))
        community_ids = sorted(list(community_dict))

        # Compute (unnormalised) average counts of each cluster
        avg_segmented_counts = np.empty(segmented_data.shape)
        condensed_avg_segmented_counts = np.empty((len(community_ids), n_regions))
        cluster_sizes = np.zeros((len(community_ids),))

        # Offset -1 if there is one
        if np.min(community_ids) == -1:
            communities = np.array(communities) + 1
            community_ids = np.array(community_ids) + 1

        for id in community_ids:
            avg_segmented_counts[np.where(communities==id)[0]] = np.mean(segmented_data[np.where(communities==id)[0], :], axis=0)
            condensed_avg_segmented_counts[id] = avg_segmented_counts[np.where(communities==id)[0][0],:]
            cluster_sizes[id] = np.where(communities==id)[0].shape[0]

        self.cluster_assignments = communities

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

    def learn_tree(self, data=None, segmented_region_sizes=None, n_reps=10, cluster=True, full=True, cluster_tree_n_iters=4000, nu_tree_n_iters=4000, full_tree_n_iters=4000, max_tries=2, robustness_thr=0.5, **kwargs):
        if segmented_region_sizes is None:
            segmented_region_sizes = self.bps['segmented_region_sizes']
        if data is None:
            data = self.data['filtered_counts']

        if data.shape[1] != segmented_region_sizes.shape[0]:
            print('Condensing regions...')
            segmented_data = self.condense_regions(data, segmented_region_sizes)
            print('Done.')
        else:
            segmented_data = data

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
            print('Learning cluster tree...')

            # Get the average read counts
            clustered_segmented_data, cluster_sizes, cluster_assignments, Q = self.condense_segmented_clusters(segmented_data)

            cnt = 0
            robustness_score = 0.
            tree = None
            while robustness_score < robustness_thr:
                if cnt >= max_tries:
                    break
                nu = tree.nu if tree is not None else 1.0
                tree, robustness_score, trees = self.learn_tree_parallel(clustered_segmented_data, segmented_region_sizes, n_reps=n_reps, nu=nu, cluster_sizes=cluster_sizes, initial_tree=tree, n_iters=cluster_tree_n_iters, verbose=self.verbose, **kwargs)
                cnt += 1

            print(f"Cluster tree finished with a robustness score of {robustness_score} after {cnt} tries")

            # Expand clusters back into cells to get the cell per bin genotype
            cluster_bin_genotypes = tree.outputs['inferred_cnvs']
            cluster_node_ids = tree.outputs['cell_node_ids']

            cell_bin_genotypes = np.empty((segmented_data.shape[0], cluster_bin_genotypes.shape[1]))
            cell_node_ids = np.empty((segmented_data.shape[0], 1))
            for id in np.unique(cluster_assignments):
                cell_bin_genotypes[np.where(cluster_assignments==id)[0]] = cluster_bin_genotypes[id]
                cell_node_ids[np.where(cluster_assignments==id)[0],0] = cluster_node_ids[id,1]

            # Update tree data to cell-level instead of cluster-level
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

            tree.read_tree_str(tree.tree_str)

            self.best_cluster_tree = tree
            self.cluster_tree_robustness_score = robustness_score
            self.clustering_score = Q
            self.cluster_tree_list = trees

            if self.region_gene_map is not None:
                self.best_cluster_tree.set_gene_event_dicts(self.region_gene_map)

        if full:
            if self.best_cluster_tree is not None:
                tree = self.best_cluster_tree

                print('Initializing nu for full tree.')
                # Update the nu on the full data (the nu on the clustered data is very different) with this tree
                nu = tree.nu
                tree = self.learn_single_tree(segmented_data, segmented_region_sizes, nu=nu, initial_tree=tree, n_iters=nu_tree_n_iters, move_probs=[0,0,0,0,0,0,0,0,0,0,0,1,0], postfix=f"nu_tree_{self.postfix}", verbose=self.verbose, **kwargs)
                print('Done. Will start from nu={}'.format(tree.nu))
                print('Learning full tree...')
                cnt = 0
                robustness_score = 0.
                while robustness_score < robustness_thr:
                    if cnt >= max_tries:
                        break
                    nu = tree.nu
                    tree, robustness_score, trees = self.learn_tree_parallel(segmented_data, segmented_region_sizes, n_reps=n_reps, nu=nu, initial_tree=tree, n_iters=full_tree_n_iters, verbose=self.verbose)
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

                if self.region_gene_map is not None:
                    self.best_full_tree.set_gene_event_dicts(self.region_gene_map)

                print(f"Full tree finished with a robustness score of {robustness_score} after {cnt} tries")
            else:
                raise Exception("Full trees require a cluster tree to start from. Please re-run with cluster=True.")

        if cluster and not full:
            return self.best_cluster_tree
        elif full:
            return self.best_full_tree
