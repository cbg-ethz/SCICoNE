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

from scicone import SCICoNE, Tree


all_n_tps = config["simulate"]["n_reps"]
n_inference_reps = config["inference"]["full_trees"]["n_reps"]
tree_rep = config["inference"]["cluster_trees"]["n_reps"]

n_nodes = config["simulate"]["n_nodes"]
n_regions = [n_nodes, 2*n_nodes, 4*n_nodes]
n_bins = config["simulate"]["n_bins"]
coverages = config["simulate"]["coverages"]
n_reads = (np.array(coverages, dtype=int) * n_bins).tolist()
nu = config["simulate"]["nu"]
window_size = config["bp_detection"]["window_size"]

n_cells = config["simulate"]["n_cells"]
n_iters = config["inference"]["full_trees"]["n_iters"]  # int(1000000*n_nodes/10)

SIM_OUTPUT= os.path.join(config["output"], "simulation")
BP_OUTPUT = os.path.join(config["output"], "bp_detection")
PHENO_OUTPUT = os.path.join(config["output"], "phenograph")
TREES_OUTPUT = os.path.join(config["output"], "inference")
HMMCOPY_OUTPUT = os.path.join(config["output"], "hmm_copy")
GINKGO_OUTPUT = os.path.join(config["output"], "ginkgo")
SCOPE_OUTPUT = os.path.join(config["output"], "scope")
HCLUST_OUTPUT = os.path.join(config["output"], "hclust")
MEDALT_OUTPUT = os.path.join(config["output"], "medalt")
MEDALT_TREE_OUTPUT = os.path.join(config["output"], "medalt_tree_distances")
CONET_OUTPUT = os.path.join(config["output"], "conet")
CONET_GINKGO_OUTPUT = os.path.join(config["output"], "conet_ginkgo")
CONET_TREE_OUTPUT = os.path.join(config["output"], "conet_tree_distances")
CONET_GINKGO_TREE_OUTPUT = os.path.join(config["output"], "conet_ginkgo_tree_distances")

try:
    other_cnv_methods = config["other_methods"]
except Exception as e:
    print(e)
    other_cnv_methods = ["phenograph", "hmm_copy", "ginkgo", "scope", "hclust", "conet"]

sim_prefix=config["simulate"]["prefix"]

output_temp = config["inference"]["output_temp"]
try:
    os.mkdir(output_temp)
except:
    pass

binaries_path = config['binaries_path']

try:
    is_overdispersed = config["simulate"]["is_overdispersed"]
except:
    is_overdispersed = 1

try:
    normalize = config["simulate"]["normalize"]
except:
    normalize = False

try:
    cf = config["inference"]["full_trees"]["cluster_fraction"]
except KeyError:
    cf = 1.0


rule all:
    input:
        simulations = expand(f'{SIM_OUTPUT}_{sim_prefix}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_d_mat.txt',
                                        regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),
        best_full_tree = expand(f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
                        regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)], tree_rep_id=[x for x in range(0,tree_rep)]),

        best_cluster_tree = expand(f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
                            regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        best_full_tree_sum = expand(f'{TREES_OUTPUT}_best_full_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
                        regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)], tree_rep_id=[x for x in range(0,tree_rep)]),

        best_cluster_tree_sum = expand(f'{TREES_OUTPUT}_best_cluster_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
                            regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)]),

        other_cnvs = expand(os.path.join(config["output"], '{method}') +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_{method}_inferred.txt',
                                regions=[x for x in n_regions], reads=[x for x in n_reads], rep_id=[x for x in range(0,all_n_tps)], method=other_cnv_methods),
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
        is_overdispersed = is_overdispersed,
        nu = nu,
        min_reg_size = window_size,
        max_regions_per_node = 10,
        normalize = normalize,
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ground_truth.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        tree = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree.txt'
    run:
        done = False
        while not done:
            if is_overdispersed:
                try:
                    cmd_output = subprocess.run([params.sim_bin, f"--n_cells={params.n_cells}", f"--n_nodes={params.n_nodes}",\
                        f"--n_regions={wildcards.regions}", f"--n_bins={params.n_bins}", f"--n_reads={wildcards.reads}", f"--nu={params.nu}",\
                        f"--min_reg_size={params.min_reg_size}", f"--max_regions_per_node={params.max_regions_per_node}",\
                        f"--ploidy=2", f"--verbosity=0", f"--postfix={wildcards.rep_id}", f"--seed={wildcards.rep_id}"])
                except subprocess.SubprocessError as e:
                    print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)
            else:
                try:
                    cmd_output = subprocess.run([params.sim_bin, f"--n_cells={params.n_cells}", f"--n_nodes={params.n_nodes}",\
                        f"--n_regions={wildcards.regions}", f"--n_bins={params.n_bins}", f"--n_reads={wildcards.reads}",\
                        f"--min_reg_size={params.min_reg_size}", f"--max_regions_per_node={params.max_regions_per_node}",\
                        f"--ploidy=2", f"--verbosity=0", f"--postfix={wildcards.rep_id}", f"--seed={wildcards.rep_id}"])
                except subprocess.SubprocessError as e:
                    print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            # Normalize?
            if params.normalize:
                mat = np.loadtxt(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.csv", delimiter=',')
                mat = (mat / np.mean(mat, axis=1)[:,np.newaxis]) * 2
                np.savetxt(f"{params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.csv", mat, delimiter=',')

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
    run:
        try:
            os.makedirs(params.binary)
        except FileExistsError:
            print("breakpoint detection directory already exists.")
        data = np.loadtxt(input.d_matrix_file, delimiter=',')
        n_cells = data.shape[0]
        n_bins = data.shape[1]

        # Optional filter for low coverages
        def filter_lr(lr_matrix, H=None):
            freq = np.fft.fftfreq(lr_matrix.shape[-1], 1) # 1 Hz sampling rate

            filtered_lr = np.empty(lr_matrix.shape)
            for c in range(lr_matrix.shape[0]):
                X = np.fft.fft(lr_matrix[c])
                Y = X * H
                y = np.fft.ifft(Y)
                filtered_lr[c] = y

            return filtered_lr

        freq = np.fft.fftfreq(n_bins, 1)
        df = 0.015
        gpl = np.exp(- ((freq-1/(2*params.window_size))/(2*df))**2)  # pos. frequencies
        gmn = np.exp(- ((freq+1/(2*params.window_size))/(2*df))**2)
        g = gpl + gmn

        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=f"{n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_{params.postfix}")

        bps = dict()
        if np.sum(data[0])/n_bins <= 10:
            bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=0.1, bp_limit=params.bp_limit, compute_sp=False, evaluate_peaks=False)
            filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
            bps = sci.detect_breakpoints(data, window_size=params.window_size, bp_limit=params.bp_limit, threshold=params.threshold, lr=filtered_lr)
        else:
            bps = sci.detect_breakpoints(data, window_size=params.window_size, threshold=params.threshold, bp_limit=params.bp_limit)

        np.savetxt(output.segmented_regions, bps['segmented_regions'], delimiter=',')
        np.savetxt(output.segmented_region_sizes, bps['segmented_region_sizes'], delimiter=',')

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

rule learn_cluster_trees:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        cluster_tree_n_reps = config["inference"]["cluster_trees"]["n_reps"],
        posfix = "ctree" + f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}",
        cov = lambda w: int(w.reads)/n_bins,
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
        gamma = 1./params.cov
        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=params.posfix + f"{sim_prefix}")

        # Run cluster trees
        sci.learn_tree(segmented_counts, segmented_region_sizes, verbosity=2, n_reps=params.cluster_tree_n_reps, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters, max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, max_scoring=True, gamma=gamma)

        # Store best cluster tree
        with open(output.cluster_tree, "w") as file:
            for line in sci.best_cluster_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.cluster_tree_inferred_cnvs, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')

rule learn_cluster_trees_sum:
    params:
        cluster_tree_n_iters = config["inference"]["cluster_trees"]["n_iters"],
        cluster_tree_n_tries = config["inference"]["cluster_trees"]["n_tries"],
        cluster_tree_n_reps = config["inference"]["cluster_trees"]["n_reps"],
        posfix = "ctree_sum" + f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        cluster_tree_sum = f'{TREES_OUTPUT}_best_cluster_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs_sum = f'{TREES_OUTPUT}_best_cluster_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree_cnvs.csv'
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        alpha = 1./segmented_counts.shape[1]

        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=params.posfix + f"{sim_prefix}")

        # Run cluster trees
        sci.learn_tree(segmented_counts, segmented_region_sizes, verbosity=2, cluster=True, full=False, cluster_tree_n_iters=params.cluster_tree_n_iters,
                                max_tries=params.cluster_tree_n_tries, robustness_thr=0.5, alpha=alpha, max_scoring=False, n_reps=params.cluster_tree_n_reps)

        # Store best cluster tree
        with open(output.cluster_tree_sum, "w") as file:
            for line in sci.best_cluster_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.cluster_tree_inferred_cnvs_sum, sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')

rule learn_full_trees:
    params:
        full_tree_n_iters = config["inference"]["full_trees"]["n_iters"],
        move_probs = config["inference"]["full_trees"]["move_probs"],
        cf = cf,
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

        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=params.posfix + f"{sim_prefix}")

        tree = Tree(binaries_path + '/inference', "", postfix=params.posfix)
        tree.outputs['region_sizes'] = segmented_region_sizes
        tree.outputs['region_neutral_states'] = np.ones((segmented_region_sizes.shape[0],)) * 2

        tree.read_from_file(input.cluster_tree)
        sci.best_cluster_tree = tree

        # Run full trees starting from cluster tree
        sci.learn_tree(segmented_counts, segmented_region_sizes, verbosity=2, n_reps=10, cluster=False, full=True, full_tree_n_iters=params.full_tree_n_iters,
                        max_tries=3, robustness_thr=0.5, alpha=alpha, move_probs=params.move_probs, max_scoring=True)

        # Store best full tree
        with open(output.full_tree, "w") as file:
            for line in sci.best_full_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.full_tree_inferred_cnvs, sci.best_full_tree.outputs['inferred_cnvs'], delimiter=',')

rule learn_full_trees_sum:
    params:
        full_tree_n_iters = config["inference"]["full_trees"]["n_iters"],
        move_probs = config["inference"]["full_trees"]["move_probs"],
        cf = cf,
        posfix = "ftreesum" + f"{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        full_tree = f'{TREES_OUTPUT}_best_full_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
        full_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree_sum/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree_cnvs.csv'
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        alpha = 1./segmented_counts.shape[1]

        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=params.posfix + f"{sim_prefix}")

        tree = Tree(binaries_path + '/inference', "", postfix=params.posfix)
        tree.outputs['region_sizes'] = segmented_region_sizes
        tree.outputs['region_neutral_states'] = np.ones((segmented_region_sizes.shape[0],)) * 2

        tree.read_from_file(input.cluster_tree)
        sci.best_cluster_tree = tree

        # Run full trees starting from cluster tree
        sci.learn_tree(segmented_counts, segmented_region_sizes, verbosity=2, n_reps=10, cluster=False, full=True, full_tree_n_iters=params.full_tree_n_iters,
                        max_tries=3, robustness_thr=0.5, alpha=alpha, move_probs=params.move_probs, max_scoring=False)

        # Store best full tree
        with open(output.full_tree, "w") as file:
            for line in sci.best_full_tree.tree_str.splitlines():
                file.write(f"{line}\n")
        np.savetxt(output.full_tree_inferred_cnvs, sci.best_full_tree.outputs['inferred_cnvs'], delimiter=',')

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

rule ginkgo_inference:
    params:
        script = config["ginkgo"]["script"],
        n_nodes = n_nodes,
        scratch = config["ginkgo"]["scratch"],
        mem = config["ginkgo"]["mem"],
        time = config["ginkgo"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_HMMcopy_inferred.txt
        ginkgo_inferred_cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ginkgo_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat} {output.ginkgo_inferred_cnvs}"

rule scope_inference:
    params:
        script = config["scope"]["script"],
        n_nodes = n_nodes,
        scratch = config["scope"]["scratch"],
        mem = config["scope"]["mem"],
        time = config["scope"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_scope_inferred.txt
        scope_inferred_cnvs = SCOPE_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_scope_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat} {output.scope_inferred_cnvs}"

rule conet_inference:
    params:
        script = config["conet"]["script"],
        n_nodes = n_nodes,
        scratch = config["conet"]["scratch"],
        mem = config["conet"]["mem"],
        time = config["conet"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}",
        intermediate_dir = CONET_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '/',
        bin_path = config["conet"]["bin_path"]
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        init_cnvs = HMMCOPY_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_hmm_copy_inferred.txt',
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_conet_inferred.txt
        conet_inferred_cnvs = CONET_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_inferred.txt',
        conet_inferred_tree = CONET_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_tree.txt',
        conet_inferred_attachments = CONET_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_attachments.txt'
    shell:
        "python {params.script} --bin_path {params.bin_path} --counts {input.d_mat} --cnvs {input.init_cnvs} --intermediate_dir {params.intermediate_dir} --seed {wildcards.rep_id} --out_cnvs {output.conet_inferred_cnvs} --out_tree {output.conet_inferred_tree} --out_attachments {output.conet_inferred_attachments}"

rule conet_ginkgo_inference:
    params:
        script = config["conet"]["script"],
        n_nodes = n_nodes,
        scratch = config["conet"]["scratch"],
        mem = config["conet"]["mem"],
        time = config["conet"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}",
        intermediate_dir = CONET_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '/',
        bin_path = config["conet"]["bin_path"]
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        init_cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ginkgo_inferred.txt',
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_conet_inferred.txt
        conet_inferred_cnvs = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_ginkgo_inferred.txt',
        conet_inferred_tree = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_ginkgo_tree.txt',
        conet_inferred_attachments = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_conet_ginkgo_attachments.txt'
    shell:
        "python {params.script} --bin_path {params.bin_path} --counts {input.d_mat} --cnvs {input.init_cnvs} --intermediate_dir {params.intermediate_dir} --seed {wildcards.rep_id} --out_cnvs {output.conet_inferred_cnvs} --out_tree {output.conet_inferred_tree} --out_attachments {output.conet_inferred_attachments}"

rule run_medalt:
    params:
        script = config["medalt"]["script"],
        n_nodes = n_nodes,
        scratch = config["medalt"]["scratch"],
        mem = config["medalt"]["mem"],
        time = config["medalt"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        cnvs = SCOPE_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_scope_inferred.txt'
    output:
        medalt_tree = f'{MEDALT_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_inferred.txt'
    shell:
        " python2 {params.script} -I {input.cnvs} -O {output.medalt_tree}"

rule run_medalt_ginkgo:
    params:
        script = config["medalt"]["script"],
        n_nodes = n_nodes,
        scratch = config["medalt"]["scratch"],
        mem = config["medalt"]["mem"],
        time = config["medalt"]["time"],
        script_inp = str(n_nodes)+"nodes_{regions}regions_{reads}reads_{rep_id}"
    input:
        cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ginkgo_inferred.txt'
    output:
        medalt_tree = f'{MEDALT_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_ginkgo_inferred.txt'
    shell:
        " python2 {params.script} -I {input.cnvs} -O {output.medalt_tree}"

rule compute_medalt_tree_distances:
    params:
        script = config["medalt_tree_distances"]["script"],
        n_nodes = n_nodes,
    input:
        cnvs = SCOPE_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_scope_inferred.txt',
        medalt_tree = f'{MEDALT_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_inferred.txt',
    output:
        medalt_tree_distance = f'{MEDALT_TREE_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_tree_distance.txt'
    shell:
        "python {params.script} {input.medalt_tree} {input.cnvs} {output.medalt_tree_distance}"

rule compute_medalt_ginkgo_tree_distances:
    params:
        script = config["medalt_tree_distances"]["script"],
        n_nodes = n_nodes,
    input:
        cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ginkgo_inferred.txt',
        medalt_tree = f'{MEDALT_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_ginkgo_inferred.txt',
    output:
        medalt_tree_distance = f'{MEDALT_TREE_OUTPUT}/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_medalt_ginkgo_tree_distance.txt'
    shell:
        "python {params.script} {input.medalt_tree} {input.cnvs} {output.medalt_tree_distance}"
