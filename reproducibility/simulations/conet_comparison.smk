"""
Here we show that the SCICoNE results in the CONET paper are due to SCICoNE
being passed normalised instead of raw counts.

We assess two cases:
1. SCICoNE with w=20 and input = normalized data to recover their results
2. SCICoNE with w=4 and input = normalized data * 50 to obtain proper results
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm as tqdm
import random
from math import sqrt
import seaborn as sns

from conet.generative_model import CNSampler, EventTree, CountsGenerator, EventTreeGenerator
from scicone import SCICoNE, Tree

all_n_tps = config["simulate"]["n_reps"]
n_cells = config["simulate"]["n_cells"]
n_bins = config["simulate"]["n_bins"]
n_nodes = config["simulate"]["n_nodes"] # values: [10,20,30]
window_size = config["bp_detection"]["window_size"]

sim_prefix=config["simulate"]["prefix"]
binaries_path = config['binaries_path']
output_path = config["output"]

SIM_OUTPUT= os.path.join(config["output"], "simulation")
BP_OUTPUT = os.path.join(config["output"], "bp_detection")
TREES_OUTPUT = os.path.join(config["output"], "inference")
GINKGO_OUTPUT = os.path.join(config["output"], "ginkgo")
CONET_GINKGO_OUTPUT = os.path.join(config["output"], "conet_ginkgo")
CONET_GINKGO_TREE_OUTPUT = os.path.join(config["output"], "conet_ginkgo_tree_distances")

output_temp = config["inference"]["output_temp"]
try:
    os.mkdir(output_temp)
except:
    pass

# the default case for cluster_fraction variable
try:
    cf = config["inference"]["full_trees"]["cluster_fraction"]
except KeyError:
    cf = 1.0

rule all:
    input:
        deltas = os.path.join(output_path, 'deltas_conetsims.csv')
    run:
        print("rule all")

rule run_sim_conet:
    params:
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' + '/'+ '{rep_id}' + '_ground_truth.txt',
    run:
        import random
        from math import sqrt
        from conet.generative_model import CNSampler, EventTree, CountsGenerator, EventTreeGenerator

        loci = params.n_bins
        tree_size = params.n_nodes
        cells = params.n_cells
        random.seed(int(wildcards.rep_id))
        np.random.seed(int(wildcards.rep_id))

        # generate event tree and cell data, this might take a while
        cn_s = CNSampler.create_default_sampler()
        t_gen = EventTreeGenerator(cn_sampler=cn_s, tree_size=tree_size, no_loci=loci)
        tree: EventTree = t_gen.generate_random_tree()
        d_gen = CountsGenerator(cn_s, tree)
        counts, attachment, corrected_counts, brkp_matrix = d_gen.generate_data(loci, cells)
        conet_cc = np.transpose(np.array(corrected_counts)[:, 5:])
        counts = np.array(counts)

        np.savetxt(output.d_mat, conet_cc, delimiter=',')
        np.savetxt(output.ground_truth, counts, delimiter=',')

rule detect_breakpoints:
    params:
        binary = config["bp_detection"]["bin"],
        window_size = config["bp_detection"]["window_size"],
        verbosity = config["bp_detection"]["verbosity"],
        threshold = config["bp_detection"]["threshold"],
        bp_limit = config["bp_detection"]["bp_limit"],
        postfix = sim_prefix
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_d_mat.txt'
    output:
        segmented_regions = BP_OUTPUT + '_' + sim_prefix +'_trees' + '_{fairness}' + '/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_segmented_regions.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees' + '_{fairness}' + '/' + str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_segmented_region_sizes.txt',
    run:
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

        sci = SCICoNE(binaries_path, output_temp, persistence=False, postfix=f"{n_nodes}nodes_{wildcards.fairness}_{wildcards.rep_id}_{params.postfix}")

        window_size = params.window_size
        if wildcards.fairness == 'fair':
            # Assuming data was simulated with coverage of 100reads per bin
            data = data * 50
            window_size = 4

        bps = dict()
        if np.sum(data[0])/n_bins <= 10:
            bps = sci.detect_breakpoints(data, window_size=window_size, threshold=0.1, bp_limit=params.bp_limit, compute_sp=False, evaluate_peaks=False)
            filtered_lr = filter_lr(bps['lr_vec'].T, H=g)
            bps = sci.detect_breakpoints(data, window_size=window_size, bp_limit=params.bp_limit, threshold=params.threshold, lr=filtered_lr)
        else:
            bps = sci.detect_breakpoints(data, window_size=window_size, threshold=params.threshold, bp_limit=params.bp_limit)

        np.savetxt(output.segmented_regions, bps['segmented_regions'], delimiter=',')
        np.savetxt(output.segmented_region_sizes, bps['segmented_region_sizes'], delimiter=',')

rule segment_regions:
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' +\
         '/' + '{rep_id}' + '_d_mat.txt',
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' +  "_segmented_counts_shape.txt"
    run:
        filtered_counts = np.loadtxt(input.d_matrix_file, delimiter=',')
        if wildcards.fairness == 'fair':
            # Assuming data was simulated with coverage of 100reads per bin
            filtered_counts = filtered_counts * 50

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
        posfix = "ctree" + f"{str(n_nodes)}" + "nodes_{fairness}_{rep_id}",
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_cluster_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree_cnvs.csv'
    threads: 10
    run:
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        segmented_region_sizes = np.loadtxt(input.segmented_region_sizes, delimiter=',')
        cov = np.mean(np.sum(segmented_counts, axis=1) / np.sum(segmented_region_sizes))
        alpha = 1./segmented_counts.shape[1]
        gamma = 1./cov
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
        posfix = "ctree_sum" + f"{str(n_nodes)}" + "nodes_{fairness}_{rep_id}"
    input:
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}' + '/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        cluster_tree_sum = f'{TREES_OUTPUT}_best_cluster_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs_sum = f'{TREES_OUTPUT}_best_cluster_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree_cnvs.csv'
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
        posfix = "ftree" + f"{str(n_nodes)}" + "nodes_{fairness}_{rep_id}"
    input:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree.txt',
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}/'+ str(n_nodes)\
         + 'nodes' + '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}/'+ str(n_nodes)\
         + 'nodes' + '/'+ '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        full_tree = f'{TREES_OUTPUT}_best_full_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_full_tree.txt',
        full_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_full_tree_cnvs.csv'
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
        posfix = "ftreesum" + f"{str(n_nodes)}" + "nodes_{fairness}_{rep_id}"
    input:
        cluster_tree = f'{TREES_OUTPUT}_best_cluster_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree.txt',
        segmented_counts = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}/'+ str(n_nodes)\
         + 'nodes' + '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_region_sizes = BP_OUTPUT + '_' + sim_prefix +'_trees_{fairness}/'+ str(n_nodes)\
         + 'nodes' + '/' + '{rep_id}' + '_segmented_region_sizes.txt'
    output:
        full_tree = f'{TREES_OUTPUT}_best_full_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_full_tree.txt',
        full_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_full_tree_cnvs.csv'
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

rule ginkgo_inference:
    params:
        script = config["ginkgo"]["script"],
        n_nodes = n_nodes,
        scratch = config["ginkgo"]["scratch"],
        mem = config["ginkgo"]["mem"],
        time = config["ginkgo"]["time"],
        script_inp = str(n_nodes)+"nodes_{rep_id}"
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_d_mat.txt'
    output:
        ginkgo_inferred_cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_ginkgo_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat} {output.ginkgo_inferred_cnvs}"


rule conet_ginkgo_inference:
    params:
        script = config["conet"]["script"],
        n_nodes = n_nodes,
        scratch = config["conet"]["scratch"],
        mem = config["conet"]["mem"],
        time = config["conet"]["time"],
        script_inp = str(n_nodes)+"nodes_{rep_id}",
        intermediate_dir = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '/',
        bin_path = config["conet"]["bin_path"]
        em_iters = config["conet"]["em_iters"]
        pt_iters = config["conet"]["pt_iters"]
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_d_mat.txt',
        init_cnvs = GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/'+ '{rep_id}' + '_ginkgo_inferred.txt',
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_conet_inferred.txt
        conet_inferred_cnvs = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/'+ '{rep_id}' + '_conet_ginkgo_inferred.txt',
        conet_inferred_tree = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/'+ '{rep_id}' + '_conet_ginkgo_tree.txt',
        conet_inferred_attachments = CONET_GINKGO_OUTPUT +'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_conet_ginkgo_attachments.txt'
    shell:
        "python {params.script} --bin_path {params.bin_path} --counts {input.d_mat} --cnvs {input.init_cnvs} --intermediate_dir {params.intermediate_dir} --seed {wildcards.rep_id} --em_iters {params.em_iters} --pt_iters {params.pt_iters} --out_cnvs {output.conet_inferred_cnvs} --out_tree {output.conet_inferred_tree} --out_attachments {output.conet_inferred_attachments}"


rule compute_deltas:
    input:
        simulations = expand(f'{SIM_OUTPUT}_{sim_prefix}/{str(n_nodes)}nodes' + '/{rep_id}_ground_truth.txt',
                                        rep_id=[x for x in range(0,all_n_tps)]),
        best_full_tree = expand(f'{TREES_OUTPUT}_best_full_tree' + '/{fairness}/' + f'{str(n_nodes)}nodes' + '/{rep_id}_full_tree_cnvs.csv',
                         rep_id=[x for x in range(0,all_n_tps)], fairness=["fair", "unfair"]),

        best_cluster_tree = expand(f'{TREES_OUTPUT}_best_cluster_tree' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree_cnvs.csv',
                            rep_id=[x for x in range(0,all_n_tps)], fairness=["fair", "unfair"]),

        best_full_tree_sum = expand(f'{TREES_OUTPUT}_best_full_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_full_tree_cnvs.csv',
                        rep_id=[x for x in range(0,all_n_tps)], fairness=["fair", "unfair"]),

        best_cluster_tree_sum = expand(f'{TREES_OUTPUT}_best_cluster_tree_sum' + '/{fairness}' + f'/{str(n_nodes)}nodes' + '/{rep_id}_cluster_tree_cnvs.csv',
                            rep_id=[x for x in range(0,all_n_tps)], fairness=["fair", "unfair"]),

        ginkgo = expand(GINKGO_OUTPUT+'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_ginkgo_inferred.txt',
                                             rep_id=[x for x in range(0,all_n_tps)]),

        conet = expand(CONET_GINKGO_OUTPUT+'/'+ str(n_nodes) + 'nodes' + '/' + '{rep_id}' + '_conet_ginkgo_inferred.txt',
                                             rep_id=[x for x in range(0,all_n_tps)]),
    output:
        out_fname = os.path.join(output_path, 'deltas_conetsims.csv')
    run:
        rows = []
        rows.append('index,rep_id,method,delta')
        i = 0
        for true_cnvs in input.simulations:
            rep_id = true_cnvs.split('/')[-1].split('_')[0]
            gt = np.loadtxt(true_cnvs, delimiter=',')
            if gt.shape[0] == n_bins: # should be cells by bins
                gt = np.transpose(gt)

            i+=1
            method = 'diploid'
            inf = np.ones(gt.shape) * 2
            error = np.sqrt(np.mean((gt-inf)**2))
            rows.append(f'{i},{rep_id},{method},{error}')

            i += 1
            method = 'ginkgo'
            inf_cnvs = f'{GINKGO_OUTPUT}/{n_nodes}nodes/{rep_id}_ginkgo_inferred.txt'
            inf = np.loadtxt(inf_cnvs, delimiter=' ')
            error = np.sqrt(np.mean((gt-inf)**2))
            rows.append(f'{i},{rep_id},{method},{error}')

            i += 1
            method = 'conet_ginkgo'
            inf_cnvs = f'{CONET_GINKGO_OUTPUT}/{n_nodes}nodes/{rep_id}_conet_ginkgo_inferred.txt'
            inf = np.loadtxt(inf_cnvs, delimiter=',')
            error = np.sqrt(np.mean((gt-inf)**2))
            rows.append(f'{i},{rep_id},{method},{error}')

            for fairness in ['unfair', 'fair']:
                method = f'cluster_tree_{fairness}'
                inf_cnvs = f'{TREES_OUTPUT}_best_cluster_tree/{fairness}/{n_nodes}nodes/{rep_id}_cluster_tree_cnvs.csv'
                inf = np.loadtxt(inf_cnvs, delimiter=',')
                error = np.sqrt(np.mean((gt-inf)**2))
                rows.append(f'{i},{rep_id},{method},{error}')

                i += 1
                method = f'cluster_tree_sum_{fairness}'
                inf_cnvs = f'{TREES_OUTPUT}_best_cluster_tree_sum/{fairness}/{n_nodes}nodes/{rep_id}_cluster_tree_cnvs.csv'
                inf = np.loadtxt(inf_cnvs, delimiter=',')
                error = np.sqrt(np.mean((gt-inf)**2))
                rows.append(f'{i},{rep_id},{method},{error}')

                i += 1
                method = f'full_tree_{fairness}'
                inf_cnvs = f'{TREES_OUTPUT}_best_full_tree/{fairness}/{n_nodes}nodes/{rep_id}_full_tree_cnvs.csv'
                inf = np.loadtxt(inf_cnvs, delimiter=',')
                error = np.sqrt(np.mean((gt-inf)**2))
                rows.append(f'{i},{rep_id},{method},{error}')

                i += 1
                method = f'full_tree_sum_{fairness}'
                inf_cnvs = f'{TREES_OUTPUT}_best_full_tree_sum/{fairness}/{n_nodes}nodes/{rep_id}_full_tree_cnvs.csv'
                inf = np.loadtxt(inf_cnvs, delimiter=',')
                error = np.sqrt(np.mean((gt-inf)**2))
                rows.append(f'{i},{rep_id},{method},{error}')

        with open(output.out_fname, 'w') as f:
            f.write('\n'.join(rows))

rule plot_deltas:
    input:
        deltas_csv = os.path.join(output_path, 'deltas_conetsims.csv')
    output:
        deltas_png = os.path.join(output_path, 'deltas_conetsims.png')
    run:
        df = pd.read_csv(input.deltas_csv, index_col=0)
        plt.figure(figsize=(10,6))
        bp = sns.boxplot(data=df, x='method', y='delta', hue='method', dodge=False)
        plt.setp(bp.get_xticklabels(), rotation=45)
        plt.ylim([0,0.5])
        plt.savefig(output.deltas_png)
