import os
import itertools
import pandas as pd
import numpy as np
from collections import Counter
from tqdm import tqdm as tqdm
import phenograph
from scgenpy.preprocessing.utils import *

"""
This workflow generates observations sampled from a tree and infers a tree from
them, using the simulated segmentation. It requires a configuration file
specifying simulation and inference parameters.
"""

# Simulation settings
sim_bin = config["simulation"]["binary"]
n_cells = config["simulation"]["n_cells"]
n_nodes = config["simulation"]["n_nodes"]
n_regions = config["simulation"]["n_regions"]
n_bins = config["simulation"]["n_bins"]
n_reads = config["simulation"]["n_reads"]
nu = config["simulation"]["nu"]
all_n_tps = config["simulation"]["n_reps"]
sim_prefix=config["simulation"]["prefix"]
SIM_OUTPUT= config["simulation"]["output"]
sim_output_file_exts = ['d_mat.txt','ground_truth.txt','region_sizes.txt', 'tree.txt']

# Tree inference settings
inf_bin = config["inference"]["binary"]
ploidy = config["inference"]["ploidy"]
verbosity = config["inference"]["verbosity"]
copy_number_limit = config["inference"]["copy_number_limit"]
seed = config["inference"]["seed"]
tree_rep = config["inference"]["cluster_trees"]["n_reps"]
TREES_OUTPUT = config["inference"]["output"]
trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt']

rule all:
    input:
        avg_counts_median_norm = expand(f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts_median.csv'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)]),

        best_cluster_tree = expand(f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)], tree_rep_id=[x for x in range(0,tree_rep)]),

        best_full_tree = expand(f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt'\
        ,regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,all_n_tps)], tree_rep_id=[x for x in range(0,tree_rep)]),



    run:
        print("DONE.")

rule run_sim:
    params:
        sim_bin = sim_bin,
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        all_n_tps = all_n_tps,
        nu = nu
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ground_truth.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        tree = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree.txt'
    shell:
        "{params.sim_bin} --n_regions {wildcards.regions} --n_reads {wildcards.reads} --n_cells {params.n_cells} --n_bins {params.n_bins} --n_nodes \
        {params.n_nodes} --verbosity 0 --ploidy 2 --postfix {wildcards.rep_id} --nu {params.nu} --seed {wildcards.rep_id}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.csv {output.d_mat}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_ground_truth.csv {output.ground_truth}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_region_sizes.txt {output.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_tree.txt {output.tree}"

rule segment_regions:
    input:
        d_matrix_file = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' +\
         '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        segmented_region_sizes = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt'
    output:
        segmented_counts = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt"
    run:
        counts = np.loadtxt(input.d_matrix_file, delimiter=',')
        n_cells = counts.shape[0]
        region_sizes = np.loadtxt(input.segmented_region_sizes)
        n_regions = len(region_sizes)
        sum_region_sizes = np.sum(region_sizes)
        condensed_mat = np.zeros((n_cells, n_regions))

        print("segmenting the bins...")
        for i in tqdm(range(n_cells)):
            region_id = 0
            region_count = 0
            # import ipdb; ipdb.set_trace() # debugging starts here
            for j in range(counts.shape[1]):
                to_add = counts[i][j]
                condensed_mat[i][region_id] += to_add
                region_count += 1
                if region_count == region_sizes[region_id]:
                    region_id += 1
                    region_count = 0

        if not np.allclose(condensed_mat.sum(axis=1), counts.sum(axis=1)):
            raise AssertionError(
                "not all values of the sums before & after "
                "segmentation are close")

        print("saving the segmented regions...")
        np.savetxt(output.segmented_counts, condensed_mat, delimiter=",")
        np.savetxt(output.segmented_counts_shape, condensed_mat.shape)

rule phenograph_clustering_trees:
    input:
        segmented_counts = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv"
    output:
        clusters_phenograph_assignment = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'
    run:
        normalised_regions = np.loadtxt(input.segmented_counts, delimiter=",")
        n_cells = normalised_regions.shape[0]
        print(f"n_cells: {str(n_cells)}")
        n_neighbours = int(n_cells / 10)
        print(f"n_neighbours to be used: {str(n_neighbours)}")
        communities, graph, Q = phenograph.cluster(data=normalised_regions, k=n_neighbours, n_jobs=threads, jaccard=True)

        print(f"Communities: {communities}")
        communities_df = pd.DataFrame(communities, columns=["cluster"])
        communities_df["cell_barcode"] = communities_df.index
        communities_df = communities_df[["cell_barcode", "cluster"]]

        communities_df.to_csv(output.clusters_phenograph_assignment, sep="\t", index=False)

rule create_averaged_region_matrix:
    """
    Creates the averaged regions by using cluster assignments
    """
    input:
        segmented_counts = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        clusters_phenograph_assignment = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_clusters_phenograph_assignment.tsv'
    output:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv'
    run:
        print("loading the segmented counts...")
        segmented_counts = np.loadtxt(input.segmented_counts, delimiter=',')
        print("normalising the regions...")

        phenograph_assignments = pd.read_csv(input.clusters_phenograph_assignment, sep='\t')
        communities = phenograph_assignments.cluster.values

        community_dict = dict((Counter(communities)))
        community_ids = sorted(list(community_dict))

        cluster_sizes = [v for (k,v) in sorted(community_dict.items())]
        print(f"cluster sizes: {cluster_sizes}")

        cells_by_cluster = []
        for cluster in community_ids:
            cells_by_cluster.append(segmented_counts[communities == cluster])

        avg_clusters = [m.mean(0) for m in cells_by_cluster]
        avg_clusters_df = pd.DataFrame(avg_clusters)
        print(f"shape of average clusters: {avg_clusters_df.shape}")

        replicated_df = pd.DataFrame(np.repeat(avg_clusters_df.values,cluster_sizes,axis=0))

        np.savetxt(output.avg_counts, replicated_df.values, delimiter=",")

rule median_normalise_averages:
    params:
        ploidy = config["inference"]["ploidy"]
    input:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv'
    output:
        avg_counts_median_norm = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts_median.csv'
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.avg_counts, header=None)
        for index, row in df.iterrows():
            row_median = row.median()
            df.loc[index] /= row_median # median normalise
            df.loc[index] *= params.ploidy # set the median to ploidy

        np.savetxt(output.avg_counts_median_norm, df.values, delimiter=',')

rule learn_empty_tree:
    params:
        binary = inf_bin,
        ploidy = ploidy,
        verbosity = verbosity,
        copy_number_limit = copy_number_limit,
        n_iters = config["inference"]["learn_nu"]["n_iters"],
        n_nodes = config["inference"]["learn_nu"]["n_nodes"],
        move_probs = config["inference"]["learn_nu"]["move_probs"],
        seed = seed,
        posfix = f"empty_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}rep"
    input:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv',
        segmented_counts_shape = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) +\
         'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt'
    output:
        empty_tree = f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_empty_tree.txt',
        empty_tree_inferred_cnvs = f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}__empty_tree_cnvs.csv'
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]

        move_probs_str = ",".join(str(p) for p in params.move_probs)

        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.avg_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--move_probs={move_probs_str}", f"--seed={params.seed}", f"--region_sizes_file={input.segmented_region_sizes}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.empty_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.empty_tree_inferred_cnvs)

        if params.verbosity > 0:
            debug_info_out_path = os.path.join(analysis_path, "tree_learning", "debug_info", "empty_tree")
            debug_info_with_ap = os.path.join(debug_info_out_path, analysis_prefix)
            if not os.path.exists(debug_info_out_path):
                os.makedirs(debug_info_out_path)
            os.rename(f"{params.posfix}_cell_node_ids.tsv", f"{debug_info_with_ap}__cell_node_ids.tsv")
            os.rename(f"{params.posfix}_cell_region_cnvs.csv", f"{debug_info_with_ap}__cell_region_cnvs.csv")
            os.rename(f"{params.posfix}_markov_chain.csv", f"{debug_info_with_ap}__markov_chain.csv")
            os.rename(f"{params.posfix}_rel_markov_chain.csv", f"{debug_info_with_ap}__rel_markov_chain.csv")
            os.rename(f"{params.posfix}_acceptance_ratio.csv", f"{debug_info_with_ap}__acceptance_ratio.csv")
            os.rename(f"{params.posfix}_gamma_values.csv", f"{debug_info_with_ap}__gamma_values.csv")

rule learn_cluster_trees:
    params:
        binary = inf_bin,
        ploidy = ploidy,
        verbosity = verbosity,
        copy_number_limit = copy_number_limit,
        n_iters = config["inference"]["cluster_trees"]["n_iters"],
        n_nodes = config["inference"]["cluster_trees"]["n_nodes"],
        move_probs = config["inference"]["cluster_trees"]["move_probs"],
        posfix = f"cluster_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}_{tree_rep_id}",
        alpha = config["inference"]["cluster_trees"]["alpha"]
    input:
        avg_counts = f'{TREES_OUTPUT}_phenograph/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_avg_counts.csv',
        segmented_counts_shape = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes) +\
         'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        empty_tree = ancient(f'{TREES_OUTPUT}_nu_on_avg/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_empty_tree.txt')
    output:
        cluster_tree = f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_cluster_tree.txt',
        cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_cluster_tree_cnvs.csv'
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]

        move_probs_str = ",".join(str(p) for p in params.move_probs)

        with open(input.empty_tree) as file:
            for l in file:
                l_parts = l.split(':')
                if l_parts[0] == 'Nu':
                    nu = l_parts[1].strip()
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.avg_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--tree_file={input.empty_tree}", f"--move_probs={move_probs_str}", f"--seed={wildcards.tree_rep_id}",\
                f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}", f"--alpha={params.alpha}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.cluster_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.cluster_tree_inferred_cnvs)

        if params.verbosity > 0:
            debug_info_out_path = os.path.join(analysis_path, "tree_learning", "debug_info", "cluster_tree")
            debug_info_with_ap = os.path.join(debug_info_out_path, analysis_prefix)
            if not os.path.exists(debug_info_out_path):
                os.makedirs(debug_info_out_path)
            os.rename(f"{params.posfix}_cell_node_ids.tsv", f"{debug_info_with_ap}__{wildcards.tree_rep}_cell_node_ids.tsv")
            os.rename(f"{params.posfix}_cell_region_cnvs.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_cell_region_cnvs.csv")
            os.rename(f"{params.posfix}_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_markov_chain.csv")
            os.rename(f"{params.posfix}_rel_markov_chain.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_rel_markov_chain.csv")
            os.rename(f"{params.posfix}_acceptance_ratio.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_acceptance_ratio.csv")
            os.rename(f"{params.posfix}_gamma_values.csv", f"{debug_info_with_ap}__{wildcards.tree_rep}_gamma_values.csv")

rule pick_best_cluster_tree:
    input:
        cluster_trees = ancient(expand(f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{{regions}}regions_{{reads}}reads/{{rep_id}}_{tree_rep_id}_cluster_tree.txt'\
        , tree_rep_id=[x for x in range(0,tree_rep)])),
        cluster_tree_inferred_cnvs = ancient(expand(f'{TREES_OUTPUT}_cluster_tree/{str(n_nodes)}nodes_' + '{{regions}}regions_{{reads}}reads/{{rep_id}}_{tree_rep_id}_cluster_tree_cnvs.csv'\
        , tree_rep_id=[x for x in range(0,tree_rep)]))
    output:
        best_tree = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        best_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree_cnvs.csv'
    run:
        import operator

        trees_sorted = sorted(input.cluster_trees)
        trees_inferred_cnvs_sorted = sorted(input.cluster_tree_inferred_cnvs)

        tree_scores = get_tree_scores(trees_sorted)
        max_index, max_score = max(enumerate(tree_scores), key=operator.itemgetter(1))

        os.symlink(trees_sorted[max_index], output.best_tree)
        os.symlink(trees_inferred_cnvs_sorted[max_index], output.best_tree_inferred_cnvs)

rule pick_best_full_tree:
    input:
        full_trees = ancient(expand(f'{TREES_OUTPUT}_full_tree/{str(n_nodes)}nodes_' + '{{regions}}regions_{{reads}}reads/{{rep_id}}_{tree_rep_id}_full_tree.txt'\
        , tree_rep_id=[x for x in range(0,tree_rep)])),
        full_tree_inferred_cnvs = ancient(expand(f'{TREES_OUTPUT}_full_tree/{str(n_nodes)}nodes_' + '{{regions}}regions_{{reads}}reads/{{rep_id}}_{tree_rep_id}_full_tree_cnvs.csv'\
        , tree_rep_id=[x for x in range(0,tree_rep)]))
    output:
        best_tree = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree.txt',
        best_tree_inferred_cnvs = f'{TREES_OUTPUT}_best_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_full_tree_cnvs.csv'
    run:
        import operator

        trees_sorted = sorted(input.full_trees)
        trees_inferred_cnvs_sorted = sorted(input.full_tree_inferred_cnvs)

        tree_scores = get_tree_scores(trees_sorted)
        max_index, max_score = max(enumerate(tree_scores), key=operator.itemgetter(1))

        os.symlink(trees_sorted[max_index], output.best_tree)
        os.symlink(trees_inferred_cnvs_sorted[max_index], output.best_tree_inferred_cnvs)

rule learn_nu_on_cluster_tree:
    params:
        binary = inf_bin,
        ploidy = ploidy,
        verbosity = verbosity,
        copy_number_limit = copy_number_limit,
        n_iters = config["inference"]["learn_nu_cluster_trees"]["n_iters"],
        move_probs = config["inference"]["learn_nu_cluster_trees"]["move_probs"],
        n_nodes = 0, # needed only for the output naming
        seed = seed,
        posfix = f"nu_on_cluster_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}rep"
    input:
        best_tree = f'{TREES_OUTPUT}_best_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_cluster_tree.txt',
        segmented_counts = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt'
    output:
        nu_on_cluster_tree = f'{TREES_OUTPUT}_nu_on_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_nu_on_cluster_tree.txt',
        nu_on_cluster_tree_inferred_cnvs = f'{TREES_OUTPUT}_nu_on_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}__nu_on_cluster_tree.csv'
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]

        move_probs_str = ",".join(str(p) for p in params.move_probs)

        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.segmented_counts}", f"--tree_file={input.best_tree}",\
             f"--n_regions={n_regions}",f"--n_nodes={params.n_nodes}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}",\
                f"--move_probs={move_probs_str}", f"--seed={params.seed}", f"--region_sizes_file={input.segmented_region_sizes}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.nu_on_cluster_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.nu_on_cluster_tree_inferred_cnvs)

rule learn_full_trees:
    params:
        binary = inf_bin,
        ploidy = ploidy,
        verbosity = verbosity,
        copy_number_limit = copy_number_limit,
        n_iters = config["inference"]["full_trees"]["n_iters"],
        n_nodes = config["inference"]["full_trees"]["n_nodes"],
        move_probs = config["inference"]["full_trees"]["move_probs"],
        cf = config["inference"]["full_trees"]["cluster_fraction"],
        posfix = f"full_tree_tree{str(n_nodes)}" + "nodes_{regions}regions_{reads}reads_{rep_id}_{tree_rep_id}",
    input:
        nu_on_cluster_tree = f'{TREES_OUTPUT}_nu_on_cluster_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_nu_on_cluster_tree.txt',
        segmented_counts = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + "_segmented_counts.csv",
        segmented_counts_shape = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' +  "_segmented_counts_shape.txt",
        segmented_region_sizes = SIM_OUTPUT + '_' + sim_prefix +'/'+ str(n_nodes)\
         + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt'
    output:
        full_tree = f'{TREES_OUTPUT}_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_full_tree.txt',
        full_tree_inferred_cnvs = f'{TREES_OUTPUT}_full_tree/{str(n_nodes)}nodes_' + '{regions}regions_{reads}reads/{rep_id}_{tree_rep_id}_full_tree_cnvs.csv'
    run:
        input_shape = np.loadtxt(input.segmented_counts_shape)
        (n_cells, n_regions) = [int(input_shape[i]) for i in range(2)]

        move_probs_str = ",".join(str(p) for p in params.move_probs)

        with open(input.nu_on_cluster_tree) as file:
            for l in file:
                l_parts = l.split(':')
                if l_parts[0] == 'Nu':
                    nu = l_parts[1].strip()
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.segmented_counts}", f"--n_regions={n_regions}",\
                f"--n_cells={n_cells}", f"--ploidy={params.ploidy}", f"--verbosity={params.verbosity}", f"--postfix={params.posfix}",\
                f"--copy_number_limit={params.copy_number_limit}", f"--n_iters={params.n_iters}", f"--n_nodes={params.n_nodes}",\
                f"--tree_file={input.nu_on_cluster_tree}", f"--cf={params.cf}"\
                f"--move_probs={move_probs_str}", f"--seed={wildcards.tree_rep_id}", f"--region_sizes_file={input.segmented_region_sizes}", f"--nu={nu}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            print(f"subprocess out: {cmd_output}")
            print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_tree_inferred.txt", output.full_tree)
        os.rename(f"{params.n_nodes}nodes_{n_regions}regions_{params.posfix}_inferred_cnvs.csv", output.full_tree_inferred_cnvs)
