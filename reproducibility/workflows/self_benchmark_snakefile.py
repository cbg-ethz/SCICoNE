import datetime

# configfile: "self_benchmark_config.json"

'''
parameters
'''

'''
Trees of size n=10, 20 and 40
n_regions = n, 2n and 4n
10,000 bins and 500 cells
1000, 10,000 and 100,000 reads per cell
'''

n_nodes = config["cnv_trees"]["n_nodes"] # values: [10,20,30]
n_regions = [n_nodes,2*n_nodes,4*n_nodes]
n_bins = 10000
n_reads = [10000, 30000, 100000] # add 300000

try:
    n_repetitions = config["simulation"]["n_reps"]
except KeyError:
    n_repetitions = 100 

try:
    n_inference_reps = config["cnv_trees"]["n_reps"]
except KeyError:
    n_inference_reps = 10

n_cells = 500
n_iters = config["cnv_trees"]["n_iterations"]  # int(1000000*n_nodes/10)

output_file_exts = ['d_mat.txt','ground_truth.txt','region_sizes.txt', 'tree.txt', 'inferred_cnvs.txt', 'tree_inferred.txt', 'HMMcopy_inferred.txt','inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt'] # , 'inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

INFERENCE_OUTPUT = config["inference_output"]
SIM_OUTPUT= config["sim_output"]

sim_prefix=config["sim_prefix"]
inference_prefix = config["inference_prefix"]

'''
rules
'''

rule all:
    input:
#        read_region_sims = expand(SIM_OUTPUT + '_' + prefix +'/'+ str(n_nodes) + 'nodes_'  + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/'+ '{rep_id}' +'_' + '{output_ext}', \
#                output_ext=output_file_exts, regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,n_repetitions)]),
        inferences_with_rep = expand(INFERENCE_OUTPUT + '_' + inference_prefix + '_' + '{iters}'  +'/'+ str(n_nodes) + 'nodes_'  + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/'+ \ 
                '{rep_id}_infrep{rep_inf}'+'_' + '{output_ext}', iters=n_iters, output_ext=trees_inf_output_exts, regions=n_regions,reads=n_reads, rep_id=[x for x in range(0,n_repetitions)], rep_inf=[x for x in range(0,n_inference_reps)])
    output:
    shell:
	    "echo STATUS:SUCCESS. All of the rules are ran through."

rule infer_trees:
    params:
        binary = config["inference_bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"]
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt'
    output:
        inferred_cnvs = INFERENCE_OUTPUT+ '_' + inference_prefix + '_' + '{iters}'  +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}_infrep{rep_inf}' + '_inferred_cnvs.txt',
        inferred_tree = INFERENCE_OUTPUT+ '_' + inference_prefix + '_' + '{iters}'  +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}_infrep{rep_inf}' + '_tree_inferred.txt'
    shell:
        "{params.binary} --n_reads {wildcards.reads} --n_regions {wildcards.regions} --n_nodes {params.n_nodes} --n_bins {params.n_bins} --n_iters {wildcards.iters} --n_cells {params.n_cells} --verbosity 0 \
        --ploidy 2  --postfix {wildcards.rep_id}_infrep{wildcards.rep_inf} --d_matrix_file {input.d_mat} --region_sizes_file {input.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_infrep{wildcards.rep_inf}_tree_inferred.txt {output.inferred_tree}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_infrep{wildcards.rep_inf}_inferred_cnvs.txt {output.inferred_cnvs}"

rule infer_trees_with_segmentation:
    params:
        binary = config["inference_bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"]
    input:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt'
    output:
        inferred_cnvs = INFERENCE_OUTPUT+ '_' + inference_prefix + '_' + '{iters}' +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}_infrep{rep_inf}' + '_inferred_cnvs_segmented.txt',
        inferred_tree = INFERENCE_OUTPUT+ '_' + inference_prefix + '_' + '{iters}' +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}_infrep{rep_inf}' + '_tree_inferred_segmented.txt'
    shell:
        "{params.binary} --n_regions {wildcards.regions}  --n_reads {wildcards.reads} --n_nodes {params.n_nodes} --n_bins {params.n_bins} --n_iters {wildcards.iters} --n_cells {params.n_cells} --verbosity 0 \
        --ploidy 2 --postfix {wildcards.rep_id}_infrep{wildcards.rep_inf} --d_matrix_file {input.d_mat}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_infrep{wildcards.rep_inf}_tree_inferred_segmented.txt {output.inferred_tree}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_infrep{wildcards.rep_inf}_inferred_cnvs_segmented.txt {output.inferred_cnvs}"
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
    threads:
        config["hmm_copy"]["threads"]
    output:
        # sample output: 10nodes_10regions_100000reads_sim1_HMMcopy_inferred.txt
        inferred_cnvs = INFERENCE_OUTPUT+ '_' + inference_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_HMMcopy_inferred.txt'
    shell:
        " Rscript {params.script} {input.d_mat}" 

rule run_sim:
    params:
        sim_bin = config["simulations_bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        n_repetitions = n_repetitions,
        n_iters = n_iters,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"]
    threads:
        config["cnv_trees"]["threads"]
    output:
        d_mat = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_d_mat.txt',
        ground_truth = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_ground_truth.txt',
        region_sizes = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_region_sizes.txt',
        tree = SIM_OUTPUT+ '_' + sim_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree.txt'
    shell:
        "{params.sim_bin} --n_regions {wildcards.regions} --n_reads {wildcards.reads} --n_iters {params.n_iters} --n_cells {params.n_cells} --n_bins {params.n_bins} --n_nodes \
        {params.n_nodes} --verbosity 0 --ploidy 2 --postfix {wildcards.rep_id}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_d_mat.txt {output.d_mat}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_ground_truth.txt {output.ground_truth}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_region_sizes.txt {output.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_tree.txt {output.tree}"

