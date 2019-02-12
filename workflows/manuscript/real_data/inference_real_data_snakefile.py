
verbosity = config["verbosity"]
# methods config
n_nodes = config["cnv_trees"]["n_nodes"]
n_iters = config["cnv_trees"]["n_iterations"]
ploidy = config["cnv_trees"]["ploidy"]
# seed = config["cnv_trees"]["seed"]
cn_limit = config["cnv_trees"]["copy_number_limit"]

# dataset config
n_bins = config["dataset"]["n_bins"]
n_cells = config["dataset"]["n_cells"]
postfix = config["dataset"]["name"]
matrix = config["dataset"]["matrix_path"]

INFERENCE_OUTPUT = config["inference_output"]
BP_OUTPUT= config["bp_output"]

inference_prefix = config["inference_prefix"]

trees_inf_output_exts = ['tree_inferred.txt', 'inferred_cnvs.txt'] # , 'inferred_cnvs_segmented.txt', 'tree_inferred_segmented.txt']

try:
    n_inference_reps = config["cnv_trees"]["n_reps"]
except KeyError:
    n_inference_reps = 10


'''
rules
'''

rule all:
    input:
        inferences_with_rep = expand(INFERENCE_OUTPUT + '/' + inference_prefix +'/'+ str(n_nodes) + 'nodes_'  + '0regions_'+ '-1'+'reads'+ '/'+ \
                'infrep{rep_inf}'+'_' + '{output_ext}', output_ext=trees_inf_output_exts, \
                rep_inf=[x for x in range(0,n_inference_reps)])
    output:
    shell:
	    "echo STATUS:SUCCESS. All of the rules are ran through."

rule breakpoint_detection:
    # --n_bins 18175 --n_cells 260 --postfix vancouver --window_size 10 --verbosity 3 --threshold 3 --d_matrix_file /Users/mtuncel/git_repos/sc-dna/data/adi_steif/read_count_tables/SA501X3F_corr_amp.txt
    params:
        binary = config["breakpoint_detection"]["bin"],
        n_bins = n_bins,
        n_cells = n_cells,
        postfix = postfix,
        window_size = config["breakpoint_detection"]["window_size"],
        verbosity = verbosity,
        threshold = config["breakpoint_detection"]["threshold"],
        scratch = config["breakpoint_detection"]["scratch"],
        mem = config["breakpoint_detection"]["mem"],
        time = config["breakpoint_detection"]["time"]
    input:
        d_mat = matrix
    output:
        region_sizes = BP_OUTPUT + '/' + postfix + "_segmented_region_sizes.txt"
    shell:
        "{params.binary} --n_bins {params.n_bins} --n_cells {params.n_cells} --postfix {params.postfix} --d_matrix_file {input.d_mat};\
        mv {params.postfix}_segmented_region_sizes.txt {output.region_sizes}"

rule infer_trees:
    # --print_precision 15 --n_nodes 10 --n_bins 18175 --n_cells 260 --postfix vancouver --n_iters 50000 --verbosity 1 --ploidy 2 --d_matrix_file /Users/mtuncel/git_repos/sc-dna/data/adi_steif/read_count_tables/SA501X3F_corr_amp.txt --region_sizes_file /Users/mtuncel/git_repos/sc-dna/cmake-build-debug/vancouver_segmented_region_sizes.txt --seed 22 --copy_number_limit 3
    params:
        binary = config["cnv_trees"]["bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        n_iters = n_iters,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"],
        cn_limit = cn_limit
    input:
        d_mat = matrix,
        region_sizes = rules.breakpoint_detection.output.region_sizes
    output:
        inferred_cnvs = INFERENCE_OUTPUT+ '/' + inference_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + 'infrep{rep_inf}' + '_inferred_cnvs.txt',
        inferred_tree = INFERENCE_OUTPUT+ '/' + inference_prefix +'/'+ str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + 'infrep{rep_inf}' + '_tree_inferred.txt'
    shell:
        "{params.binary} --n_reads {wildcards.reads} --n_regions {wildcards.regions}  --copy_number_limit {params.cn_limit}  --n_nodes {params.n_nodes} --n_bins {params.n_bins} --n_iters {params.n_iters} --n_cells {params.n_cells} --verbosity 0 \
        --ploidy 2  --postfix infrep{wildcards.rep_inf} --d_matrix_file {input.d_mat} --region_sizes_file {input.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_infrep{wildcards.rep_inf}_tree_inferred.txt {output.inferred_tree}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_infrep{wildcards.rep_inf}_inferred_cnvs.txt {output.inferred_cnvs}"
