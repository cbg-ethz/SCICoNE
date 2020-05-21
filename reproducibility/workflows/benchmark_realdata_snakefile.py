
# parameters
# --n_regions 314 --n_reads 38218  --print_precision 16 --n_nodes 10 --n_bins 18175 --n_iters 100 --n_cells 260 --verbosity 2 --ploidy 2 --seed 42 --postfix adi_real
# --d_matrix_file ../data/adi_steif/read_count_tables/SA501X3F_corr_amp.txt  --region_sizes_file ../data/adi_steif/read_count_tables/region_sizes_real.txt

n_nodes = config["cnv_trees"]["n_nodes"]
n_regions = 126
n_bins = 18175
n_cells = 260
n_reads = 38218
n_iters = [10000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000]  # config["cnv_trees"]["n_iters"] # 10000000
ploidy = [2]
n_repetitions = 10

SIM_OUTPUT = "/cluster/work/bewi/members/tuncel/data/dna/recomb_data/recomb_real_data" 
output_file_exts = ['inferred_cnvs.txt', 'tree_inferred.txt']
# sample input: 10nodes_314regions_38218reads_adi_real_inferred_cnvs.txt
rule all:
    input:
        read_region_sims = expand(SIM_OUTPUT + '_' + 'ploidy_{plo}' +'/'+'{iters}' + 'n_iters_'  + str(n_nodes) + 'nodes_'  + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/'+ '{rep_id}' +'_' +  \
                '{output_ext}', plo=ploidy, iters=n_iters,  output_ext=output_file_exts, regions=n_regions,reads=n_reads, rep_id=[str(x)+"_new_prior" for x in range(0,n_repetitions)])
    output:
    shell:
        "echo STATUS:SUCCESS. All of the rules are ran through."
rule infer_trees:
    params:
        binary = config["inference_bin"],
        n_nodes = n_nodes,
        n_bins = n_bins,
        n_cells = n_cells,
        n_regions = n_regions,
        n_reads = n_reads,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"]
    input:
        d_mat = config["d_mat"],
        region_sizes = config["regions_file"]
    output:
        inferred_cnvs = SIM_OUTPUT+ '_' + 'ploidy_{plo}' +'/' + '{iters}' + 'n_iters_' + str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_inferred_cnvs.txt',
        inferred_tree = SIM_OUTPUT+ '_' + 'ploidy_{plo}' +'/' + '{iters}' + 'n_iters_' + str(n_nodes) + 'nodes_' + '{regions}'+'regions_'+ '{reads}'+'reads'+ '/' + '{rep_id}' + '_tree_inferred.txt'
    shell:
        "{params.binary} --n_reads {wildcards.reads} --n_nodes {params.n_nodes} --n_regions {params.n_regions} --n_bins {params.n_bins} --n_iters {wildcards.iters} --n_cells \
        {params.n_cells} --verbosity 0 \
        --ploidy {wildcards.plo} --postfix {wildcards.rep_id} --d_matrix_file {input.d_mat} --region_sizes_file {input.region_sizes}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_tree_inferred.txt {output.inferred_tree}; \
        mv {params.n_nodes}nodes_{wildcards.regions}regions_{wildcards.reads}reads_{wildcards.rep_id}_inferred_cnvs.txt {output.inferred_cnvs}"
