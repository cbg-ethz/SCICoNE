
n_rep = config['params']['n_rep']
sim_data = ['simdata_'+ str(x) for x in range(0,n_rep)]
SIM_OUTPUT = 'results/'


n_bins = config['params']['n_bins']
n_regions = config['params']['n_regions']
n_nodes =  config['params']['n_nodes']
n_reads =  config['params']['n_reads']
n_cells =  config['params']['n_cells']
window_size =  config['params']['window_size']

rule all:
    input:
        bps_comparison = expand(SIM_OUTPUT + '{dataset_id}' + '_all_bps_comparison.csv', dataset_id=sim_data),
        ground_truth = expand(SIM_OUTPUT + str(n_nodes) + 'nodes_' + str(n_regions) + 'regions_' + str(n_reads) + 'reads_' + '{dataset_id}' + '_effective_regions.txt', dataset_id=sim_data)
    output:
    shell:
        "echo Successfully generated all of the data."


rule simulate:
    params:
        binary = config["bin"]+'simulation',
        n_bins = n_bins,
        n_regions = n_regions,
        n_nodes = n_nodes,
        n_reads = n_reads,
        n_cells = n_cells,
        dmatrix_name = str(n_nodes) + 'nodes_' + str(n_regions) + 'regions_' + str(n_reads) + 'reads_' + '{dataset_id}_d_mat.txt',
        effective_regions_name = str(n_nodes) + 'nodes_' + str(n_regions) + 'regions_' + str(n_reads) + 'reads_' + '{dataset_id}_effective_regions.txt'
    input:
    output:
        dmatrix = SIM_OUTPUT + str(n_nodes) + 'nodes_' + str(n_regions) + 'regions_' + str(n_reads) + 'reads_' + '{dataset_id}_d_mat.txt',
        effective_regions = SIM_OUTPUT + str(n_nodes) + 'nodes_' + str(n_regions) + 'regions_' + str(n_reads) + 'reads_' + '{dataset_id}_effective_regions.txt'
        # 10nodes_40regions_10000reads_sim_tiny2_d_mat.txt
    shell:
        "{params.binary} --n_bins {params.n_bins} --n_regions {params.n_regions} --n_nodes {params.n_nodes} --n_reads {params.n_reads} \
        --n_cells {params.n_cells} --postfix {wildcards.dataset_id} --verbosity 2 --ploidy 2; mv {params.dmatrix_name} {output.dmatrix}; \
        mv {params.effective_regions_name} {output.effective_regions}"


rule bp_detect:
    params:
        binary = config["bin"]+'breakpoint_detection',
        n_bins = n_bins,
        n_cells = n_cells,
        window_size = window_size,
        bps_comparisons_name = '{dataset_id}_all_bps_comparison.csv'
    input:
        dmatrix = rules.simulate.output.dmatrix
    output:
        SIM_OUTPUT + '{dataset_id}_all_bps_comparison.csv'
    shell:
        "{params.binary} --n_bins {params.n_bins} --n_cells {params.n_cells} --window_size {params.window_size} \
                --threshold 1 --verbosity 3 --postfix {wildcards.dataset_id} --d_matrix_file {input.dmatrix}; \
                mv {params.bps_comparisons_name} {output}"


