configfile: "config.json"


'''
parameters
'''

n_regions = [25,50,100]
n_reads = [10000, 100000, 1000000]
n_repetitions = 100
n_iters = 1000000 # 1 million iters for each setting


SIM_OUTPUT= "simulations_output"


'''
rules
'''
# sample input: 50regions_10000reads_1538315330937578_deltas.csv

rule all:
    input:
        expand(SIM_OUTPUT+'/{regions}'+'regions_'+ '{reads}'+'reads'+'_deltas.csv', regions=n_regions,reads=n_reads)
    params:
        sim_bin = config["simulations_bin"]
    output:
    shell:
	    "echo STATUS:SUCCESS. All of the rules are ran through."

rule run_sim:
    input:
    output:
        files = SIM_OUTPUT+'/{regions}'+'regions_'+ '{reads}'+'reads'+'_deltas.csv'
    shell:
        "touch {output.files}"
