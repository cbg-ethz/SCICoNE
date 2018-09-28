configfile: "config.json"


'''
parameters
'''

n_regions = [25,50,100]
n_reads = [10000, 100000, 1000000]
n_repetitions = 100

'''
rules
'''


rule all:
    input:
    params:
	    sim_bin = config["simulations_bin"]
    output:
    shell:
	    "{params.sim_bin}"


