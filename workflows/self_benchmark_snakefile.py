import datetime

configfile: "self_benchmark_config.json"


'''
parameters
'''

n_regions = [25,50,100]
n_reads = [10000, 100000, 1000000]
n_repetitions = 100
n_iters = 1000000 # 1 million iters for each setting


SIM_OUTPUT= "simulations_output"

time = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")

'''
rules
'''
# sample input: 50regions_10000reads_1538315330937578_deltas.csv

rule all:
    input:
        read_region_sims = expand(SIM_OUTPUT+'/'+time+'_'+'{regions}'+'regions_'+ '{reads}'+'reads'+'_deltas.csv', regions=n_regions,reads=n_reads)
    output:
    shell:
	    "echo STATUS:SUCCESS. All of the rules are ran through."

rule run_sim:
    params:
        sim_bin = config["simulations_bin"],
        n_repetitions = n_repetitions,
        n_iters = n_iters,
        scratch = config["cnv_trees"]["scratch"],
        mem = config["cnv_trees"]["mem"],
        time = config["cnv_trees"]["time"]
    threads:
        config["cnv_trees"]["threads"]
    output:
        # files = rules.all.input.read_region_sims
        SIM_OUTPUT+'/'+time+'_'+'{regions}'+'regions_'+ '{reads}'+'reads'+'_deltas.csv'
    shell:
        "{params.sim_bin} --n_regions {wildcards.regions} --n_reads {wildcards.reads} --n_iters {params.n_iters} --n_rep {params.n_repetitions}; ls -tp {wildcards.regions}regions_{wildcards.reads}reads_*_deltas.csv | grep -v / | head -n1 | xargs -I '{{}}' mv '{{}}' {output}"



