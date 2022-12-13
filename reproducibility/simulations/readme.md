# Benchmarking SCICoNE
This directory contains a Snakemake workflow to benchmark SCICoNE against other single-cell copy number calling methods. The main workflow is defined in `benchmark.smk` and requires a configuration file containing the parameters used for all methods, `config.json`. After following the installation instructions, you can run this workflow on a HPC cluster using slurm with this command (with appropriately set `<your-path-to>` both in the command and in `config.json`):

```
sbatch -n 10 --time 119:00:00 --wrap="snakemake -s <your-path-to>/SCICoNE/reproducibility/workflows/benchmark.smk --configfile <your-path-to>/SCICoNE/reproducibility/workflows/config.json -j 500 -k --cluster 'sbatch -n 10 --time 23:57:00 --mem-per-cpu=10000' "
```

## Running SCICoNE independently
If you wish to run SCICoNE in Python but outside of the Snakemake workflow while keeping the same parameters as in the benchmarking, you can either copy the code from the  `detect_breakpoints`, `segment_regions`, `learn_cluster_trees` and `learn_full_trees` rules in sequence from `benchmark.smk`, or equivalently run the `run_scicone_recommended_params.py` script (this follows the same steps as the tutorial [tutorial](!https://github.com/cbg-ethz/SCICoNE/blob/master/docs/notebooks/tutorial.ipynb), but with the parameters we used in the full benchmarking).
