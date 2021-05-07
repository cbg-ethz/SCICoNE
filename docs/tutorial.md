# Command line tutorial
Here is an example run from the command line on simulated data. Please see [reference.md](../scicone/reference.md) for details on the arguments to the SCICoNE binaries. This tutorial assumes that you have already installed the C++ program as described in the [README](../README.md). All the files generated from this run are available at [example_data/](example_data/).

__Step 1__: First we generate the data:
```
./build/simulation --n_cells=200 --n_nodes=5 --n_regions=10 --n_bins=1000 --n_reads=10000 --nu=10.0 --min_reg_size=10 --max_regions_per_node=1 --ploidy=2 --postfix=test --seed=42
```

__Step 2__: Then we detect potential breakpoints from the cells by bins counts matrix:
```
./build/breakpoint_detection --d_matrix_file=5nodes_10regions_10000reads_test_d_mat.csv --n_bins=1000 --n_cells=200 --window_size=30 --threshold=6.0 --postfix=test
```
This command creates files containing the region stops and their sizes.

__Step 3__: To proceed with copy number calling and tree inference, we need to use the detected breakpoints to produce a cells by regions counts matrix. We can use a utility Python script for this:
```
python ./scripts/segment_counts.py 5nodes_10regions_10000reads_test_d_mat.csv test_segmented_region_sizes.txt
```

__Step 4__: We can now use the output file to infer a copy number event tree:
```
./build/inference --d_matrix_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt --region_sizes_file=test_segmented_region_sizes.txt --n_regions=8 --n_cells=200 --ploidy=2 --copy_number_limit=2 --n_iters=4000 --seed=42 --postfix=test
```
This will create a file containing the tree in an in-house format. To visualise the tree, convert it to `graphviz` and generate the figure:
```
python ./scripts/convert_to_graphviz.py test_tree_inferred.txt
dot -Tpdf test_tree_inferred.graphviz -o test_tree_inferred.pdf
```

## Finding a tree on clustered data
Often it is a good idea to first cluster the data and find an initial tree describing those clusters. We start by clustering the segmented data obtained from __Step 3__. This can be done any way you like, but we provide a script to run PhenoGraph clustering (which you must install to before running the script):
```
python ./scripts/run_phenograph.py 5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt
```

This outputs a file containing the cluster assignments for each cell. We use this to condense the clusters and get their sizes:
```
python ./scripts/condense_clusters.py 5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt 5nodes_10regions_10000reads_test_d_mat_segmented_counts_cluster_assignments.txt
```

We are now ready to run the tree inference on the condensed data. The command is similar to the one in __Step 4__, except that now we use the `--cluster_sizes_file` parameter as well and we adjust `--n_cells` to the number of clusters:
```
./build/inference --d_matrix_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts_condensed.txt  --region_sizes_file=test_segmented_region_sizes.txt --cluster_sizes_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts_cluster_sizes.txt --n_regions=8 --n_cells=6 --ploidy=2 --copy_number_limit=2 --n_iters=4000 --seed=42 --postfix=cluster
```

We can again visualise this tree:
```
python ./scripts/convert_to_graphviz.py cluster_tree_inferred.txt
dot -Tpdf cluster_tree_inferred.graphviz -o cluster_tree_inferred.pdf
```
