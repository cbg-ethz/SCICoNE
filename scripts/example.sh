#!/bin/bash
../build/simulation --n_cells=200 --n_nodes=5 --n_regions=10 --n_bins=1000 --n_reads=10000 --nu=10.0 --min_reg_size=10 --max_regions_per_node=1 --ploidy=2 --postfix=test --seed=42

../build/breakpoint_detection --d_matrix_file=5nodes_10regions_10000reads_test_d_mat.csv --n_bins=1000 --n_cells=200 --window_size=30 --threshold=6.0 --postfix=test

python ./segment_counts.py 5nodes_10regions_10000reads_test_d_mat.csv test_segmented_region_sizes.txt

# Full tree
../build/inference --d_matrix_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt --region_sizes_file=test_segmented_region_sizes.txt --n_regions=8 --n_cells=200 --ploidy=2 --verbosity=1 --copy_number_limit=2 --n_iters=4000 --seed=42 --postfix=test

python ./convert_to_graphviz.py test_tree_inferred.txt
dot -Tpdf test_tree_inferred.graphviz -o test_tree_inferred.pdf

# Cluster tree
python ./run_phenograph.py 5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt

python ./condense_clusters.py 5nodes_10regions_10000reads_test_d_mat_segmented_counts.txt 5nodes_10regions_10000reads_test_d_mat_segmented_counts_cluster_assignments.txt

../build/inference --d_matrix_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts_condensed.txt  --region_sizes_file=test_segmented_region_sizes.txt --cluster_sizes_file=5nodes_10regions_10000reads_test_d_mat_segmented_counts_cluster_sizes.txt --n_regions=8 --n_cells=6 --ploidy=2 --verbosity=1 --copy_number_limit=2 --n_iters=4000 --seed=42 --postfix=cluster

python ./convert_to_graphviz.py cluster_tree_inferred.txt
dot -Tpdf cluster_tree_inferred.graphviz -o cluster_tree_inferred.pdf
