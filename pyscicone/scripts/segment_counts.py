import argparse
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("d_matrix_file", help="The cell by bins read counts matrix CSV file.")
parser.add_argument("segmented_region_sizes", help="The segmented_region_sizes.txt file output from breakpoint_detection.")
args = parser.parse_args()

d_matrix_file = args.d_matrix_file
segmented_region_sizes = args.segmented_region_sizes

bin_counts = np.loadtxt(d_matrix_file, delimiter=',')
n_cells = bin_counts.shape[0]
region_sizes = np.loadtxt(segmented_region_sizes)
n_regions = len(region_sizes)
sum_region_sizes = np.sum(region_sizes)
region_counts = np.zeros((n_cells, n_regions))

print("Segmenting the bins into regions...")
for i in range(n_cells):
    region_id = 0
    region_count = 0
    for j in range(bin_counts.shape[1]):
        to_add = bin_counts[i][j]
        region_counts[i][region_id] += to_add
        region_count += 1
        if region_count == region_sizes[region_id]:
            region_id += 1
            region_count = 0

if not np.allclose(region_counts.sum(axis=1), region_counts.sum(axis=1)):
    raise AssertionError("Not all values of the sums before & after segmentation are close")

d_matrix_filename = os.path.splitext(d_matrix_file)[0]
out_file = d_matrix_filename + '_segmented_counts.txt'
np.savetxt(out_file, region_counts, delimiter=",")

print(f"Saved the cells by regions counts matrix into {out_file}, with {n_cells} rows and {n_regions} columns.")
