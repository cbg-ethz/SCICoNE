import argparse
import os
import numpy as np
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("input_data_file")
parser.add_argument("cluster_assignments_file")
args = parser.parse_args()

input_data_file = args.input_data_file
cluster_assignments_file = args.cluster_assignments_file

input_data = np.loadtxt(input_data_file, delimiter=',')
cluster_assignments = np.loadtxt(cluster_assignments_file, delimiter=',').astype(int)

N, P = input_data.shape

cluster_assignments_df = pd.DataFrame(cluster_assignments, columns=["cluster"])
cluster_assignments_df["cell_barcode"] = cluster_assignments_df.index
cluster_assignments_df = cluster_assignments_df[["cell_barcode", "cluster"]]
cluster_dict = dict((Counter(cluster_assignments)))
cluster_ids = sorted(list(cluster_dict))

# Compute average counts of each cluster
avg_data = np.empty(input_data.shape)
condensed_avg_data = np.empty((len(cluster_ids), P))
cluster_sizes = np.zeros((len(cluster_ids),))

# Offset -1 if there is one
if np.min(cluster_ids) == -1:
    cluster_assignments = np.array(cluster_assignments) + 1
    cluster_ids = np.array(cluster_ids) + 1

for id in cluster_ids:
    avg_data[np.where(cluster_assignments==id)[0]] = np.mean(input_data[np.where(cluster_assignments==id)[0], :], axis=0)
    condensed_avg_data[id] = avg_data[np.where(cluster_assignments==id)[0][0],:]
    cluster_sizes[id] = np.where(cluster_assignments==id)[0].shape[0]

print(f"Cluster sizes: {cluster_sizes}")

input_data_filename = os.path.splitext(input_data_file)[0]
out_file = input_data_filename + '_condensed.txt'
np.savetxt(out_file, condensed_avg_data, delimiter=",")

print(f"Saved the condensed clusters matrix into {out_file}.")

out_file = input_data_filename + '_cluster_sizes.txt'
np.savetxt(out_file, cluster_sizes, delimiter=",")

print(f"Saved the cluster sizes into {out_file}.")
