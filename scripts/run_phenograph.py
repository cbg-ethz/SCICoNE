import argparse
import os
import numpy as np
import phenograph


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_data_file")
    parser.add_argument("--normalise", default=True, type=bool)
    parser.add_argument("--n_neighbours", type=int)
    args = parser.parse_args()

    input_data_file = args.input_data_file
    normalise = args.normalise
    n_neighbours = args.n_neighbours
    input_data = np.loadtxt(input_data_file, delimiter=',')

    # Cluster the segmented counts
    N, P = input_data.shape
    K = max(int(N / 10), 2) # avoid errors
    if n_neighbours:
        K = int(n_neighbours)
        
    print(f"n_neighbours to be used: {str(K)}")

    # Cluster the normalised segmented data
    if normalise:
        print("Will cluster normalised data.")
        input_data = input_data/np.sum(input_data, axis=1)[:,np.newaxis] * P
    communities, graph, Q = phenograph.cluster(data=input_data, k=K, n_jobs=1, jaccard=True, min_cluster_size=1, seed=42)

    print(f"Found {len(np.unique(communities))} clusters.")

    input_data_file = os.path.splitext(input_data_file)[0]
    out_file = input_data_file + '_cluster_assignments.txt'
    np.savetxt(out_file, communities, delimiter=",")

    print(f"Saved the cluster assignments into {out_file}.")
