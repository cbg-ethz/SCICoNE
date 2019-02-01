import h5py
import argparse
import pandas as pd
import numpy as np

def merge_chromosomes(h5):

    n_cells = h5['cell_barcodes'][:].shape[0]
    all_chromosomes = list(h5['normalized_counts'].keys())
    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(h5['normalized_counts'][chr][:][0:n_cells,:]) # select only the cells, not cell groups

    cell_all_chrs = np.concatenate(cnv_matrices, axis=1)
    return cell_all_chrs


parser = argparse.ArgumentParser()
parser.add_argument("-h5", "--hdf5", required=True, help="cellranger-dna hdf5 output")
parser.add_argument("-b", "--bins", required=True, help="list of 10xgenomics artifacts (always low quality) bins to exclude")
parser.add_argument("-o","--output_path",required=False, default="./", help="path to the output")
parser.add_argument("-s", "--sample_name",required=False, default="", help="name of the sample")

args = parser.parse_args()

h5f = h5py.File(args.hdf5, 'r')

mat = merge_chromosomes(h5f)

bin_size = h5f["constants"]["bin_size"][()]
n_bins = mat.shape[1]
bin_ids = [x for x in range(0,n_bins)]
bin_df = pd.DataFrame(bin_ids, columns=["bin_ids"])

bin_df["start"] = bin_df["bin_ids"] * bin_size
bin_df["end"] = bin_df["start"] + bin_size
print(bin_df.head())

# exclude 10x artifact bins
artifact_bins = np.loadtxt(args.bins, delimiter='\t').astype(bool)

assert(artifact_bins.shape[0] == mat.shape[1])

print("artifact bins mask len")
print(len(artifact_bins))
print("artifact bins mask sum")
print(sum(artifact_bins))

print("matrix shape before & after filtering")
print(mat.shape)
mat = mat[:,~artifact_bins]
print(mat.shape)

print("bin_df shape before & after filtering")
print(bin_df.shape)
bin_df = bin_df[~artifact_bins]
print(bin_df.shape)

np.savetxt(args.output_path + '/' + args.sample_name +"_filtered_counts.tsv", mat, delimiter='\t')

bin_df.to_csv(args.output_path + '/' + args.sample_name + "_bins_genome.tsv",sep='\t',index=False)

print("Output written to: " + args.output_path)

h5f.close()
