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

n_cells = h5f['cell_barcodes'].value.shape[0]
all_chromosomes = list(h5f['normalized_counts'].keys())

number_chromosomes = sorted([int(x) for x in sorted(all_chromosomes)[:-2]])
ordered_chromosomes = [str(x) for x in number_chromosomes] + sorted(all_chromosomes)[-2:]

chr_lengths = []
for ch in ordered_chromosomes:
    chr_lengths.append(h5f['normalized_counts'][ch][0:n_cells,:].shape[1])

chr_ends = np.cumsum(chr_lengths)

chr_start_positions = [None for x in range(0,chr_ends[-1])]
for idx, pos in enumerate(chr_ends):
    chr_start_positions[pos-1] = ordered_chromosomes[idx] # -1 because it is a size info


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

print("filtering chromosome start positions")
filtered_chr_starts = np.array(chr_start_positions)[~artifact_bins]

df_chr_starts = pd.DataFrame(columns=["chr"])
for idx,val in enumerate(filtered_chr_starts):
    if(val != None):
        #print((idx,val))
        df_chr_starts.loc[idx] = val

print("writing output...")

np.savetxt(args.output_path + '/' + args.sample_name +"_filtered_counts.tsv", mat, delimiter='\t')

bin_df.to_csv(args.output_path + '/' + args.sample_name + "_bins_genome.tsv",sep='\t',index=False)

df_chr_starts.to_csv(args.output_path + '/' + args.sample_name + "_chr_starts.tsv",sep='\t')

print("Output written to: " + args.output_path)

h5f.close()
