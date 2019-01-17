import h5py
import argparse
import pandas as pd
import numpy as np

def merge_chromosomes(h5):
    
    n_cells = h5['cell_barcodes'].value.shape[0]
    all_chromosomes = list(h5['cnvs'].keys())
    # list of all cnv arrays
    cnv_matrices = []
    for chr in all_chromosomes:
        cnv_matrices.append(h5['cnvs'][chr].value[0:n_cells,:]) # select only the cells, not cell groups
        
    cell_all_chrs = np.concatenate(cnv_matrices, axis=1)
    return cell_all_chrs

def filter_negative_bins(mat):
    df_arr = pd.DataFrame(mat)
    df_arr[df_arr < 0] = None
    column_filter_mask = ~df_arr.isnull().any()
    return column_filter_mask


parser = argparse.ArgumentParser()
parser.add_argument("-h5", "--hdf5", required=True, help="cellranger-dna hdf5 output")
parser.add_argument("-b", "--bins", required=False, help="list of low quality bins to exclude")

args = parser.parse_args()

if args.bins == None:
    print("no bins are provided")
