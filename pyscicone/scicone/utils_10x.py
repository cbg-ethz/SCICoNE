import scicone.utils as utils
import h5py
import numpy as np

def read_hdf5(h5f_path, bins_to_exclude=None):
    extracted_data = dict()
    with h5py.File(h5f_path) as h5f:
        filtered_res = extract_filtered_corrected_counts_matrix(h5f, bins_to_exclude=bins_to_exclude)
        extracted_data['filtered_counts'] = filtered_res['filtered_counts']
        extracted_data['excluded_bins'] = filtered_res['excluded_bins']
        extracted_data['filtered_chromosome_stops'] = extract_chromosome_stops(h5f, bins_to_exclude=extracted_data['excluded_bins'])

    return extracted_data

def merge_data_by_chromosome(h5f, key="normalized_counts"):
    n_cells = h5f["cell_barcodes"][:].shape[0]
    sorted_chromosome_list = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    matrix_list = []
    for chr in sorted_chromosome_list:
        matrix_list.append(
            h5f[key][chr][:][0:n_cells, :]
        )  # select only the cells, not cell groups

    merged_matrix = np.concatenate(matrix_list, axis=1)
    return merged_matrix


def extract_filtered_corrected_counts_matrix(h5f, bins_to_exclude=None):
    filtered_counts = merge_data_by_chromosome(h5f, key='normalized_counts')
    sorted_chromosomes = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    # Keep only single cells
    n_cells = h5f["cell_barcodes"].shape[0]
    filtered_counts = filtered_counts[:n_cells,:]

    # Exclude unmappable bins
    is_mappable = []
    for ch in sorted_chromosomes:
        is_mappable = np.concatenate(
            [is_mappable, h5f["genome_tracks"]["is_mappable"][ch][()]]
        )

    is_excluded = ~np.array(is_mappable, dtype=bool)
    excluded_bins = np.where(is_excluded)[0]
    if bins_to_exclude:
        bins_to_exclude = np.array(bins_to_exclude)
        excluded_bins = np.unique(np.concatenate((excluded_bins, bins_to_exclude),0))
        is_excluded[excluded_bins] = True

    filtered_counts = filtered_counts[:, ~is_excluded]

    return dict(filtered_counts=filtered_counts, excluded_bins=excluded_bins)

def extract_chromosome_stops(h5f, bins_to_exclude=None):
    sorted_chromosomes = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    chr_ends = np.cumsum(h5f["constants"]["num_bins_per_chrom"][()])

    chr_stops = dict()
    if bins_to_exclude is not None:
        bins_to_exclude = np.array(bins_to_exclude)
        for idx, pos in enumerate(chr_ends):
            chr_stops[sorted_chromosomes[idx]] = pos-1 - len(bins_to_exclude[np.where(bins_to_exclude < pos)[0]])
    else:
        for idx, pos in enumerate(chr_ends):
            chr_stops[sorted_chromosomes[idx]] = pos-1

    return chr_stops
