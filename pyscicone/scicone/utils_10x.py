import scicone.utils as utils
import h5py
import numpy as np

DEFAULT_BIN_SIZE_KB=20 # the 10x Genomics setting

def read_hdf5(h5f_path, bins_to_exclude=None, downsampling_factor=1):
    extracted_data = dict()
    with h5py.File(h5f_path, 'r') as h5f:
        res = extract_corrected_counts_matrix(h5f, downsampling_factor=downsampling_factor, filter=False)
        extracted_data['unfiltered_counts'] = res['unfiltered_counts']
        extracted_data['unfiltered_chromosome_stops'] = extract_chromosome_stops(h5f, downsampling_factor=downsampling_factor)
        filtered_res = extract_corrected_counts_matrix(h5f, bins_to_exclude=bins_to_exclude, downsampling_factor=downsampling_factor, filter=True)
        extracted_data['filtered_counts'] = filtered_res['filtered_counts']
        extracted_data['excluded_bins'] = filtered_res['excluded_bins']
        res = extract_cnvs(h5f, bins_to_exclude=bins_to_exclude, downsampling_factor=downsampling_factor, filter=True)
        extracted_data['filtered_cnvs'] = res['filtered_cnvs']
        extracted_data['filtered_chromosome_stops'] = extract_chromosome_stops(h5f, bins_to_exclude=extracted_data['excluded_bins'], downsampling_factor=downsampling_factor)
        extracted_data['bin_size'] = DEFAULT_BIN_SIZE_KB*downsampling_factor*10**3

    return extracted_data

def merge_data_by_chromosome(h5f, key="normalized_counts", downsampling_factor=1, method='sum'):
    if method not in ['sum', 'median']:
        raise Exception('Method must be sum or median.')

    downsampling_factor = np.max([1, downsampling_factor])
    n_cells = h5f["cell_barcodes"][:].shape[0]
    sorted_chromosome_list = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    matrix_list = []
    for ch in sorted_chromosome_list:
        chr_matrix = []
        mat = h5f[key][ch][:][0:n_cells]
        n_bins = mat.shape[1]
        if downsampling_factor > 1:
            for j in range(0, n_bins, downsampling_factor):
                start = j
                end = j+downsampling_factor
                if end > n_bins:
                    end = n_bins
                    start = end-downsampling_factor
                if method=='sum':
                    chr_matrix.append(np.nansum(mat[:, start:end], axis=1).reshape(-1,1)) # ignore NaNs
                elif method=='median':
                    median = np.nanmedian(mat[:, start:end], axis=1).reshape(-1,1) # ignore NaNs
                    if key == 'cnvs':
                        idx = median > 2
                        median[idx] = np.floor(median[idx]).astype(int)
                        idx = median < 2
                        median[idx] = np.ceil(median[idx]).astype(int)
                    chr_matrix.append(median)

            chr_matrix = np.concatenate(chr_matrix, axis=1)
        else:
            chr_matrix = mat # select only the cells, not cell groups

        matrix_list.append(chr_matrix)

    merged_matrix = np.concatenate(matrix_list, axis=1)
    return merged_matrix


def extract_corrected_counts_matrix(h5f, bins_to_exclude=None, downsampling_factor=1, filter=True):
    downsampling_factor = np.max([1, downsampling_factor])
    unfiltered_counts = merge_data_by_chromosome(h5f, key='normalized_counts', downsampling_factor=downsampling_factor, method='sum')
    sorted_chromosomes = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    # Keep only single cells
    n_cells = h5f["cell_barcodes"].shape[0]
    filtered_counts = unfiltered_counts[:n_cells,:]

    if filter:
        # Exclude unmappable bins
        is_mappable = []
        for ch in sorted_chromosomes:
            chr_is_mappable = []
            vec = h5f["genome_tracks"]["is_mappable"][ch][()]
            n_bins = vec.size
            if downsampling_factor > 1:
                for j in range(0, n_bins, downsampling_factor):
                    start = j
                    end = j+downsampling_factor
                    if end > n_bins:
                        end = n_bins
                        start = end-downsampling_factor
                    bin_is_mappable = np.any(vec[start:end])
                    chr_is_mappable.append(bin_is_mappable)
            else:
                chr_is_mappable = vec

            is_mappable = np.concatenate([is_mappable, chr_is_mappable])

        is_excluded = ~np.array(is_mappable, dtype=bool)
        excluded_bins = np.where(is_excluded)[0]
        if bins_to_exclude is not None:
            bins_to_exclude = np.array(bins_to_exclude)
            excluded_bins = np.unique(np.concatenate((excluded_bins, bins_to_exclude),0))
            is_excluded[excluded_bins] = True

        filtered_counts = unfiltered_counts[:, ~is_excluded]

        return dict(filtered_counts=filtered_counts, excluded_bins=excluded_bins)
    else:
        return dict(unfiltered_counts=filtered_counts)

def extract_chromosome_stops(h5f, bins_to_exclude=None, downsampling_factor=1):
    downsampling_factor = np.max([1, downsampling_factor])
    sorted_chromosomes = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    if downsampling_factor > 1:
        chr_ends = np.cumsum([len(np.arange(0, n_bins, downsampling_factor)) for n_bins in h5f["constants"]["num_bins_per_chrom"]])
    else:
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

def extract_cnvs(h5f, bins_to_exclude=None, downsampling_factor=1, filter=True):
    downsampling_factor = np.max([1, downsampling_factor])
    unfiltered_cnvs = merge_data_by_chromosome(h5f, key='cnvs', downsampling_factor=downsampling_factor, method='median')
    sorted_chromosomes = utils.sort_chromosomes(h5f["constants"]["chroms"][()].astype(str))

    # Keep only single cells
    n_cells = h5f["cell_barcodes"].shape[0]
    filtered_cnvs = unfiltered_cnvs[:n_cells,:]

    if filter:
        # Exclude unmappable bins
        is_mappable = []
        for ch in sorted_chromosomes:
            chr_is_mappable = []
            vec = h5f["genome_tracks"]["is_mappable"][ch][()]
            n_bins = vec.size
            if downsampling_factor > 1:
                for j in range(0, n_bins, downsampling_factor):
                    start = j
                    end = j+downsampling_factor
                    if end > n_bins:
                        end = n_bins
                        start = end-downsampling_factor
                    bin_is_mappable = np.any(vec[start:end])
                    chr_is_mappable.append(bin_is_mappable)
            else:
                chr_is_mappable = vec

            is_mappable = np.concatenate([is_mappable, chr_is_mappable])

        is_excluded = ~np.array(is_mappable, dtype=bool)
        excluded_bins = np.where(is_excluded)[0]
        if bins_to_exclude is not None:
            bins_to_exclude = np.array(bins_to_exclude)
            excluded_bins = np.unique(np.concatenate((excluded_bins, bins_to_exclude),0))
            is_excluded[excluded_bins] = True

        filtered_cnvs = unfiltered_cnvs[:, ~is_excluded]

        return dict(filtered_cnvs=filtered_cnvs, excluded_bins=excluded_bins)
    else:
        return dict(unfiltered_cnvs=filtered_cnvs)
