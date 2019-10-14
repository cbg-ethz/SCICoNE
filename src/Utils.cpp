//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#include "Utils.h"

uint64_t Utils::calculate_hash(const void *buffer, size_t length) {

    /*
     * Calculates the hash using xxhash.
     * */

    unsigned long long const seed = SingletonRandomGenerator::get_instance().seed;
    unsigned long long const hash = XXH64(buffer, length, seed);
    return hash;

}

bool Utils::is_empty_map(std::map<u_int, int> &dict){
    /*
     * Returns true if the map is empty, e.g. all zero
     * */
    for (auto const &it : dict)
    {
        if(it.second != 0)
            return false;
    }
    return true;
}

void Utils::random_initialize_labels_map(std::map<u_int, int> &distinct_regions, int n_regions)
{
    /*
     * Initializes an empty map to represent the c_change attribute of a node.
     * Modifies the distinct_regions map
     * Throws out_of_range exception.
     * */

    assert(is_empty_map(distinct_regions));

    // sample the number of regions to be affected with Poisson(lambda_r)+1
    std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;

    // n_regions from Poisson(lambda_R)+1
    // boost::random::poisson_distribution<int> poisson_r(lambda_r); // the param is to be specified later

    // n_copies from Poisson(lambda_c)+1
    // boost::random::poisson_distribution<int> poisson_c(lambda_c); // the param is to be specified later
    // sign
    boost::random::bernoulli_distribution<double> bernoulli_05(0.5);

    int r = 1; // poisson_r(generator) + 1; //n_regions to sample
    // sample r distinct regions uniformly
    int regions_sampled = 0;
    // if r>n_regions then reject the move. n_regions: max region index
    if (r > n_regions)
        throw std::out_of_range("There cannot be more than n_regions amount of distinct regions");

    while (regions_sampled < r)
    {
        int uniform_val = MathOp::random_uniform(0, n_regions-1);

        if (distinct_regions.find(uniform_val) == distinct_regions.end()) // not found
        {
            int n_copies = 1; // poisson_c(generator) + 1;
            bool sign = bernoulli_05(generator);
            distinct_regions[uniform_val] = (sign? n_copies : -n_copies); // key:region id, val: copy number change
            regions_sampled++;
        }
    }

}

void Utils::read_counts(vector<vector<double>> &mat, const string &path) {

    /*
     * Parses the input data into a default filled double vector of vector.
     * */

    ifstream filein(path);

    unsigned i = 0;
    for (std::string line; std::getline(filein, line); )
    {

        std::istringstream fline(line);
        std::string val;
        int j = 0;
        while (std::getline(fline, val, ','))
        {
            mat[i][j] = stod(val);
            j++;
        }
        // assert(j == mat[i].size()); probabilities don't always sum up to one, we can get 1 column of empty bins here
        i++;
    }

    assert(i == mat.size());



}
template <class T>
vector<vector<T>> Utils::condense_matrix(vector<vector<T>>& D, vector<int>& region_sizes) {

    /*
     * Groups the matrix by the given region_sizes.
     * E.g. if region_sizes[0]=5 then sum up the first 5 columns and that will be one column in the resulting matrix
     * Biological example: The bins will be grouped by the regions
     * */

    size_t n_rows = D.size();
    size_t n_regions = region_sizes.size();

    int sum_region_sizes = accumulate( region_sizes.begin(), region_sizes.end(), 0);

    vector<vector<T>> condensed_mat(n_rows, vector<T>(n_regions));

    for (unsigned i = 0; i < n_rows; ++i)
    {
        int region_id = 0;
        int region_count = 0;
        for (int j = 0; j < sum_region_sizes; ++j) {
            T to_add = D[i][j];
            condensed_mat[i][region_id] += to_add;
            region_count++;
            if(region_count == region_sizes[region_id])
            {
                region_id++;
                region_count = 0;
            }
        }

    }


    return condensed_mat;
}

void Utils::read_vector(vector<int> &vec, const string &path) {

    /*
     * Reads a 1 dimensional vector file at path path to reference vec.
     * */

    ifstream filein(path);

    for (std::string line; std::getline(filein, line); )
    {
        istringstream fline(line);
        std::string val;
        while (std::getline(fline, val))
        {
            vec.push_back(stod(val)); // push_back is fine since this file is much smaller
        }
    }

    assert(vec.size() != 0);

}

vector<vector<int>> Utils::regions_to_bins_cnvs(vector<vector<int>> &cnvs, vector<int> &region_sizes) {
    /*
     * Creates a cells-bins cnvs matrix from the regions-bins cnvs matrix by duplicating the region columns by their sizes.
     * */

    // compute the total n_bins by summing up the region sizes
    double n_bins = accumulate( region_sizes.begin(), region_sizes.end(), 0.0);
    double n_rows = cnvs.size();
    vector<vector<int>> bins_vec(n_rows, vector<int>(n_bins));

    for (unsigned i = 0; i < n_rows; ++i) {
        int region_offset = 0;

        for (unsigned j = 0; j < cnvs[0].size(); ++j) {

            for (int k = 0; k < region_sizes[j]; ++k) {
                bins_vec[i][k + region_offset] = cnvs[i][j];
            }
            region_offset += region_sizes[j];
        }
    }

    return bins_vec;
}

map<u_int, int> Utils::map_diff(map<u_int, int> a, map<u_int, int> b) {
    /*
     * Computes the difference between all of the elements existing in either of the maps.
     * */

    map<u_int, int> difference;

    // store all the keys in a set
    set<u_int> all_keys;

    for (auto const &e : a)
        all_keys.insert(e.first);
    for (auto const &e : b)
        all_keys.insert(e.first);

    for (auto const &k : all_keys)
    {
        int diff = a[k] - b[k];
        if (diff != 0)
            difference[k] = a[k] - b[k];
    }

    return difference;
}

template vector<vector<int>> Utils::condense_matrix(vector<vector<int>>& D, vector<int>& region_sizes);
template vector<vector<double>> Utils::condense_matrix(vector<vector<double>>& D, vector<int>& region_sizes);