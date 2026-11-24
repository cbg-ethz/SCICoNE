//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#include "Utils.h"

#include <algorithm>
#include <stdexcept>

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

void Utils::random_initialize_labels_map(std::map<u_int, int> &distinct_regions, int n_regions, int max_regions_per_node, double lambda_r)
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
    boost::random::poisson_distribution<int> poisson_r(lambda_r);

    // n_copies from Poisson(lambda_c)+1
    // boost::random::poisson_distribution<int> poisson_c(lambda_c); // the param is to be specified later
    // sign
    boost::random::bernoulli_distribution<double> bernoulli_05(0.5);

    int r = std::min(max_regions_per_node, poisson_r(generator) + 1); //n_regions to sample
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
     * Parses a comma-separated matrix into mat, which the caller has already sized to
     * the shape it expects: one line of the file per row of mat.
     *
     * Every value is bounds-checked before it is written. A file with more rows or
     * columns than mat used to be written past the end of it, which is undefined
     * behaviour: depending on heap layout it either tripped an assertion or crashed
     * outright. A mismatch now throws std::runtime_error naming the expected and the
     * found shape. Callers know what the rows and columns mean, so they add that.
     *
     * A row with fewer columns than expected is tolerated, as before -- probabilities do
     * not always sum to one, so a column of empty bins can be missing -- and its
     * remaining entries keep their default value.
     * */

    ifstream filein(path);
    if (!filein.is_open())
        throw std::runtime_error("could not open " + path);

    const size_t n_rows = mat.size();
    const size_t n_cols = n_rows > 0 ? mat[0].size() : 0;

    size_t i = 0;         // rows found
    size_t max_cols = 0;  // widest row found
    for (std::string line; std::getline(filein, line); )
    {
        if (line.find_first_not_of(" \t\r") == std::string::npos)
            continue; // a blank line, such as a trailing one, is not a row

        std::istringstream fline(line);
        std::string val;
        size_t j = 0;
        while (std::getline(fline, val, ','))
        {
            if (i < n_rows && j < n_cols)
            {
                try {
                    mat[i][j] = std::stod(val);
                } catch (const std::logic_error&) {
                    throw std::runtime_error(
                        path + ": row " + std::to_string(i + 1) + ", column " + std::to_string(j + 1) +
                        " is '" + val + "', which is not a number. The matrix must hold only numbers, "
                        "with no header row and no row names.");
                }
            }
            j++;
        }
        max_cols = std::max(max_cols, j);
        i++;
    }

    if (i != n_rows || max_cols > n_cols)
    {
        std::string msg = path + ": expected " + std::to_string(n_rows) + " rows by " +
                          std::to_string(n_cols) + " columns, but found " + std::to_string(i) +
                          " rows by " + std::to_string(max_cols) + " columns.";
        if (n_rows != n_cols && i == n_cols && max_cols == n_rows)
            msg += " It looks transposed.";
        throw std::runtime_error(msg);
    }
}

vector<vector<double>> Utils::condense_matrix(vector<vector<double>>& D, vector<int>& region_sizes) {

    /*
     * Groups the matrix by the given region_sizes.
     * E.g. if region_sizes[0]=5 then sum up the first 5 columns and that will be one column in the resulting matrix
     * Biological example: The bins will be grouped by the regions
     * */

    size_t n_rows = D.size();
    size_t n_regions = region_sizes.size();

    int sum_region_sizes = accumulate( region_sizes.begin(), region_sizes.end(), 0);

    vector<vector<double>> condensed_mat(n_rows, vector<double>(n_regions));

    for (unsigned i = 0; i < n_rows; ++i)
    {
        int region_id = 0;
        int region_count = 0;
        for (int j = 0; j < sum_region_sizes; ++j) {
            double to_add = D[i][j];
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

template<class T>
void Utils::read_vector(vector<T> &vec, const string &path) {

    /*
     * Reads a 1 dimensional vector file at path to the reference vec.
     * Throws std::logic_error
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

    if (vec.size() == 0)
        throw std::logic_error("The size of the vector should not be zero!");

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

template void Utils::read_vector(vector<int> &vec, const string &path);
template void Utils::read_vector(vector<double> &vec, const string &path);
