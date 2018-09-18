//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#include "Utils.h"

uint64_t Utils::calcul_hash(const void *buffer, size_t length) {

    /*
     * Calculates the hash using xxhash.
     * */

    unsigned long long const seed = 0;   /* or any other value */
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

void Utils::initialize_labels_map(std::map<u_int, int> &distinct_regions, int n_regions, double lambda_r, double lambda_c)
{
    /*
     * Initializes an empty map to represent the c_change attribute of a node.
     * Modifies the distinct_regions map
     * Throws out_of_range exception.
     * */

    assert(is_empty_map(distinct_regions));

    // sample the number of regions to be affected with Poisson(lambda_r)+1
    std::mt19937 &generator = SingletonRandomGenerator::get_generator();

    // n_regions from Poisson(lambda_R)+1
    std::poisson_distribution<int> poisson_r(lambda_r); // the param is to be specified later

    // n_copies from Poisson(lambda_c)+1
    std::poisson_distribution<int> poisson_c(lambda_c); // the param is to be specified later
    // sign
    std::bernoulli_distribution bernoulli_05(0.5);

    int r = poisson_r(generator) + 1; //n_regions to sample
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
            int n_copies = poisson_c(generator) + 1;
            bool sign = bernoulli_05(generator);
            distinct_regions[uniform_val] = (sign? n_copies : -n_copies); // key:region id, val: copy number change
            regions_sampled++;
        }
    }

}

