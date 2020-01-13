//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#ifndef SC_DNA_SIMULATION_H
#define SC_DNA_SIMULATION_H

#include <vector>
#include "Inference.h"
#include "SingletonRandomGenerator.h"

#include <boost/random/discrete_distribution.hpp>

using namespace std;


class Simulation {


public:

    Tree tree; // true tree that generates the data
    int ploidy;
    int n_bins;
    int n_regions;
    int n_nodes;
    int n_cells;
    int n_reads;
    int max_region_size;

    vector<vector<double>> D;
    vector<int> region_sizes;
    vector<int> cluster_sizes;
    vector<vector<int>> ground_truth; // default init value: ploidy

public:
    // constructor
    Simulation(int n_regions, int n_bins, int n_nodes, int n_cells, int n_reads, int max_region_size, int ploidy)
            : n_regions(n_regions),
              n_bins(n_bins),
              n_nodes(n_nodes),
              ploidy(ploidy),
              n_cells(n_cells),
              n_reads(n_reads), max_region_size(max_region_size), tree(ploidy, n_regions),
              ground_truth(n_cells, vector<int>(n_regions, ploidy)),
              D(n_cells, vector<double>(n_bins)), region_sizes(n_regions), cluster_sizes(n_cells, 1)
    {

    }

    void sample_region_sizes(int n_bins, unsigned min_width = 1)
    {
        /*
         * Uniformly samples the region sizes and returns the vector of region sizes.
         * */

        vector<long double> rs_probs = MathOp::dirichlet_sample(n_regions);

        n_bins -= min_width * n_regions;

        long double sum_prob = 0.0;

        for (int i = 0; i < n_regions; ++i) {
            region_sizes[i] = static_cast<int>((sum_prob+rs_probs[i])*n_bins) - static_cast<int>(sum_prob*n_bins);
            sum_prob += rs_probs[i];
            region_sizes[i] += min_width;
        }

        long double sum = accumulate( region_sizes.begin(), region_sizes.end(), 0.0); //for debug purposes
        n_bins += min_width * n_regions;

        if (sum != n_bins)
            cout << "WARNING: sum of all of the region sizes is not equal to n_bins since probabilities don't always sum up to 1 due to precision loss"
                    " (e.g. they can sum up to 0.999999999999997). As a result the last bin may have zero counts." << endl;

    }

    void simulate_count_matrix(bool is_neutral, double nu)
    {
        /*
         * Simulates the count matrix and the ground_truth.
         * If is_neutral, then there are no mutations and all values of the ground truth will be equal to ploidy.
         * Otherwise, a mutation tree will be inferred and that tree will simulate the data.
         * nu is the overdispersion parameter.
         */

        Inference mcmc(n_regions, ploidy, verbosity);
        vector<vector<double>> p_read_bin_cell(n_cells, vector<double>(n_bins)); // initialize with the default value

        std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;

        if (not is_neutral) // tree will generate the data
        {
            mcmc.random_initialize(n_nodes, n_regions, 10000); // creates a random tree, mcmc.t

            // assign cells uniformly to the nodes
            for (int i = 0; i < n_cells; ++i) {
                int uniform_val = MathOp::random_uniform(0, mcmc.t.get_n_nodes());

                Node *uniform_node = mcmc.t.all_nodes_vec[uniform_val];

                for (auto const &x : uniform_node->c) // iterate over map, fill the existing region values, the others will be zero by default
                {
                    ground_truth[i][x.first] = x.second + ploidy;
                }
            }
        }
        mcmc.t.nu = mcmc.t_prime.nu = nu; // set nu

        // Split regions into bins to generate counts
        vector<vector<int>> ground_truth_bins(n_cells, vector<int>(n_bins));

        for (std::size_t i = 0; i < n_cells; ++i) // for each cell
        {
            int region_offset = 0;
            for (std::size_t j = 0; j < n_regions; ++j) // for each region
            {
                for (int l = 0; l < region_sizes[j]; ++l) {
                    ground_truth_bins[i][l + region_offset] = ground_truth[i][j];
                }

                region_offset += region_sizes[j];
            }
        }
        ground_truth = ground_truth_bins;

        // create the unnormalised p_bin_region_cell values
        for (std::size_t i = 0; i < p_read_bin_cell.size(); ++i) { // for each cell
            for (std::size_t j = 0; j < p_read_bin_cell[i].size(); ++j) { // for each bin
                if (not is_neutral)
                    p_read_bin_cell[i][j] = ground_truth[i][j];
                else
                    p_read_bin_cell[i][j] = ploidy;
            }
        }

        // create the normalised p_bin_region_cell values
        if (not is_overdispersed)
        {
            // normalize the p_read_region cell per cell to get probabilities. e.g. prob of a read belonging to a bin in a cell
            for (std::size_t k = 0; k < p_read_bin_cell.size(); ++k) {
                // find the sum value
                double sum_per_cell = accumulate(p_read_bin_cell[k].begin(), p_read_bin_cell[k].end(), 0.0);

                for (std::size_t i = 0; i < p_read_bin_cell[k].size(); ++i) {
                    if (p_read_bin_cell[k][i] != 0.0)
                        p_read_bin_cell[k][i] /= sum_per_cell; // divide all the values by the sum
                }
                assert(abs(1.0 - accumulate(p_read_bin_cell[k].begin(), p_read_bin_cell[k].end(), 0.0)) <= 0.01); // make sure probs sum up to 1
            }
        }
        else
        {
            // set alphas
            vector<vector<double>> alphas(n_cells, vector<double>(n_bins));
            for (std::size_t i = 0; i < alphas.size(); ++i) { // each cell
                for (std::size_t j = 0; j < alphas[i].size(); ++j) { // each bin
                    alphas[i][j] = nu * p_read_bin_cell[i][j];
                }
            }

            for (std::size_t i = 0; i < p_read_bin_cell.size(); ++i) {
                p_read_bin_cell[i] = MathOp::dirichlet_sample(alphas[i]);
            }
        }


        // count number of reads in each cell
        for (std::size_t i = 0; i < D.size(); ++i) // for each cell
        {
            // assign the read to region by sampling from the dist
            boost::random::discrete_distribution<> d(p_read_bin_cell[i].begin(), p_read_bin_cell[i].end()); // distribution will be different for each cell

            for (int j = 0; j < n_reads; ++j) // distribute the reads to regions
            {
                unsigned sample = d(generator);
                D[i][sample]++;
            }
        }

        // compute score of tree for generated data
        vector<vector<double>> D_regions(n_cells, vector<double>(n_regions));
        if (not is_neutral) // do not compute the tree for the null model
        {

            // Aggregate reads within regions to score tree on generated data
            for (std::size_t i = 0; i < n_cells; ++i) { // for each cell
              int region_id = 0;
              int region_count = 0;
              for (std::size_t j = 0; j < n_bins; ++j) { // for each bin
                D_regions[i][region_id] += D[i][j];
                region_count += 1;
                if (region_count == region_sizes[region_id]) {
                  region_id += 1;
                  region_count = 0;
                }
              }
            }

            // compute the tree and store it in this->tree
            mcmc.compute_t_table(D_regions,region_sizes,cluster_sizes);
            mcmc.compute_t_od_scores(D_regions, region_sizes,cluster_sizes);
        }

        this->tree = mcmc.t;

        // Unit test
        double validation_sum = 0.0;
        for (int m = 0; m < region_sizes[1]; ++m) {
            validation_sum += D[0][m + region_sizes[0]]; // region_sizes[0] is the offset
        }
        assert(validation_sum == D_regions[0][1]);

    }

    void write_output(const string& f_name_postfix)
    {
        /*
         * Writes the simulation outputs to text files.
         * The outputs are: D matrix, region sizes, ground truth and the tree that generated the data
         * */

        // write the region sizes
        std::ofstream region_sizes_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_region_sizes.txt");
        for (const auto &e : this->region_sizes) region_sizes_file << e << "\n";

        // write the ground truth
        std::ofstream ground_truth_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_ground_truth.csv");
        for (auto const &v1: this->ground_truth) {
            for (size_t i = 0; i < v1.size(); i++)
            {
                if (i == v1.size()-1) // the last element
                    ground_truth_file << v1[i];
                else // add comma
                    ground_truth_file << v1[i] << ',';
            }
            ground_truth_file << std::endl;
        }

        // write the D matrix
        std::ofstream D_mat_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_d_mat.csv");
        for (auto const &v1: this->D) {
            for (size_t i = 0; i < v1.size(); i++)
            {
                if (i == v1.size()-1) // the last element
                    D_mat_file << v1[i];
                else // add comma
                    D_mat_file << v1[i] << ',';
            }
            D_mat_file << std::endl;
        }

        // write the tree that generated the data
        std::ofstream tree_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_tree.txt");
        tree_file << this->tree;

    }

};


#endif //SC_DNA_SIMULATION_H
