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
    vector<vector<double>> d_bins;
    vector<int> region_sizes;
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
              D(n_cells, vector<double>(n_regions)), d_bins(n_cells, vector<double>(n_bins)), region_sizes(n_regions)
    {

    }


    void simulate_count_matrix(bool is_neutral, double nu)
    {
        /*
         * Simulates the count matrix and the ground_truth.
         * If is_neutral, then there are no mutations.
         * nu is the overdispersion parameter.
         * All values of the ground truth will be equal to ploidy and the count matrix reads distributions will only be based on region sizes.
         * Else, a mutation tree will be inferred and that tree will simulate the data.
         * */



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

        //  ground truth: convert regions->bins
        vector<vector<int>> ground_truth_bins = Utils::regions_to_bins_cnvs(this->ground_truth, this->region_sizes);
        this->ground_truth = ground_truth_bins;

        // create the p_read_bin_cell values, not normalized yet.
        for (std::size_t i = 0; i < p_read_bin_cell.size(); ++i) {
            for (std::size_t j = 0; j < p_read_bin_cell[i].size(); ++j) {
                if (not is_neutral)
                    p_read_bin_cell[i][j] = ground_truth[i][j];
                else
                    p_read_bin_cell[i][j] = ploidy;
            }
        }

        if (not is_overdispersed)
        {
            // normalize the p_read_region cell per cell to get probabilities. e.g. prob of a read belonging to a region in a cell
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
            for (std::size_t i = 0; i < alphas.size(); ++i) {
                for (std::size_t j = 0; j < alphas[i].size(); ++j) {
                    alphas[i][j] = nu * p_read_bin_cell[i][j];
                }
            }

            for (std::size_t i = 0; i < p_read_bin_cell.size(); ++i) {
                p_read_bin_cell[i] = MathOp::dirichlet_sample(alphas[i]);
            }
        }


        for (std::size_t i = 0; i < d_bins.size(); ++i) // for each cell
        {
            // assign the read to region by sampling from the dist
            boost::random::discrete_distribution<> d(p_read_bin_cell[i].begin(), p_read_bin_cell[i].end()); // distribution will be different for each cell

            for (int j = 0; j < n_reads; ++j) // distribute the reads to regions
            {
                unsigned sample = d(generator);
                d_bins[i][sample]++;
            }
        }

        D = Utils::condense_matrix(d_bins, region_sizes);

        if (not is_neutral) // do not compute the tree for the null model
        {
            // compute the tree and store it in this->tree
            mcmc.compute_t_table(D,region_sizes);
            mcmc.compute_t_od_scores(D, region_sizes);
        }

        this->tree = mcmc.t;

    }

    void sample_region_sizes(int n_bins, unsigned min_width = 1)
    {
        /*
         * Uniformly samples the region sizes and returns the vector of region sizes.
         * */

        vector<long double> dirichlet = MathOp::dirichlet_sample(n_regions);

        n_bins -= min_width * n_regions;

        long double sum_prob = 0.0;

        for (int i = 0; i < n_regions; ++i) {
            region_sizes[i] = static_cast<int>((sum_prob+dirichlet[i])*n_bins) - static_cast<int>(sum_prob*n_bins);
            sum_prob += dirichlet[i];
            region_sizes[i] += min_width;
        }

        long double sum = accumulate( region_sizes.begin(), region_sizes.end(), 0.0); //for debug purposes
        n_bins += min_width * n_regions;

        if (sum != n_bins)
            cout << "WARNING: sum of all of the region sizes is not equal to n_bins since probabilities don't always sum up to 1 due to precision loss"
                    " (e.g. they can sum up to 0.999999999999997) as a result the last bin may have zero count." <<endl;

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
        for (auto const &v1: this->d_bins) {
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
