//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#ifndef SC_DNA_SIMULATION_H
#define SC_DNA_SIMULATION_H

#include <vector>
#include "Inference.h"
#include "SingletonRandomGenerator.h"

using namespace std;


class Simulation {


public:

    Tree tree;
    int ploidy;
    int n_regions;
    int n_nodes;
    double lambda_r;
    double lambda_c;
    int n_cells;
    int n_reads;
    int max_region_size;
    double delta;

    vector<vector<double>> D;
    vector<int> region_sizes;
    vector<vector<int>> ground_truth; // default init value: ploidy
    vector<vector<int>> inferred_cnvs;

public:
    // constructor
    Simulation(int n_regions, int n_nodes, double lambda_r, double lambda_c, int n_cells, int n_reads,
               int max_region_size, int ploidy, int verbosity)
            : n_regions(n_regions),
              n_nodes(n_nodes),
              lambda_r(lambda_r),
              lambda_c(lambda_c),
              ploidy(ploidy),
              n_cells(n_cells),
              n_reads(n_reads), max_region_size(max_region_size), tree(ploidy, n_regions), ground_truth(n_cells, vector<int>(n_regions,ploidy)), inferred_cnvs(n_cells, vector<int>(n_regions,ploidy)), D(n_cells, vector<double>(n_regions)), region_sizes(n_regions)
    {

    }

    void simulate_count_matrix(bool is_neutral, int verbosity)
    {
        /*
         * Simulates the count matrix and the ground_truth.
         * If is_neutral, then there are no mutations.
         * All values of the ground truth will be equal to ploidy and the count matrix reads distributions will only be based on region sizes.
         * Else, a mutation tree will be inferred and that tree will simulate the data.
         * */



        Inference mcmc(n_regions, ploidy, verbosity);
        vector<vector<double>> p_read_region_cell(n_cells, vector<double>(n_regions)); // initialize with the default value

        std::mt19937 &generator = SingletonRandomGenerator::get_generator();

        if (not is_neutral) // tree will generate the data
        {
            mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree, mcmc.t

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

            // create the p_read_region_cell values, not normalized yet.
            for (int i = 0; i < p_read_region_cell.size(); ++i) {
                for (int j = 0; j < p_read_region_cell[i].size(); ++j) {

                    if (not is_neutral)
                        p_read_region_cell[i][j] = ground_truth[i][j] * region_sizes[j];
                    else
                        p_read_region_cell[i][j] = ploidy * region_sizes[j];
                }
            }


            // normalize the p_read_region cell per cell to get probabilities. e.g. prob of a read belonging to a region in a cell
            for (int k = 0; k < p_read_region_cell.size(); ++k) {
                // find the sum value
                double sum_per_cell = accumulate(p_read_region_cell[k].begin(), p_read_region_cell[k].end(), 0.0);

                for (int i = 0; i < p_read_region_cell[k].size(); ++i) {
                    if (p_read_region_cell[k][i] != 0.0)
                        p_read_region_cell[k][i] /= sum_per_cell; // divide all the values by the sum
                }
                assert(abs(1.0 - accumulate(p_read_region_cell[k].begin(), p_read_region_cell[k].end(), 0.0)) <= 0.01); // make sure probs sum up to 1
            }
            for (int i = 0; i < D.size(); ++i) // for each cell
            {
                // assign the read to region by sampling from the dist
                std::discrete_distribution<> d(p_read_region_cell[i].begin(), p_read_region_cell[i].end()); // distribution will be different for each cell

                for (int j = 0; j < n_reads; ++j) // distribute the reads to regions
                {
                    unsigned sample = d(generator);
                    D[i][sample]++;
                }
            }
        }


    void split_regions_to_bins()
    {

        // compute the total n_bins by summing up the region sizes
        double n_bins = accumulate( region_sizes.begin(), region_sizes.end(), 0.0);



        vector<vector<double>> D_bins(n_cells, vector<double>(n_bins));
        vector<vector<int>> ground_truth_bins(n_cells, vector<int>(n_bins));

        for (int i = 0; i < D.size(); ++i) // for each cell
        {
            int region_offset = 0;
            for (int j = 0; j < D[0].size(); ++j) // for each region
            {

                for (int k = 0; k < D[i][j]; ++k) // for the count value
                {
                    int val =  MathOp::random_uniform(0, region_sizes[j]-1);
                    D_bins[i][val + region_offset]++;
                }

                for (int l = 0; l < region_sizes[j]; ++l) {
                    ground_truth_bins[i][l + region_offset] = ground_truth[i][j];
                }

                region_offset += region_sizes[j];
            }
        }

        // Unit test
        double validation_sum = 0.0;
        for (int m = 0; m < region_sizes[1]; ++m) {
            validation_sum += D_bins[0][m + region_sizes[0]]; // region_sizes[0] is the offset
        }
        assert(validation_sum == D[0][1]);

        D = D_bins;
        ground_truth = ground_truth_bins;
    }

    void sample_region_sizes()
    {
        /*
         * Uniformly samples the region sizes and returns the vector of region sizes.
         * */


        // random uniform sample the region sizes
        for (int j = 0; j < n_regions; ++j) {
            int region_size = MathOp::random_uniform(1,max_region_size);
            region_sizes[j] = region_size;
        }
    }



    void infer_cnvs(int n_iters, int verbosity)
    {

        sample_region_sizes();
        simulate_count_matrix(false, verbosity);

        Inference mcmc(n_regions, ploidy, verbosity);
        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // re-creates a random tree as mcmc.t

        // compute the initial tree using D and region sizes
        mcmc.compute_t_table(D,region_sizes);
        // move probabilities
        vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f};
        mcmc.infer_mcmc(D,region_sizes, move_probs, n_iters);

        inferred_cnvs = mcmc.assign_cells_to_nodes(D, region_sizes);


        // compute the Frobenius avg. of the difference of the inferred CNVs and the ground truth
        delta = MathOp::frobenius_avg(inferred_cnvs, ground_truth);

        if (verbosity > 1)
            mcmc.write_best_tree();
    }



    void write_d_vector(string f_name_postfix)
    {
        /*
         * Writes the d vector to file.
         * */

        std::ofstream delta_file(to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_" + f_name_postfix +  "_deltas.csv");
        delta_file << delta;
    }
};


#endif //SC_DNA_SIMULATION_H
