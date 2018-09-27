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
    int n_regions;
    int n_nodes;
    double lambda_r;
    double lambda_c;
    int n_cells;
    int n_reads;
    int max_region_size;
    vector<vector<double>> D;
    vector<int> region_sizes;
    vector<vector<int>> ground_truth;



public:
    // constructor
    Simulation(int n_regions, int n_nodes, double lambda_r, double lambda_c, int n_cells, int n_reads,
               int max_region_size, int ploidy)
            : n_regions(n_regions),
              n_nodes(n_nodes),
              lambda_r(lambda_r),
              lambda_c(lambda_c),
              n_cells(n_cells),
              n_reads(n_reads), max_region_size(max_region_size), tree(ploidy, n_regions), ground_truth(n_cells, vector<int>(n_regions,ploidy)), D(n_cells, vector<double>(n_regions)),
              region_sizes(n_regions)
    {
        Inference mcmc(n_regions);
        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree

        vector<vector<double>> p_read_region_cell(n_cells, vector<double>(n_regions)); // initialize with the default value

        // random uniform sample the region sizes
        for (int j = 0; j < n_regions; ++j) {
            int region_size = MathOp::random_uniform(1,max_region_size);
            region_sizes[j] = region_size;
        }

        // assign cells uniformly to the nodes
        for (int i = 0; i < n_cells; ++i) {
            int uniform_val = MathOp::random_uniform(0, mcmc.t.n_regions);

            Node* uniform_node = mcmc.t.all_nodes_vec[uniform_val];

            for (auto const& x : uniform_node->c) // iterate over map, fill the existing region values, the others will be zero by default
            {
                ground_truth[i][x.first] = x.second + ploidy;
            }
        }

        // create the p_read_region_cell values, not normalized yet.
        for (int i = 0; i < p_read_region_cell.size(); ++i) {
            for (int j = 0; j < p_read_region_cell[i].size(); ++j) {
                p_read_region_cell[i][j] = ground_truth[i][j] * region_sizes[j];
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


        std::mt19937 &generator = SingletonRandomGenerator::get_generator();


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

    void infer_cnvs(int n_iters)
    {
        // move probabilities
        vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f};
        Inference mcmc(n_regions);
        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree

        mcmc.compute_t_table(D,region_sizes);
        mcmc.infer_mcmc(D,region_sizes, move_probs, n_iters);


    }


};


#endif //SC_DNA_SIMULATION_H
