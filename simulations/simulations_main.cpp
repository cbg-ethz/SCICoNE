//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#include "Inference.h"
#include <iostream>
#include "SingletonRandomGenerator.h"

using namespace std;



void simulate(int n_regions, int n_nodes, double lambda_r, double lambda_c, int n_cells, int n_reads, int ploidy=2)
{

    int max_region_size = 10;

    // generate the random tree

    Inference mcmc(n_regions);
    mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree

    vector<vector<double>> p_read_region_cell(n_cells, vector<double>(n_regions)); // initialize with the default value

    vector<int> region_sizes(n_regions); // sampling the region sizes
    for (int j = 0; j < n_regions; ++j) {
        int region_size = MathOp::random_uniform(1,max_region_size);
        region_sizes[j] = region_size;
    }

    // assign cells uniformly to the nodes
    vector<vector<int>> cell_regions(n_cells, vector<int>(n_regions, ploidy)); //fill constructor, ploidy is the default value

    for (int i = 0; i < n_cells; ++i) {
        int uniform_val = MathOp::random_uniform(0, mcmc.t.n_regions);

        Node* uniform_node = mcmc.t.all_nodes_vec[uniform_val];

        for (auto const& x : uniform_node->c) // iterate over map, fill the existing region values, the others will be zero by default
        {
            cell_regions[i][x.first] = x.second + ploidy;
        }
    }

    // create the p_read_region_cell values, not normalized yet.
    for (int i = 0; i < p_read_region_cell.size(); ++i) {
        for (int j = 0; j < p_read_region_cell[i].size(); ++j) {
            p_read_region_cell[i][j] = cell_regions[i][j] * region_sizes[j];
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

    // create the D matrix
    vector<vector<double>> D(n_cells, vector<double>(n_regions)); // initialize with the default value
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



//
//    std::discrete_distribution<> d(weights.begin(), weights.end());
//
//    unsigned sample = d(generator);




    cout << "debug";



}

int main(int argc, char* argv[])
{
    // define n_regions
    int n_regions = 50;
    int n_nodes = 50;
    double lambda_r = 0.1;
    double lambda_c = 0.2;
    int n_cells = 500;
    int n_reads = 10000;

    simulate(n_regions, n_nodes, lambda_r, lambda_c, n_cells, n_reads);



}

