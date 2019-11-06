//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#include <iostream>
#include <string>
#include "Simulation.h"
#include <cxxopts.hpp>

#include "globals.cpp"

// globals
int print_precision;
int copy_number_limit;
double lambda_r;
double lambda_c;
double lambda_s;
double cf;
double c_penalise;
unsigned is_overdispersed;
string f_name_posfix;
int verbosity;
double eta;
// endof globals

using namespace std;

int main(int argc, char* argv[]) {

    // the default values
    int n_regions = 50;
    int n_nodes = 50;
    lambda_r = 0.1;
    lambda_c = 0.4;
    cf = 0.0;
    c_penalise = 1.0;
    int n_cells = 500;
    int n_bins = 10000;
    int n_reads = 10000;
    int n_iters = 5000;
    int max_region_size = 25;
    int ploidy = 2;
    verbosity = 0;
    int seed = -1;
    copy_number_limit = 15;
    is_overdispersed = 0;
    eta = 1e-4;
    double nu = 1.0;

    // minimum region size should be bigger than window_size
    unsigned min_region_size = 10;
    string f_name_postfix = "";

    print_precision = 15;

    cxxopts::Options options("Mcmc simulations", "Simulates the count matrix. Outputs the count matrix, region sizes, ground truth and the tree that generated the data.");
    options.add_options()
            ("n_bins", "Number of bins of the input matrix", cxxopts::value(n_bins))
            ("n_cells", "Number of cells", cxxopts::value(n_cells))
            ("n_nodes", "Number of nodes of the tree", cxxopts::value(n_nodes))
            ("n_regions", "Number of regions", cxxopts::value(n_regions))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_reads", "Number of reads per cell", cxxopts::value(n_reads))
            ("ploidy", "The ploidy information", cxxopts::value(ploidy))
            ("verbosity", "verbosity of the program", cxxopts::value(verbosity))
            ("seed", "Seed", cxxopts::value(seed))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("print_precision", "The precision points of the score values to be printed", cxxopts::value(print_precision))
            ("copy_number_limit", "the maximum copy number profile one bin or region can have", cxxopts::value(copy_number_limit))
            ("min_reg_size", "the minimum size that a region can have", cxxopts::value(min_region_size))
            ("c_penalise","term that penalises trees containing cancelling events to be added to tree event prior",cxxopts::value(c_penalise))
            ("nu","nu parameter, the overdispersion variable",cxxopts::value(nu))
            ;

    auto result = options.parse(argc, argv);


    if (result.count("seed"))
    {
        seed = result["seed"].as<int>();
        //set a seed number for reproducibility
        SingletonRandomGenerator::get_instance(seed);
    }



    Simulation sim(n_regions, n_bins, n_nodes, n_cells, n_reads, max_region_size, ploidy);

    if(result.count("nu"))
    {
        std::cout<<"Simulating with overdispersion, coefficient: " << nu << std::endl;
        is_overdispersed = 1;
        sim.tree.nu = nu;
    }
    else
    {
        std::cout<<"Simulating without overdispersion" << std::endl;
        is_overdispersed = 0;
    }
    sim.sample_cluster_sizes(n_cells);

    sim.sample_region_sizes(n_bins, min_region_size);

    sim.simulate_count_matrix(false, nu);

    sim.split_regions_to_bins();

    sim.write_output(f_name_postfix);

    std::cout<<"Successfully simulated." << std::endl;

    return EXIT_SUCCESS;

}
