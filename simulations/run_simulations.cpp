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
double lambda_r;
double lambda_c;
double lambda_s;
// endof globals

using namespace std;

int main(int argc, char* argv[]) {

    // the default values
    int n_regions = 50;
    int n_nodes = 50;
    lambda_r = 0.1;
    lambda_c = 0.2;
    int n_cells = 500;
    int n_bins = 10000;
    int n_reads = 10000;
    int n_iters = 5000;
    int max_region_size = 10;
    int ploidy = 2;
    int verbosity = 0;
    int seed = 0;
    string f_name_postfix = "";

    print_precision = 16;

    cxxopts::Options options("Mcmc simulations", "Simulates the count matrix. Outputs the count matrix, region sizes, ground truth and the tree that generated the data.");
    options.add_options()
            ("n_bins", "Number of bins of the input matrix", cxxopts::value(n_bins))
            ("n_cells", "Number of cells", cxxopts::value(n_cells))
            ("n_nodes", "Number of nodes of the tree", cxxopts::value(n_nodes))
            ("n_regions", "Number of regions", cxxopts::value(n_regions))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_reads", "Number of reads per cell", cxxopts::value(n_reads))
            ("ploidy", "The ploidy information", cxxopts::value(ploidy))
            ("verbosity", "Verbosity of the programme, 1 provides standard output, 2 also writes several files", cxxopts::value(verbosity))
            ("seed", "Seed", cxxopts::value(seed))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("print_precision", "The precision points of the score values to be printed", cxxopts::value(print_precision))
            ("lambda_r","lambda param for the poisson that generates the number of regions", cxxopts::value(lambda_r))
            ("lambda_c","lambda param for the poisson that generates the copy number state of a region", cxxopts::value(lambda_c))
            ;

    auto result = options.parse(argc, argv);

    if (result.count("n_bins")) {
        n_bins = result["n_bins"].as<int>();
    }
    if (result.count("n_nodes")) {
        n_nodes = result["n_nodes"].as<int>();
    }
    if (result.count("n_cells")) {
        n_cells = result["n_cells"].as<int>();
    }
    if (result.count("n_regions")) {
        n_regions = result["n_regions"].as<int>();
    }
    if (result.count("n_reads")) {
        n_reads = result["n_reads"].as<int>();
    }
    if (result.count("n_iters")) {
        n_iters = result["n_iters"].as<int>();
    }
    if (result.count("ploidy")) {
        ploidy = result["ploidy"].as<int>();
    }
    if (result.count("verbosity")) {
        verbosity = result["verbosity"].as<int>();
    }
    if (result.count("postfix")) {
        f_name_postfix = result["postfix"].as<string>();
    }
    if (result.count("seed"))
    {
        seed = result["seed"].as<int>();
        //set a seed number for reproducibility
        SingletonRandomGenerator::get_generator(seed);
    }
    if (result.count("print_precision")) {
        print_precision = result["print_precision"].as<int>();
    }
    if (result.count("lambda_r")) {
        lambda_r = result["lambda_r"].as<double>();
    }
    if (result.count("lambda_c")) {
        lambda_c = result["lambda_c"].as<double>();
    }


    Simulation sim(n_regions, n_bins, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy,
                   verbosity);

    sim.sample_region_sizes(n_bins, 1);
    sim.simulate_count_matrix(false, verbosity);
    sim.split_regions_to_bins();

    sim.write_output(f_name_postfix);

    return EXIT_SUCCESS;

}
