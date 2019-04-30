//
// Created by Tuncel  Mustafa Anil on 11/5/18.
//
#include "globals.cpp"
#include "Inference.h"
#include <vector>
#include <iostream>
#include <string>
#include <cxxopts.hpp>
#include "Utils.h"


// globals
int print_precision;
double lambda_s;
double lambda_r;
double lambda_c;
double c_penalise;
int copy_number_limit;
double v;
unsigned tree_prior_in_chi;
unsigned is_overdispersed;

// endof globals

using namespace std;




int main( int argc, char* argv[]) {

    // set the globals
    print_precision = 15;
    lambda_s = 0.5;
    lambda_r = 2.0;
    lambda_c = 1.0;
    c_penalise = 1.0;
    copy_number_limit = 5;
    is_overdispersed = 1;
    int n_cells;
    int n_bins = 10000;
    size_t n_regions;
    int ploidy=2;
    int verbosity=0;
    string file;
    string region_sizes_file = "";
    string d_matrix_file = "";
    string f_name_postfix = ""; //posfix

    v = std::nan("");
    tree_prior_in_chi = 1;

    cxxopts::Options options("Score tree", "Scores the tree written in a file");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            ("v","v value used in size changing moves",cxxopts::value(v))
            ("tree_prior_chi", "whether to include the tree prior in X", cxxopts::value(tree_prior_in_chi))
            ("file", "file", cxxopts::value(file))
            ;
    auto result = options.parse(argc, argv);


    // read the input d_matrix
    vector<vector<double>> d_bins(n_cells, vector<double>(n_bins));
    Utils::read_counts(d_bins, d_matrix_file);

    // read the region_sizes file
    vector<int> region_sizes;
    vector<vector<double>> d_regions;

    Utils::read_vector(region_sizes, region_sizes_file);
    // Merge the bins into regions
    d_regions = Utils::condense_matrix(d_bins, region_sizes);
    n_regions = region_sizes.size();

    Inference mcmc(n_regions, ploidy, verbosity);
    mcmc.initialize_from_file(file);

    //mcmc.compute_neutral_table(d_regions, region_sizes);
    mcmc.compute_t_table(d_regions,region_sizes);

    // write the tree
    // std::ofstream tree_file("./" +f_name_postfix+"_tree_rescored" + ".txt");
    // tree_file << mcmc.t;
    std::cout << mcmc.t;


    return EXIT_SUCCESS;
}