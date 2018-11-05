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

// endof globals

using namespace std;




int main( int argc, char* argv[]) {


    print_precision = 17;
    lambda_s = 0.5;
    lambda_r = 2.0;
    lambda_c = 1.0;


    int n_cells;
    int n_bins = 10000;
    int n_regions;
    int ploidy=2;
    int verbosity=2;
    string file;
    string region_sizes_file = "";
    string d_matrix_file = "";
    string f_name_postfix = ""; //posfix

    cxxopts::Options options("Score tree", "Scores the tree written in a file");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("file", "file", cxxopts::value(file));
    auto result = options.parse(argc, argv);

    if (result.count("region_sizes_file"))
    {
        region_sizes_file = result["region_sizes_file"].as<string>();
    }

    if (result.count("d_matrix_file"))
    {
        d_matrix_file = result["d_matrix_file"].as<string>();
    }
    if (result.count("postfix")) {
        f_name_postfix = result["postfix"].as<string>();
    }
    if (result.count("file"))
    {
        file = result["file"].as<string>();
        //set a seed number for reproducibility
    }

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

    mcmc.compute_neutral_table(d_regions, region_sizes);
    mcmc.compute_t_table(d_regions,region_sizes);

    // write the tree
    std::ofstream tree_file("./" +f_name_postfix+"_tree_rescored" + ".txt");
    tree_file << mcmc.t;



    return EXIT_SUCCESS;
}