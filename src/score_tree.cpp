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
unsigned is_overdispersed;
string f_name_posfix;
double eta;
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
    int n_regions;
    int ploidy=2;
    int verbosity=0;
    string file;
    string region_sizes_file = "";
    string d_matrix_file = "";
    string f_name_postfix = ""; //posfix
    double nu = 1.0;
    eta = 1e-4;

    cxxopts::Options options("Score tree", "Scores the tree written in a file");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ',', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_regions", "Number of regions in the input matrix", cxxopts::value(n_regions))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("file", "file", cxxopts::value(file))
            ("is_overdispersed", "multinomial or dirichlet multinomial in the likelihood", cxxopts::value(is_overdispersed))
            ("nu","nu parameter, the overdispersion variable",cxxopts::value(nu))
            ;
    auto result = options.parse(argc, argv);

    std::cout << "Reading the input matrix..." << std::endl;
    vector<vector<double>> d_regions(n_cells, vector<double>(n_regions));
    Utils::read_counts(d_regions, d_matrix_file);

    std::cout << "Reading the region_sizes file..." << std::endl;
    vector<int> region_sizes;
    Utils::read_vector(region_sizes, region_sizes_file);

    n_regions = region_sizes.size();

    Inference mcmc(n_regions, ploidy, verbosity);
    mcmc.initialize_from_file(file);
    mcmc.t.nu = nu;
    //mcmc.compute_neutral_table(d_regions, region_sizes);
    mcmc.compute_t_table(d_regions,region_sizes);
    mcmc.compute_t_od_scores(d_regions, region_sizes);
    mcmc.update_t_prime(); // set t_prime to t

    // write the tree
    // std::ofstream tree_file("./" +f_name_posfix+"_tree_rescored" + ".txt");
    // tree_file << mcmc.t;
    std::cout << mcmc.t;


    return EXIT_SUCCESS;
}