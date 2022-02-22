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
double cf;
double c_penalise;
int copy_number_limit;
unsigned is_overdispersed;
unsigned smoothed;
string f_name_posfix;
int verbosity;
double eta;
// endof globals

using namespace std;




int main( int argc, char* argv[]) {

    // set the globals
    print_precision = 15;
    lambda_s = 0.5;
    lambda_r = 2.0;
    lambda_c = 1.0;
    cf = 1.0;
    c_penalise = 1.0;
    copy_number_limit = 5;
    is_overdispersed = 1;
    int n_cells;
    int n_regions;
    int ploidy=2;
    verbosity=0;
    string file;
    string region_sizes_file = "";
    string d_matrix_file = "";
    string cluster_sizes_file = "";
    string f_name_postfix = ""; //posfix
    double nu = 1.0;
    eta = 1e-4;

    // region neutral states
    string region_neutral_states_file;

    cxxopts::Options options("Score tree", "Scores the tree written in a file");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ',', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("cluster_sizes_file", "Path to the file containing the cluster sizes, each line contains one cluster size", cxxopts::value(cluster_sizes_file))
            ("n_regions", "Number of regions in the input matrix", cxxopts::value(n_regions))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("tree_file", "the tree file to load", cxxopts::value(file))
            ("is_overdispersed", "multinomial or dirichlet multinomial in the likelihood", cxxopts::value(is_overdispersed))
            ("nu","nu parameter, the overdispersion variable",cxxopts::value(nu))
            ("cf", "cluster fraction variable between 0 and 1 to affect the tree prior coefficient", cxxopts::value(cf))
            ("region_neutral_states_file", "Path to the file containing the neutral state of each region to use as the root of the tree", cxxopts::value(region_neutral_states_file))
            ;
    auto result = options.parse(argc, argv);

    std::cout << "Reading the input matrix..." << std::endl;
    vector<vector<double>> d_regions(n_cells, vector<double>(n_regions));
    Utils::read_counts(d_regions, d_matrix_file);

    std::cout << "Reading the region_sizes file..." << std::endl;
    vector<int> region_sizes;
    Utils::read_vector(region_sizes, region_sizes_file);

    n_regions = region_sizes.size();

    vector<int> cluster_sizes;
    if (result.count("cluster_sizes_file"))
    {
      std::cout << "Reading the cluster_sizes file..." << std::endl;
      Utils::read_vector(cluster_sizes, cluster_sizes_file);
    }
    else
    {
      if (n_cells < 20)
        std::cout << "Warning: there are only " << n_cells <<  " observations. If these are clusters, the cluster_sizes_file parameter should be specified for accurate tree scoring.";

      cluster_sizes = std::vector<int>(n_cells, 1);
    }

    vector<int> region_neutral_states;
    if (result.count("region_neutral_states_file")) {
      std::cout << "Reading the region_neutral_states file..." << std::endl;
      Utils::read_vector(region_neutral_states, region_neutral_states_file);
    }
    else {
      std::cout << "Assuming root to have copy number state " << ploidy << " in all regions" << std::endl;
      region_neutral_states = std::vector<int>(n_regions, ploidy);
    }

    Inference mcmc(n_regions, n_cells, region_neutral_states, ploidy, verbosity);
    mcmc.initialize_from_file(file);
    mcmc.t.nu = nu;
    //mcmc.compute_neutral_table(d_regions, region_sizes);
    mcmc.compute_t_table(d_regions,region_sizes,cluster_sizes);
    mcmc.compute_t_root_scores(d_regions, region_sizes,cluster_sizes);
    mcmc.update_t_prime(); // set t_prime to t

    // write the tree
    // std::ofstream tree_file("./" +f_name_posfix+"_tree_rescored" + ".txt");
    // tree_file << mcmc.t;
    std::cout << mcmc.t;


    return EXIT_SUCCESS;
}
