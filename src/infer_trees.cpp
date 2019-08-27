#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <limits>
#include "MathOp.h"
#include "Tree.h"
#include "Inference.h"

#include <chrono> // for measuring the execution time

#include <cxxopts.hpp>
#include "globals.cpp"



// globals
int print_precision;
int copy_number_limit;
double lambda_s;
double lambda_r;
double lambda_c;
double c_penalise;
unsigned is_overdispersed;
double eta;

// endof globals

using namespace std;
using namespace std::chrono;


int main( int argc, char* argv[]) {

    int n_iters = 5000; // the default value is 5000 iterations.
    int n_cells;
    int n_bins = 10000;
    int ploidy = 2;
    int verbosity = 0;
    int seed = 0;
    string f_name_postfix;
    string region_sizes_file;
    string d_matrix_file;

    unsigned size_limit = std::numeric_limits<unsigned>::max();

    size_t n_regions;
    size_t n_regions_initial = 0; // used for naming the output for I/O workflow purposes

    // random tree parameters
    int n_nodes = 3;
    lambda_r = 0.1;
    lambda_c = 0.2;
    c_penalise = 1.0;
    is_overdispersed = 1;
    eta = 1e-4;

    int n_reads = -1; // -1 means not specified

    // set the globals
    print_precision = 15;
    lambda_s = 0.5;
    copy_number_limit = 5;

    double nu = 1.0;
    bool random_init = true;

    // move probabilities
    vector<float> move_probs;

    cxxopts::Options options("Single cell CNV inference", "finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_regions", "Number of regions to be contained in the output file", cxxopts::value(n_regions))
            ("n_reads", "Number of reads to be contained in the output file", cxxopts::value(n_reads))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("verbosity", "verbosity", cxxopts::value(verbosity))
            ("seed", "seed", cxxopts::value(seed))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            ("size_limit", "the limitation on the max size of the tree", cxxopts::value(size_limit))
            ("copy_number_limit", "the maximum copy number profile one bin or region can have", cxxopts::value(copy_number_limit))
            // random tree parameters
            ("n_nodes","the number of nodes in the random initialised tree", cxxopts::value(n_nodes))
            ("lambda_r","lambda param for the poisson that generates the number of regions", cxxopts::value(lambda_r))
            ("lambda_c","lambda param for the poisson that generates the copy number state of a region", cxxopts::value(lambda_c))
            ("c_penalise","term that penalises trees containing cancelling events to be added to tree event prior",cxxopts::value(c_penalise))
            // ("is_overdispersed", "multinomial or dirichlet multinomial in the likelihood", cxxopts::value(is_overdispersed))
            ("nu","nu parameter, the overdispersion variable",cxxopts::value(nu))
            ("random_init","Boolean parameter to enable random initialisation of the tree", cxxopts::value(random_init))
            ("move_probs","The vector of move probabilities",cxxopts::value(move_probs))
            ;

    auto result = options.parse(argc, argv);

    if (not result.count("move_probs"))
    {
        move_probs = {0.0f,1.0f,0.0f,1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f};
        //-------------pr--w-pr--sw--w-sw---ar----w-ar--id---w-id---cs---w-cs--geno---od--
    }
    if (result.count("region_sizes_file"))
    {
        region_sizes_file = result["region_sizes_file"].as<string>();
    }

    if (not result.count("d_matrix_file"))
    {
        cerr << "the D matrix file is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("n_bins"))
    {
        cerr << "the number of bins is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (result.count("n_regions"))
    {
        // used only for naming the output
        n_regions_initial = result["n_regions"].as<size_t>();
    }
    if (not result.count("n_cells"))
    {
        cerr << "the number of cells is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (result.count("seed"))
    {
        //set a seed number for reproducibility
        SingletonRandomGenerator::get_instance(seed);
    }
    if (not result.count("random_init"))
    {
        random_init = true; // default value
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

    // run mcmc inference
    Inference mcmc(n_regions, ploidy, verbosity);

    if (random_init)
    {
        try {
            mcmc.random_initialize(n_nodes, n_regions, 10000); // creates a random tree
        }catch (const std::runtime_error& e)
        {
            std::cerr << " a runtime error was caught during the random tree initialize function, with message '"
                      << e.what() << '\'' << std::endl;
            return EXIT_FAILURE; // reject the move
        }
    }

    bool learn_nu = static_cast<bool>(move_probs[11]); // if move 11 is probable

    if(learn_nu)
    {
        if(result.count("nu"))
        {
            mcmc.t.nu = mcmc.t_prime.nu = nu;
            std::cout<<"Nu is given and going to be updated further by the chain"<<std::endl;
        }
        else
        {
            std::cout<<"Nu is not given and going to be learned from the data"<<std::endl;
        }
    }
    else
    {
        if(result.count("nu"))
        {
            mcmc.t.nu = mcmc.t_prime.nu = nu;
            std::cout<<"Nu is given and not going to be changed"<<std::endl;
        }
        else
        {
            //non-overdispersed version
            is_overdispersed = 0;
            std::cout<<"Non-overdispersed tree setting"<<std::endl;
        }
    }

    mcmc.compute_t_table(d_regions,region_sizes);
    mcmc.compute_t_od_scores(d_regions, region_sizes);

    mcmc.update_t_prime(); // set t_prime to t

    // Get starting timepoint
    auto start = high_resolution_clock::now();
    mcmc.infer_mcmc(d_regions, region_sizes, move_probs, n_iters, size_limit);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by infer_mcmc function: "
         << duration.count() << " microseconds" << endl;

    vector<vector<int>> inferred_cnvs = mcmc.assign_cells_to_nodes(d_regions, region_sizes); // returns the inferred CNVs

    vector<vector<int>> inferred_cnvs_bins = Utils::regions_to_bins_cnvs(inferred_cnvs, region_sizes);

    // write the inferred(best) tree
    std::ofstream tree_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions_initial) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_tree_inferred" + ".txt");
    tree_file << mcmc.best_tree;


    // write the inferred CNVs
    std::ofstream inferred_cnvs_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions_initial) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_inferred_cnvs" + ".txt");
    for (auto const &v1: inferred_cnvs_bins) {
        for (auto const &v2: v1)
            inferred_cnvs_file << v2 << ' ';
        inferred_cnvs_file << endl;
    }

    return EXIT_SUCCESS;
}
