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
double cf;
double c_penalise;
unsigned is_overdispersed;
string f_name_posfix;
int verbosity;
double eta;
double beta;

// endof globals

using namespace std;
using namespace std::chrono;


int main( int argc, char* argv[]) {

    int n_iters = 5000; // the default value is 5000 iterations.
    int n_cells;
    int ploidy = 2;
    verbosity = 0;
    int seed = 0;
    string region_sizes_file;
    string d_matrix_file;
    string cluster_sizes_file;
    string tree_file;

    unsigned size_limit = std::numeric_limits<unsigned>::max();

    size_t n_regions;

    // random tree parameters
    int n_nodes = 3;
    lambda_r = 0.1;
    lambda_c = 0.2;
    cf = 0.0;
    c_penalise = 10.0;
    is_overdispersed = 1;
    eta = 1e-4;
    // Penalty on number of genotypes
    beta = 0.0;

    // set the globals
    print_precision = 15;
    lambda_s = 0.5;
    copy_number_limit = 5;

    // overdispersion parameter
    double nu = 1.0;
    // learning rate params
    double alpha = 0.0;
    double gamma = 1.0;
    bool random_init = true;

    // move probabilities
    vector<float> move_probs;
    string move_probs_str = "0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.01,0.1,0.01,1.0,0.01";

    // scoring mode
    bool max_scoring = true;

    // region neutral states
    string region_neutral_states_file;

    cxxopts::Options options("Single cell CNV inference", "finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ',', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("cluster_sizes_file", "Path to the file containing the cluster sizes, each line contains one cluster size", cxxopts::value(cluster_sizes_file))
            ("tree_file", "the tree file to load", cxxopts::value(tree_file))
            ("n_regions", "Number of regions in the input matrix", cxxopts::value(n_regions))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("verbosity", "verbosity of the program", cxxopts::value(verbosity))
            ("seed", "seed", cxxopts::value(seed))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_posfix))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            ("size_limit", "the limitation on the max size of the tree", cxxopts::value(size_limit))
            ("copy_number_limit", "the maximum copy number profile one bin or region can have", cxxopts::value(copy_number_limit))
            // random tree parameters
            ("n_nodes","the number of nodes in the random initialised tree", cxxopts::value(n_nodes))
            ("lambda_r","lambda param for the poisson that generates the number of regions", cxxopts::value(lambda_r))
            ("lambda_c","lambda param for the poisson that generates the copy number state of a region", cxxopts::value(lambda_c))
            ("cf", "cluster fraction variable between 0 and 1 to affect the tree prior coefficient", cxxopts::value(cf))
            ("c_penalise","term that penalises trees containing cancelling events to be added to tree event prior",cxxopts::value(c_penalise))
            // ("is_overdispersed", "multinomial or dirichlet multinomial in the likelihood", cxxopts::value(is_overdispersed))
            ("nu","nu parameter, the overdispersion variable",cxxopts::value(nu))
            ("alpha","alpha parameter, the constant in the learning rate",cxxopts::value(alpha))
            ("gamma","gamma parameter, the initial learning rate value",cxxopts::value(gamma))
            ("random_init","Boolean parameter to enable random initialisation of the tree", cxxopts::value(random_init))
            ("move_probs","The vector of move probabilities",cxxopts::value(move_probs)->default_value(move_probs_str))
            ("max_scoring","Boolean parameter to decide whether to take the maximum score or to marginalize over all assignments during inference",cxxopts::value<bool>(max_scoring)->default_value("true"))
            ("region_neutral_states_file", "Path to the file containing the neutral state of each region to use as the root of the tree", cxxopts::value(region_neutral_states_file))
            ("beta", "term that penalises the number of genotypes in the tree", cxxopts::value(beta))
            ;

    auto result = options.parse(argc, argv);

    if (not result.count("move_probs")) {
        move_probs = {0.0f,1.0f,0.0f,1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.01f, 0.1f, 0.01f, 1.0f, 0.01f};
        //-------------pr--w-pr--sw--w-sw---ar---w-ar---id---w-id---cs---w-cs--geno---es----ca-----od----dl---
    }
    if (result.count("region_sizes_file"))
    {
        region_sizes_file = result["region_sizes_file"].as<string>();
    }
    else
    {
        cerr << "the region sizes file is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("d_matrix_file"))
    {
        cerr << "the D matrix file is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("n_regions"))
    {
        cerr << "the number of regions is not provided."<<endl;
        return EXIT_FAILURE;
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
    if (result.count("tree_file")) // undo random_init
        random_init = false;

    if (verbosity > 0) {
        if(result.count("alpha"))
            std::cout<<"the learning rate alpha constant is specified to be " << alpha <<std::endl;
        else
            std::cout<<"the learning rate alpha constant is " << alpha << " since it is not specified" <<std::endl;
        if(result.count("gamma"))
            std::cout<<"initial gamma is specified to be " << gamma <<std::endl;
        else
            std::cout<<"gamma value is not specified, the default value is: " << gamma <<std::endl;
        if (result.count("cf"))
            std::cout<<"Cluster fraction in tree prior is: " << cf << std::endl;
        if (result.count("beta"))
            std::cout<<"beta in tree prior is: " << beta << std::endl;
    }

    if (verbosity > 0)
        std::cout << "Reading the input matrix..." << std::endl;
    vector<vector<double>> d_regions(n_cells, vector<double>(n_regions));
    Utils::read_counts(d_regions, d_matrix_file);

    if (verbosity > 0)
        std::cout << "Reading the region_sizes file..." << std::endl;
    vector<int> region_sizes;
    Utils::read_vector(region_sizes, region_sizes_file);

    vector<int> cluster_sizes;
    bool read_cluster_sizes = false;
    if (result.count("cluster_sizes_file")) {
      if (cluster_sizes_file.compare("") != 0) {
        if (verbosity > 0)
            std::cout << "Reading the cluster_sizes file..." << std::endl;
        Utils::read_vector(cluster_sizes, cluster_sizes_file);
        read_cluster_sizes = true;
      }
    }
    if (not read_cluster_sizes) {
      cluster_sizes = std::vector<int>(n_cells, 1);
      if (n_cells < 20)
        std::cout << "Warning: there are only " << n_cells <<  " observations. If these are clusters, the cluster_sizes_file parameter should be specified for accurate tree scoring." << std::endl;
    }

    if (not max_scoring) {
       if (verbosity > 0)
           std::cout << "Will perform sum scoring." << std::endl;
       cf = 1.0;
	     max_scoring = false;
       // ES and CA moves are only available in max_scoring mode
       move_probs[11] = 0.0f;
       move_probs[12] = 0.0f;
   } else {
       if (verbosity > 0)
           std::cout << "Will perform maximum scoring." << std::endl;
       move_probs.back() = 0.0f; // no need to prune tree;
   }

    vector<int> region_neutral_states;
    bool read_region_neutral_states = false;
    if (result.count("region_neutral_states_file")) {
      if (region_neutral_states_file.compare("") != 0) {
        if (verbosity > 0)
            std::cout << "Reading the region_neutral_states file..." << std::endl;
        Utils::read_vector(region_neutral_states, region_neutral_states_file);
        read_region_neutral_states = true;
      }
    }
    if (not read_region_neutral_states) {
      if (verbosity > 0)
        std::cout << "Assuming root to have copy number state " << ploidy << " in all regions" << std::endl;
      region_neutral_states = std::vector<int>(n_regions, ploidy);
    }

    // run mcmc inference
    Inference mcmc(n_regions, n_cells, region_neutral_states, ploidy, verbosity, max_scoring);

    if (random_init)
    {
        int max_iters = 10000;
        if (verbosity > 0)
            std::cout<<"Trying to randomly initialise a valid tree with " << max_iters << " iterations." << std::endl;
        try {
            mcmc.random_initialize(n_nodes, n_regions, max_iters); // creates a random tree
        }catch (const std::runtime_error& e)
        {
            std::cerr << "A runtime error was caught during the random tree initialize function, with message '"
                      << e.what() << '\'' << std::endl;
            return EXIT_FAILURE; // reject the move
        }
    }
    if (result.count("tree_file")) { // starting tree is specified
        mcmc.initialize_from_file(tree_file);
        nu = mcmc.t.nu;
    }

    bool learn_nu = static_cast<bool>(move_probs[13]); // if move is probable

    if(learn_nu)
    {
        if(result.count("nu"))
        {
            mcmc.t.nu = mcmc.t_prime.nu = nu;
            if (verbosity > 0)
                std::cout<<"Nu is given and going to be updated further by the chain"<<std::endl;
        }
        else
        {
            if (verbosity > 0)
                std::cout<<"Nu is not given and going to be learned from the data"<<std::endl;
        }
    }
    else
    {
        if(result.count("nu"))
        {
            mcmc.t.nu = mcmc.t_prime.nu = nu;
            if (verbosity > 0)
                std::cout<<"Nu is given and not going to be changed"<<std::endl;
        }
        else
        {
            //non-overdispersed version
            is_overdispersed = 0;
            if (verbosity > 0)
                std::cout<<"Non-overdispersed tree setting"<<std::endl;
        }
    }

    if (verbosity > 0)
        std::cout << "Computing the t table and overdispersion scores for the initial tree..." << std::endl;
    mcmc.compute_t_table(d_regions,region_sizes,cluster_sizes);
    mcmc.compute_t_od_scores(d_regions, region_sizes,cluster_sizes);

    mcmc.update_t_prime(); // set t_prime to t

    // Get starting timepoint
    auto start = high_resolution_clock::now();

    if (verbosity > 0)
        std::cout << "Inferring the tree using MCMC with " << n_iters << " iterations..." << std::endl;
    mcmc.infer_mcmc(d_regions, region_sizes, move_probs, n_iters, size_limit, alpha, gamma, cluster_sizes);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);

    if (verbosity > 0)
        std::cout << "Time taken by infer_mcmc function: "
             << duration.count() << " microseconds" << std::endl;

    if (verbosity > 0)
        std::cout << "Writing the inferred tree and cnvs..." <<std::endl;

    vector<vector<int>> inferred_cnvs = mcmc.assign_cells_to_nodes(d_regions, region_sizes, cluster_sizes); // returns the inferred CNVs

    vector<vector<int>> inferred_cnvs_bins = Utils::regions_to_bins_cnvs(inferred_cnvs, region_sizes);

    std::vector<std::map<int, double>> attachment_scores = mcmc.t_scores; // cell-by-node score matrix

    // write the inferred(best) tree
    std::ofstream tree_output_file("./" + f_name_posfix + "_tree_inferred" + ".txt");
    tree_output_file << mcmc.best_tree;


    // write the inferred CNVs
    std::ofstream inferred_cnvs_file("./" + f_name_posfix + "_inferred_cnvs" + ".csv");
    for (auto const &v1: inferred_cnvs_bins) {
        for (size_t i = 0; i < v1.size(); i++)
        {
            if (i == v1.size()-1) // the last element
                inferred_cnvs_file << v1[i];
            else // add comma
                inferred_cnvs_file << v1[i] << ',';
        }
        inferred_cnvs_file << endl;
    }

    if (verbosity > 1)
    {
      // write score matrix
      std::ofstream attachment_scores_file("./" + f_name_posfix + "_attachment_scores" + ".csv");
      for (int j=0;j<n_cells;j++) {
        map<int,double>::const_iterator map_iter;
        int i = 0;
        for(map_iter=attachment_scores[j].begin(); map_iter!=attachment_scores[j].end(); ++map_iter)
        {
          if (i == attachment_scores[j].size()-1) // the last element
              attachment_scores_file << map_iter->second;
          else // add comma
              attachment_scores_file << map_iter->second << ',';
          i++;
        }
        attachment_scores_file << endl;
      }
    }

    if (verbosity > 0)
        std::cout << "Tree inference is successfully completed!" <<std::endl;
    return EXIT_SUCCESS;
}
