#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <string>

#include <iterator>


#include "MathOp.h"
#include "Tree.h"
#include "Inference.h"

#include <chrono> // for measuring the execution time

#include <cxxopts.hpp>



using namespace std;
using namespace std::chrono;


int main( int argc, char* argv[]) {

    int n_iters = 5000; // the default value is 10000 iterations.
    int n_cells = 0;
    int n_bins = 0;
    int ploidy = 2;
    int verbosity = 0;


    cxxopts::Options options("Single cell CNV inference", "finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified");
    options.add_options()
            ("n_bins", "Number of regions", cxxopts::value(n_bins))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_cells", "Number of cells", cxxopts::value(n_cells))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("verbosity", "verbosity", cxxopts::value(verbosity));

    auto result = options.parse(argc, argv);

    if (result.count("n_bins"))
    {
        n_bins = result["n_bins"].as<int>();
    }
    if (result.count("n_cells"))
    {
        n_cells = result["n_cells"].as<int>();
    }
    if (result.count("n_iters"))
    {
        n_iters = result["n_iters"].as<int>();
    }
    if (result.count("verbosity"))
    {
        verbosity = result["verbosity"].as<int>();
    }
    if (result.count("ploidy"))
    {
        ploidy = result["ploidy"].as<int>();
    }


    // set a seed number for reproducibility
    //SingletonRandomGenerator::get_generator(42);

    // parse input, using the fill constructor
    vector<vector<double>> mat(n_cells, vector<double>(n_bins));
    Utils::read_counts(mat, "../input_data/CCGL1ANXX_1_chr1_norm_counts.tsv");

    // compute the AIC scores
    u_int window_size = 5;
    auto start = std::chrono::high_resolution_clock::now(); // start the clock
    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);
    auto stop = high_resolution_clock::now();

    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "\n\nTime taken by the MathOp::likelihood_ratio function on first sample lane 6 dna data: "
         << duration.count() << " microseconds" << endl;


    vector<vector<double>> sigma;

    int n_breakpoints = aic_vec.size();
    cout <<"n_breakpoints: " << n_breakpoints << " n_cells: " << n_cells <<endl;

    for (auto &vec: aic_vec) // compute sigma matrix
    {
        auto res = MathOp::combine_scores(vec);
        sigma.push_back(res);
    }

    vector<double> log_priors;
    for (int j = 0; j < n_cells; ++j) {
        log_priors.push_back(MathOp::breakpoint_log_prior(j, n_cells,0.5));
    }


    vector<vector<long double>> log_posterior;

    for (int k = 0; k < n_breakpoints; ++k) {
        log_posterior.push_back(vector<long double>());
        for (int j = 0; j < n_cells; ++j) {
            long double val = log_priors[j] + sigma[k][j];
            log_posterior[k].push_back(val);
        }
    }

    vector<vector<long double>> posterior;
    int k_star = 4;

    vector<long double> s_p;
    vector<bool> is_breakpoint; // boolean mask for the breakpoints
    for (int l = 0; l < n_breakpoints; ++l)
    {
        posterior.push_back(vector<long double>());

        long double max_num = *max_element(log_posterior[l].begin(), log_posterior[l].begin()+k_star-1);
        long double max_denom = *max_element(log_posterior[l].begin(), log_posterior[l].end());

        for (int j = 0; j < k_star - 1; ++j) {
            long double val =exp(log_posterior[l][j] - max_num);
            posterior[l].push_back(val);
        }
        for (int k = k_star -1 ; k < log_posterior[l].size(); ++k) {
            long double val =exp(log_posterior[l][k] - max_denom);
            posterior[l].push_back(val);
        }


        long double sp_num = std::accumulate(posterior[l].begin(), posterior[l].begin()+k_star-1, 0.0);
        sp_num  = log(sp_num) + max_num;
        long double sp_denom = std::accumulate(posterior[l].begin(), posterior[l].end(), 0.0);
        sp_denom = log(sp_denom) + max_denom;

        long double sp_val = sp_denom-sp_num;

        s_p.push_back(sp_val);

        double breakpoint_threshold = 150.0;
        is_breakpoint.push_back(sp_val > breakpoint_threshold);
    }

    std::ofstream output_file("./breast_tissue_E_2k_s_p_window_size_5.txt"); // TODO: build this string dynamicly and have a param to write this file, perhaps verbosity=2
    for (const auto &e : s_p) output_file << e << "\n";


    // create D, r and N matrices
    vector<vector<double>> D_real; // the D matrix created from the real data

    for (vector<double> &row : mat)
    {
        vector<double> new_row;
        double cum_val = 0.0;
        for (int i = 0; i < row.size(); ++i)
        {
            cum_val += row[i];
            if (is_breakpoint[i]) // breakpoint
            {
                new_row.push_back(cum_val);
                cum_val = 0.0;
            }
        }
        if (cum_val != 0.0) // consider the last bin as well
            new_row.push_back(cum_val);
        D_real.push_back(new_row);
    }

    vector<int> r_real; // the r matrix from the real data
    int region_size = 0;
    for (bool elem : is_breakpoint)
    {
        region_size++;
        if (elem) //breakpoint
        {
            r_real.push_back(region_size);
            region_size = 0;
        }
    }
    if (region_size != 0)
        r_real.push_back(region_size);

    int sum_r_real = std::accumulate(r_real.rbegin(), r_real.rend(), 0);
    assert(sum_r_real == is_breakpoint.size());

    vector<double> N_real; // the N matrix

    for (vector<double> &row : D_real)
    {
        double sum_row = std::accumulate(row.begin(),row.end(),0.0);
        N_real.push_back(sum_row);
    }


    // move probabilities
    vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f};

    Inference mcmc(r_real.size(), ploidy, verbosity);

//    mcmc.initialize_worked_example();
    u_int n_nodes = 15;
    double lambda_r = 0.1;
    double lambda_c = 0.2;
    int n_regions = D_real[0].size()-1;
    try {
        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree
    }catch (const std::runtime_error& e)
    {
        std::cerr << " a runtime error was caught during the random tree initialize function, with message '"
                  << e.what() << "'\n";
        return EXIT_FAILURE; // reject the move
    }


    mcmc.compute_t_table(D_real,r_real);

    mcmc.infer_mcmc(D_real, r_real, move_probs, n_iters);
    mcmc.assign_cells_to_nodes(D_real, r_real); // returns the inferred CNVs

    if (verbosity > 1)
        mcmc.write_best_tree();


    return EXIT_SUCCESS;
}
