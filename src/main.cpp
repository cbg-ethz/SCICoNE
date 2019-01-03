#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <queue>
#include <string>

#include <iterator>
#include <algorithm>

#include "MathOp.h"
#include "Tree.h"
#include "Inference.h"

#include <chrono> // for measuring the execution time

#include <cxxopts.hpp>

#include "Simulation.h"

#include "globals.cpp"

#include "SignalProcessing.h"


// globals
int print_precision;
int copy_number_limit;
double lambda_s;
double lambda_r;
double lambda_c;

// endof globals

using namespace std;
using namespace std::chrono;


vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size = 5)
{
    /*
     * Performs the breakpoint detection
     * */

    size_t n_cells = mat.size();



    // compute the AIC scores

    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);

    vector<vector<double>> sigma;

    size_t n_breakpoints = aic_vec.size();
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

    vector<double> s_p;

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

        double sp_val = sp_denom-sp_num;

        s_p.push_back(sp_val);

    }

    return s_p;

}


int main( int argc, char* argv[]) {

    int n_iters = 5000; // the default value is 5000 iterations.
    int n_cells;
    int n_bins = 10000;
    int ploidy = 2;
    int verbosity = 0;
    int seed = 0;
    int window_size = 10;
    string f_name_postfix;
    string region_sizes_file;
    string d_matrix_file;
    bool to_segment = true; // if true then segmentation occurs

    int size_limit = -1;

    size_t n_regions;
    size_t n_regions_initial = 0; // used for naming the output for I/O workflow purposes

    // random tree parameters
    int n_nodes = 50;
    lambda_r = 0.1;
    lambda_c = 0.2;

    int max_region_size = 10;

    int n_reads = -1; // -1 means not specified

    // set the globals
    print_precision = 16;
    lambda_s = 0.5;
    copy_number_limit = 5;


    cxxopts::Options options("Single cell CNV inference", "finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_regions", "Number of regions to be contained in the output file", cxxopts::value(n_regions))
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
            ("max_region_size","the maximum size that a region can have", cxxopts::value(max_region_size))
            ("window_size", "the size of the window used in breakpoint detection", cxxopts::value(window_size))
            ;

    auto result = options.parse(argc, argv);

    if (result.count("region_sizes_file"))
    {
        region_sizes_file = result["region_sizes_file"].as<string>();
        to_segment = false;
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


    // read the input d_matrix
    vector<vector<double>> d_bins(n_cells, vector<double>(n_bins));
    Utils::read_counts(d_bins, d_matrix_file);

    // read the region_sizes file
    vector<int> region_sizes;
    vector<vector<double>> d_regions;

    if (to_segment)
    {

        // perform segmentation and peak detection, then define the region sizes
        SignalProcessing dsp;

        // create an sp_null signal
        size_t _n_regions = n_bins; // assume n_regions=n_bins for the null model
        // compute n_reads by summing the cell (a row of the matrix) up
        n_reads = static_cast<int>(accumulate(d_bins[0].begin(), d_bins[0].end(), 0.0));

        assert(n_reads != -1); // n_reads is needed for the null model

        vector<double> s_p = breakpoint_detection(d_bins, window_size);
//        std::ofstream output_file1("./"+ to_string(window_size) + "ws" + to_string(n_regions) + "regions_" + to_string(n_nodes) + "nodes_"+"sim_sp.txt");
//        for (const auto &e : s_p) output_file1 << e << endl;
//        cout <<"sp created"<<endl;

        vector<double> sp_cropped = dsp.crop(s_p, window_size);

        // median normalise sp_cropped
        dsp.median_normalise(sp_cropped);
        vector<double> sp_cropped_copy(sp_cropped); // the copy sp vector that'll contain the NaN values

        int lb = 0;
        size_t ub = sp_cropped.size();

        vector<int> all_max_ids;
        map<double, pair<unsigned,unsigned>> q_map; // a map that serves as a queue
        // use a map because you want the values to be always sorted
        // start with the 0.0 value, it will be removed from the map before other values are inserted
        q_map.emplace(0.0 , std::make_pair(lb,ub)); // initial boundries

        while(!q_map.empty())
        {
            pair<unsigned,unsigned> elem = q_map.rbegin()->second; // use the rbegin to get the largest val
            q_map.erase(q_map.rbegin()->first); // remove the largest val by it's key

            if (elem.second - elem.first > window_size)
            {
                int threshold_coef = 4;
                int max_idx = dsp.evaluate_peak(sp_cropped, sp_cropped_copy, elem.first, elem.second,
                                                threshold_coef);
                if (max_idx != -1) // TODO: instead of checking for -1 use proper exception handling
                {

                    // replace the nearby bins by nan
                    int start_idx, stop_idx;
                    start_idx = max_idx - window_size;
                    stop_idx = max_idx + window_size;

                    // check the boundries
                    if (start_idx < 0)
                        start_idx = 0;
                    if (stop_idx > sp_cropped.size())
                        stop_idx = sp_cropped.size();
                    // set the nearby bins to nan
                    for (int i = start_idx; i < stop_idx; ++i) {
                        sp_cropped_copy[i] = std::nan("");
                    }


                    // compute the left and right medians
                    vector<double> left_vec(sp_cropped.begin() + elem.first, sp_cropped.begin() + max_idx);
                    vector<double> right_vec(sp_cropped.begin() + max_idx + 1, sp_cropped.begin() + elem.second);

                    if (!left_vec.empty())
                    {
                        double median_left = MathOp::median(left_vec);
                        // normalise bins on the left by left median
                        for (int i = elem.first; i < max_idx; ++i) {
                            sp_cropped[i] /= median_left;
                        }
                        if (max_idx - elem.first > window_size)
                        {
                            int max_left_idx = dsp.find_highest_peak(sp_cropped, elem.first, max_idx);
                            q_map.emplace(sp_cropped[max_left_idx] , std::make_pair(elem.first,max_idx));
                        }
                    }
                    if (!right_vec.empty())
                    {
                        double median_right = MathOp::median(right_vec);
                        // normalise bins on the right by right median
                        for (int j = max_idx + 1; j < elem.second; ++j) {
                            sp_cropped[j] /= median_right;
                        }
                        if (elem.second - (max_idx+1) > window_size)
                        {
                            int max_right_idx = dsp.find_highest_peak(sp_cropped, max_idx+1, elem.second);
                            q_map.emplace(sp_cropped[max_right_idx] , std::make_pair(max_idx + 1,elem.second));
                        }
                    }

                    all_max_ids.push_back(max_idx);

                }
            }
        }

        std::sort(all_max_ids.begin(), all_max_ids.end());

        for (int k = 0; k < all_max_ids.size(); ++k) {
            cout << all_max_ids[k] + window_size << endl;
        }

        vector<bool> sp_breakpoints(sp_cropped.size());

        for (int i = 0; i < all_max_ids.size(); ++i) {
            sp_breakpoints[all_max_ids[i]] = true;
        }

        // add back the cropped window sizes
        for (int j = 0; j < window_size; ++j) {
            sp_breakpoints.push_back(false);
            sp_breakpoints.insert(sp_breakpoints.begin(), false);
        }

        region_sizes = dsp.create_region_sizes(sp_breakpoints);
        int sum_region_sizes = accumulate( region_sizes.begin(), region_sizes.end(), 0);
        assert(sum_region_sizes == n_bins);


    }
    else
    {
        Utils::read_vector(region_sizes, region_sizes_file);
    }
    // Merge the bins into regions
    d_regions = Utils::condense_matrix(d_bins, region_sizes);
    n_regions = region_sizes.size();

    // run mcmc inference

    // move probabilities
    vector<float> move_probs = {1.0f,0.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f};
    //-------------------------------w-pr------------------------------w-id--------w-cs-------

    Inference mcmc(n_regions, ploidy, verbosity);

    try {
        mcmc.random_initialize(n_nodes, n_regions, 10000); // creates a random tree
    }catch (const std::runtime_error& e)
    {
        std::cerr << " a runtime error was caught during the random tree initialize function, with message '"
                  << e.what() << "'\n";
        return EXIT_FAILURE; // reject the move
    }

    mcmc.compute_t_table(d_regions,region_sizes);

    // Get starting timepoint
    auto start = high_resolution_clock::now();
    mcmc.infer_mcmc(d_regions, region_sizes, move_probs, n_iters, size_limit);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by infer_mcmc function: "
         << duration.count() << " microseconds" << endl;

    vector<vector<int>> inferred_cnvs = mcmc.assign_cells_to_nodes(d_regions, region_sizes); // returns the inferred CNVs

    vector<vector<int>> inferred_cnvs_bins = Utils::regions_to_bins_cnvs(inferred_cnvs, region_sizes);

    string segmented_posfix = "";

    if (to_segment)
    {
        segmented_posfix = "_segmented";
    }


    // write the inferred(best) tree
    std::ofstream tree_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions_initial) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_tree_inferred" + segmented_posfix + ".txt");
    tree_file << mcmc.best_tree;


    // write the inferred CNVs
    std::ofstream inferred_cnvs_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions_initial) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_inferred_cnvs" + segmented_posfix + ".txt");
    for (auto const &v1: inferred_cnvs_bins) {
        for (auto const &v2: v1)
            inferred_cnvs_file << v2 << ' ';
        inferred_cnvs_file << endl;
    }

    return EXIT_SUCCESS;
}
