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

    int n_cells = mat.size();



    // compute the AIC scores

    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);

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

    int n_iters = 5000; // the default value is 10000 iterations.
    int n_cells;
    int n_bins = 10000;
    int ploidy = 2;
    int verbosity = 0;
    int seed = 0;
    int window_size = 5;
    string f_name_postfix = "";
    string region_sizes_file = "";
    string d_matrix_file = "";
    bool to_segment = true; // if true then segmentation occurs


    int n_regions;

    // random tree parameters
    int n_nodes = 50;
    lambda_r = 0.1;
    lambda_c = 0.2;

    int max_region_size = 10;

    int n_reads = -1; // -1 means not specified




    cxxopts::Options options("Single cell CNV inference", "finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("verbosity", "verbosity", cxxopts::value(verbosity))
            ("seed", "seed", cxxopts::value(seed))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ("print_precision", "the precision of the score printing", cxxopts::value(print_precision))
            // random tree parameters
            ("n_nodes","the number of nodes in the random initialised tree", cxxopts::value(n_nodes))
            ("lambda_r","lambda param for the poisson that generates the number of regions", cxxopts::value(lambda_r))
            ("lambda_c","lambda param for the poisson that generates the copy number state of a region", cxxopts::value(lambda_c))
            ("n_reads","the number of reads per cell", cxxopts::value(n_reads))
            ("max_region_size","the maximum size that a region can have", cxxopts::value(max_region_size));

    auto result = options.parse(argc, argv);

    if (result.count("region_sizes_file"))
    {
        region_sizes_file = result["region_sizes_file"].as<string>();
        to_segment = false;
    }

    if (result.count("d_matrix_file"))
    {
        d_matrix_file = result["d_matrix_file"].as<string>();
    }
    else
    {
        cerr << "the D matrix file is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (result.count("n_bins"))
    {
        n_bins = result["n_bins"].as<int>();
    }
    else
    {
        cerr << "the number of bins is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (result.count("n_reads"))
    {
        n_reads = result["n_reads"].as<int>();
    }

    if (result.count("n_cells"))
    {
        n_cells = result["n_cells"].as<int>();
    }
    else
    {
        cerr << "the number of cells is not provided."<<endl;
        return EXIT_FAILURE;
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
    if (result.count("seed"))
    {
        seed = result["seed"].as<int>();
        //set a seed number for reproducibility
        SingletonRandomGenerator::get_generator(seed);
    }
    if (result.count("postfix")) {
        f_name_postfix = result["postfix"].as<string>();
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



    // set the globals
    print_precision = 16;
    lambda_s = 0.5;


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
        int _n_regions = n_bins; // assume n_regions=n_bins for the null model
        assert(n_reads != -1); // n_reads is needed for the null model
        Simulation sim_null(_n_regions, n_bins, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy,
                            verbosity);
        sim_null.sample_region_sizes(n_bins, 1);
        sim_null.simulate_count_matrix(true, verbosity);
        sim_null.split_regions_to_bins();
        vector<double> sp_null = breakpoint_detection(sim_null.D, window_size);

//        std::ofstream output_file0("./"+ to_string(n_regions) + "regions_" + to_string(n_nodes) + "nodes_"+"sim_sp_null.txt");
//        for (const auto &e : sp_null) output_file0 << e << endl;

        sp_null = dsp.crop(sp_null, window_size); // crop the window sizes



        // find the threshold on null signal
        dsp.median_normalise(sp_null);
        double std = MathOp::st_deviation(sp_null);
        double threshold = (1+9*std);



        // create s_p from the input data
        vector<double> s_p = breakpoint_detection(d_bins, window_size);

//        std::ofstream output_file1("./"+ to_string(n_regions) + "regions_" + to_string(n_nodes) + "nodes_"+"sim_sp.txt");
//        for (const auto &e : s_p) output_file1 << e << endl;

        vector<double> sp_cropped = dsp.crop(s_p, window_size);

        int lb = 0;
        int ub = sp_cropped.size()-1;

        vector<int> all_max_ids;
        deque<pair<int,int>> wait_list;
        wait_list.push_back({lb,ub}); // initial boundries

        while(!wait_list.empty())
        {
            pair<int,int> elem = wait_list.back();
            wait_list.pop_back();

            if (elem.second - elem.first > 2*window_size)
            {
                int max_idx = dsp.find_highest_peak(sp_cropped, elem.first, elem.second, threshold);
                if (max_idx != -1)
                {
                    all_max_ids.push_back(max_idx);
                    wait_list.push_back({elem.first,max_idx});
                    wait_list.push_back({max_idx,elem.second});
                }
            }
        }

        std::sort(all_max_ids.begin(), all_max_ids.end());

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
        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 50000); // creates a random tree
    }catch (const std::runtime_error& e)
    {
        std::cerr << " a runtime error was caught during the random tree initialize function, with message '"
                  << e.what() << "'\n";
        return EXIT_FAILURE; // reject the move
    }


    mcmc.compute_t_table(d_regions,region_sizes);

    mcmc.infer_mcmc(d_regions, region_sizes, move_probs, n_iters);
    vector<vector<int>> inferred_cnvs = mcmc.assign_cells_to_nodes(d_regions, region_sizes); // returns the inferred CNVs

    vector<vector<int>> inferred_cnvs_bins = Utils::regions_to_bins_cnvs(inferred_cnvs, region_sizes);

    string segmented_posfix = "";

    if (to_segment)
    {
        segmented_posfix = "_segmented";
    }


    // write the inferred(best) tree
    std::ofstream tree_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_tree_inferred" + segmented_posfix + ".txt");
    tree_file << mcmc.best_tree;


    // write the inferred CNVs
    std::ofstream inferred_cnvs_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_inferred_cnvs" + segmented_posfix + ".txt");
    for (auto const &v1: inferred_cnvs_bins) {
        for (auto const &v2: v1)
            inferred_cnvs_file << v2 << ' ';
        inferred_cnvs_file << endl;
    }



//
//    Simulation sim0(n_regions, n_bins, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy, verbosity);
//    Simulation sim1(n_regions, n_bins, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy, verbosity);
//
//
////    sim0.sample_region_sizes(n_bins);
////    sim0.simulate_count_matrix(true, verbosity);
////    sim0.split_regions_to_bins();
//
//    sim1.region_sizes = sim0.region_sizes; // use the same region sizes
//    sim1.simulate_count_matrix(false, verbosity);
//    sim1.split_regions_to_bins();

//    vector<long double> s_p0 = breakpoint_detection(sim0.D);
//    vector<long double> s_p1 = breakpoint_detection(sim1.D);
//
//
//    std::ofstream output_file0("./"+ to_string(n_regions) + "regions_" + to_string(n_nodes) + "nodes_"+"sim_sp0.txt");
//    std::ofstream output_file1("./"+ to_string(n_regions) + "regions_" + to_string(n_nodes) + "nodes_"+"sim_sp1.txt");
//    for (const auto &e : s_p0) output_file0 << e << "\n";
//    for (const auto &e : s_p1) output_file1 << e << "\n";





//    // parse input, using the fill constructor
//    vector<vector<double>> mat(n_cells, vector<double>(n_bins));
//    Utils::read_counts(mat, "../input_data/CCGL1ANXX_1_chr1_norm_counts.tsv");
//
//
//    vector<long double> s_p = breakpoint_detection(mat);
//
//
//
//    vector<bool> is_breakpoint(n_bins); // boolean mask for the breakpoints
//
//    double breakpoint_threshold = 150.0;
//
//    for (int j = 0; j < s_p.size(); ++j) {
//        is_breakpoint[j] = (s_p[j] > breakpoint_threshold);
//    }
//
//    std::ofstream output_file("./breast_tissue_E_2k_s_p_window_size_5.txt"); // TODO: build this string dynamicly and have a param to write this file, perhaps verbosity=2
//    for (const auto &e : s_p) output_file << e << "\n";
//
//
//    // create D, r and N matrices
//    vector<vector<double>> D_real; // the D matrix created from the real data
//
//    for (vector<double> &row : mat) // Segmentation is performed here
//    {
//        vector<double> new_row;
//        double cum_val = 0.0;
//        for (int i = 0; i < row.size(); ++i)
//        {
//            cum_val += row[i];
//            if (is_breakpoint[i]) // breakpoint
//            {
//                new_row.push_back(cum_val);
//                cum_val = 0.0;
//            }
//        }
//        if (cum_val != 0.0) // consider the last bin as well
//            new_row.push_back(cum_val);
//        D_real.push_back(new_row);
//    }
//
//    vector<int> r_real; // the r matrix from the real data
//    int region_size = 0;
//    for (bool elem : is_breakpoint)
//    {
//        region_size++;
//        if (elem) //breakpoint
//        {
//            r_real.push_back(region_size);
//            region_size = 0;
//        }
//    }
//    if (region_size != 0)
//        r_real.push_back(region_size);
//
//    int sum_r_real = std::accumulate(r_real.rbegin(), r_real.rend(), 0);
//    assert(sum_r_real == is_breakpoint.size());
//
//    vector<double> N_real; // the N matrix
//
//    for (vector<double> &row : D_real)
//    {
//        double sum_row = std::accumulate(row.begin(),row.end(),0.0);
//        N_real.push_back(sum_row);
//    }
//
//
//    // move probabilities
//    vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f};
//
//    Inference mcmc(r_real.size(), ploidy, verbosity);
//
////    mcmc.initialize_worked_example();
//    u_int n_nodes = 15;
//    double lambda_r = 0.1;
//    double lambda_c = 0.2;
//    int n_regions = D_real[0].size()-1;
//    try {
//        mcmc.random_initialize(n_nodes, n_regions, lambda_r, lambda_c, 10000); // creates a random tree
//    }catch (const std::runtime_error& e)
//    {
//        std::cerr << " a runtime error was caught during the random tree initialize function, with message '"
//                  << e.what() << "'\n";
//        return EXIT_FAILURE; // reject the move
//    }
//
//
//    mcmc.compute_t_table(D_real,r_real);
//
//    mcmc.infer_mcmc(D_real, r_real, move_probs, n_iters);
//    mcmc.assign_cells_to_nodes(D_real, r_real); // returns the inferred CNVs
//
//    if (verbosity > 1)
//        mcmc.write_best_tree();


    return EXIT_SUCCESS;
}
