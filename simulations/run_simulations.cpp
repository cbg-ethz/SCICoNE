//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#include <iostream>
#include <string>
#include "Simulation.h"
#include <cxxopts.hpp>

using namespace std;



int main(int argc, char* argv[]) {

    // the default values
    int n_regions = 50;
    int n_nodes = 50;
    double lambda_r = 0.1;
    double lambda_c = 0.2;
    int n_cells = 500;
    int n_bins = 10000;
    int n_reads = 10000;
    int n_iters = 5000;
    int n_repetitions = 1;
    int max_region_size = 10;
    int ploidy = 2;
    int verbosity = 0;
    int seed = 0;
    string f_name_postfix = "";


    cxxopts::Options options("Mcmc simulations", "simulates cnv values, infers them and benchmarks");
    options.add_options()
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("n_cells", "Number of cells", cxxopts::value(n_cells))
            ("n_nodes", "Number of nodes of the tree", cxxopts::value(n_nodes))
            ("n_regions", "Number of regions", cxxopts::value(n_regions))
            ("n_iters", "Number of iterations", cxxopts::value(n_iters))
            ("n_rep", "Number of repetitions", cxxopts::value(n_repetitions))
            ("n_reads", "Number of reads per cell", cxxopts::value(n_reads))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("verbosity", "verbosity", cxxopts::value(verbosity))
            ("seed", "seed", cxxopts::value(seed))
            ("postfix", "postfix", cxxopts::value(f_name_postfix));

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
    if (result.count("n_rep")) {
        n_repetitions = result["n_rep"].as<int>();
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


    // create the D matrix
//    vector<vector<double>> D(n_cells, vector<double>(n_regions)); // initialize with the default value
//
//    vector<int> region_sizes(n_regions); // sampling the region sizes

    Simulation sim(n_regions, n_bins, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy,
                   verbosity);
    //double delta_random_init = sim.random_cnvs_inference();

    sim.sample_region_sizes(n_bins, 1);
    sim.simulate_count_matrix(false, verbosity);
    sim.split_regions_to_bins();



    // write the region sizes
    std::ofstream region_sizes_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_region_sizes.txt");
    for (const auto &e : sim.region_sizes) region_sizes_file << e << "\n";

    // write the ground truth
    std::ofstream ground_truth_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_ground_truth.txt");
    for (auto const &v1: sim.ground_truth) {
        for (auto const &v2: v1)
            ground_truth_file << v2 << ' ';
        ground_truth_file << '\n';
    }

    // write the D matrix
    std::ofstream D_mat_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_d_mat.txt");
    for (auto const &v1: sim.D) {
        for (auto const &v2: v1)
            D_mat_file << v2 << ' ';
        D_mat_file << '\n';
    }

    // write the tree that generated the data
    std::ofstream tree_file("./"+ to_string(n_nodes)+ "nodes_" + to_string(n_regions) + "regions_" + to_string(n_reads) + "reads_"+f_name_postfix+"_tree.txt");
    tree_file << sim.tree;



//    sim.infer_cnvs(n_iters, verbosity); // n_iters: 50000
//
//    // compute the Frobenius avg. of the difference of the inferred CNVs and the ground truth
//    double delta = MathOp::frobenius_avg(sim.inferred_cnvs, sim.ground_truth);
//    cout << "delta from our method: " << delta << endl;
////
//    // cout << "delta from random method: " << delta_random_init << endl;
//
//    sim.write_d_vector(f_name_postfix);

    //simulate(n_regions, n_nodes, lambda_r, lambda_c, n_cells, n_reads, D, region_sizes); // initializes D and region_sizes

    // initialize the tree and infer the CNV profiles
    return EXIT_SUCCESS;

}
