//
// Created by Tuncel  Mustafa Anil on 7/4/19.
//

#include "Inference.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include "globals.cpp"
#include <cxxopts.hpp>


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
// endof globals

struct double_iota
{
    /*
     * std::iota for double data type
     * */

    double_iota(double inc, double init_value = 0.0) : _value(init_value), _inc(inc) {}

    operator double() const { return _value; }
    double_iota& operator++() { _value += _inc; return *this; }
    double _value;
    double _inc;
};


void save_root_by_nu(int ploidy, size_t n_regions, int n_cells, const vector<vector<double>> &d_regions, const vector<int> &region_sizes, vector<int> &region_neutral_states)
{
    /*
     * Writes the overdispersed root scores for various values of nu.
     * */

    double min = 10e-6;
    double max = 0.05;
    double step = 0.001;

    std::vector<double> nu(std::size_t(((max - min + step - std::numeric_limits<double>::epsilon())) / step));
    std::iota(nu.begin(), nu.end(), double_iota(step, min));

    Tree t(ploidy,n_regions,region_neutral_states);

    std::ofstream root_nu_scores_file("root_scores_per_nu.csv");

    // header
    root_nu_scores_file << "nu score, root score" << endl;

    for (auto const& nu_val : nu)
    {
        t.nu = nu_val;

        vector<double> scores(n_cells);
        for (u_int i = 0; i < n_cells; ++i) {
            double sum_d = std::accumulate(d_regions[i].begin(), d_regions[i].end(), 0.0);
            scores[i] = t.get_od_root_score(region_sizes, sum_d, d_regions[i]);
        }
        double root_score = std::accumulate(scores.begin(), scores.end(), 0.0);

        root_nu_scores_file << nu_val << ',' << setprecision(15)  << root_score << endl;

    }


}

int main(int argc, char* argv[])
{
    string region_sizes_file;
    string d_matrix_file;
    int ploidy = 2;
    int n_cells;
    int n_bins;
    verbosity = 0;
    is_overdispersed = 1;
    eta = 1e-4;
    cf = 1.0;
    string region_neutral_states_file;

    cxxopts::Options options("Statistical tests", "for testing the statistical properties of the programme");
    options.add_options()
            ("region_sizes_file", "Path to the file containing the region sizes, each line contains one region size", cxxopts::value(region_sizes_file))
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("ploidy", "ploidy", cxxopts::value(ploidy))
            ("region_neutral_states_file", "Path to the file containing the neutral state of each region to use as the root of the tree", cxxopts::value(region_neutral_states_file));

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
    size_t n_regions = region_sizes.size();

    vector<int> region_neutral_states;
    if (result.count("region_neutral_states_file")) {
      std::cout << "Reading the region_neutral_states file..." << std::endl;
      Utils::read_vector(region_neutral_states, region_neutral_states_file);
    }
    else {
      std::cout << "Assuming root to have copy number state " << ploidy << " in all regions" << std::endl;
      region_neutral_states = std::vector<int>(n_regions, ploidy);
    }

    save_root_by_nu(ploidy, n_regions, n_cells, d_regions, region_sizes, region_neutral_states);


}
