//
// Created by Tuncel  Mustafa Anil on 1/3/19.
//

#include <cxxopts.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

#include "SignalProcessing.h"
#include "Utils.h"
#include "globals.cpp"

// globals
int verbosity;
string f_name_posfix;
// end of globals

using namespace std;

vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size, int k_star);

int main( int argc, char* argv[]) {

    int n_cells;
    int n_bins = 10000;
    int window_size = 10;
    double threshold_coefficient = 3.0;
    f_name_posfix = "";
    string region_sizes_file;
    string d_matrix_file;
    verbosity = 0;
    int evidence_min_cells = 4;
    unsigned breakpoints_limit = 300;

    cxxopts::Options options("Breakpoint detection executable", "detects the breakpoints in the genome across all cells.");
    options.add_options()
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ',', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("min_cells", "Minimum number of cells to consider for a bin being a breakpoint", cxxopts::value(evidence_min_cells))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("window_size", "the size of the window used in breakpoint detection", cxxopts::value(window_size))
            ("threshold", "the coefficient of the breakpoint threshold", cxxopts::value(threshold_coefficient))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_posfix))
            ("verbosity", "verbosity", cxxopts::value(verbosity))
            ("bp_limit","the maximum number of breakpoints to be returned. The breakpoints get sorted and the top ones are returned",cxxopts::value(breakpoints_limit))
            ;

    auto result = options.parse(argc, argv);

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
    if (not result.count("n_cells"))
    {
        cerr << "the number of cells is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("min_cells"))
    {
        std::cout << "min_cells is not specified, " << evidence_min_cells << " cells are going to be considered" <<std::endl;
    }

    std::cout<<"Reading the input matrix..."<<std::endl;

    vector<vector<double>> d_bins(n_cells, vector<double>(n_bins));
    Utils::read_counts(d_bins, d_matrix_file);

    std::cout<<"Input matrix is read."<<std::endl;

    // create the region_sizes
    vector<int> region_sizes;
    vector<vector<double>> d_regions;

    // perform segmentation and peak detection, then define the region sizes
    SignalProcessing dsp;

    std::cout<<"Computing the probability of a region being a breakpoint..."<<std::endl;
    vector<double> s_p = dsp.breakpoint_detection(d_bins, window_size, evidence_min_cells);
    std::cout<<"Computed probabilities for all regions."<<std::endl;

    vector<double> sp_cropped = dsp.crop(s_p, window_size);

    // median normalise sp_cropped
    dsp.median_normalise(sp_cropped);

    /* if s_p contains zero, then replace those with the previous
     * Otherwise log(0) creates -INF that messes up mean, std, even median
     * Replacing zeros with a minimum value is not a good idea, because a breakpoint can be introduced by this imputation.
     */
    std::cout<<"Replacing zeroes by the previous value, lest them to be effective in the breakpoint detection."<<std::endl;
    for (int l = 0; l < sp_cropped.size(); ++l)
        if(sp_cropped[l] == 0.0)
        {
            if (l > 0)
                sp_cropped[l] = sp_cropped[l - 1];
            else // if zero is the first element
                sp_cropped[l] = 1e-8; // a small positive number
        }
    std::cout<<"Zeroes are replaced by the previous value."<<std::endl;

    vector<double> sp_cropped_copy(sp_cropped); // the copy sp vector that'll contain the NaN values

    int lb = 0;
    size_t ub = sp_cropped.size();

    vector<int> all_max_ids;
    map<double, pair<unsigned,unsigned>> q_map; // a map that serves as a queue
    // use a map because you want the values to be always sorted
    // start with the 0.0 value, it will be removed from the map before other values are inserted
    q_map.emplace(0.0 , std::make_pair(lb,ub)); // initial boundaries

    int smaller_window_size = static_cast<int>(window_size * 0.4);

    while(!q_map.empty())
    {
        pair<int,int> elem = q_map.rbegin()->second; // use the rbegin to get the largest val
        q_map.erase(q_map.rbegin()->first); // remove the largest val by it's key

        if (elem.second - elem.first > smaller_window_size)
        {
            int max_idx = dsp.evaluate_peak(sp_cropped, sp_cropped_copy, elem.first, elem.second, threshold_coefficient,
                                            window_size); // here the big window size otherwise the output ids would not match the effective_break_points (ground truth)
            if (max_idx != -1) // -1 means rejected
            {

                // replace the nearby bins by nan
                int start_idx, stop_idx;
                start_idx = max_idx - smaller_window_size;
                stop_idx = max_idx + smaller_window_size + 1; // +1 because if max_id = 100 and window_size = 4, then 101,102,103,104 must be NaN

                // check the boundries
                if (start_idx < 0)
                    start_idx = 0;
                if (stop_idx > sp_cropped.size())
                    stop_idx = sp_cropped.size();
                // set the nearby bins to nan
                for (int i = start_idx; i < stop_idx; ++i) {
                    sp_cropped_copy[i] = std::nan("");
                }

                if ((max_idx - smaller_window_size) - elem.first > smaller_window_size + 1) // +1 to make sure the size of vector is at least 2
                {
                    // compute the left median
                    vector<double> left_vec(sp_cropped.begin() + elem.first, sp_cropped.begin() + max_idx - smaller_window_size);

                    double median_left = MathOp::median(left_vec);
                    // normalise bins on the left by left median
                    for (int i = elem.first; i < max_idx - smaller_window_size; ++i) {
                        sp_cropped[i] /= median_left;
                    }
                    try
                    {
                        int max_left_idx = dsp.find_highest_peak(sp_cropped, elem.first, max_idx - smaller_window_size);
                        q_map.emplace(sp_cropped[max_left_idx] , std::make_pair(elem.first,max_idx - smaller_window_size));
                    }
                    catch (const std::runtime_error& e)
                    {
                        std::cerr << e.what() << std::endl;
                    }
                }


                if (elem.second - (max_idx+1+smaller_window_size) > smaller_window_size + 1)
                {
                    // the right median
                    vector<double> right_vec(sp_cropped.begin() + max_idx + 1  + smaller_window_size, sp_cropped.begin() + elem.second);
                    double median_right = MathOp::median(right_vec);
                    // normalise bins on the right by right median
                    for (int j = max_idx + smaller_window_size; j < elem.second; ++j) {
                        sp_cropped[j] /= median_right;
                    }

                    try
                    {
                        int max_right_idx = dsp.find_highest_peak(sp_cropped, max_idx + 1 + smaller_window_size, elem.second);
                        q_map.emplace(sp_cropped[max_right_idx] , std::make_pair(max_idx + 1 + smaller_window_size,elem.second));
                    }
                    catch (const std::runtime_error& e)
                    {
                        std::cerr << e.what() << std::endl;
                    }
                }


                all_max_ids.push_back(max_idx);
                std::cout << "Index of the maximum " << max_idx << " is added to all_max_ids." << std::endl;

                if (all_max_ids.size() >= breakpoints_limit)
                {
                    std::cout<<"More than " << breakpoints_limit << " breakpoints are detected."
                             << " Keeping only the top " << breakpoints_limit << " ones." << std::endl;
                    break;
                }


            }
            else // max_idx = -1, empty the q_map
            {
                q_map.clear();
            }

        }
    }

    if (all_max_ids.empty())
        throw std::logic_error("No breakpoints are detected. Perhaps change the configurations and try again.");

    std::cout<<"Sorting the all_max_ids..." <<std::endl;
    std::sort(all_max_ids.begin(), all_max_ids.end());
    std::cout<<"All max ids are sorted" <<std::endl;

    std::cout<<"Writing segmented regions to file..."<<std::endl;
    std::ofstream tree_file("./" + f_name_posfix + "_segmented_regions.txt");

    for (size_t k = 0; k < all_max_ids.size(); ++k) {
        // write it to file
        tree_file << all_max_ids[k] + window_size << endl;
    }
    std::cout<<"Segmented regions are written."<<std::endl;

    std::cout<<"Writing ordered segmented region sizes to file..."<<std::endl;
    std::ofstream reg_sizes_file("./" + f_name_posfix + "_segmented_region_sizes.txt");
    int cum_sum = 0;

    for (size_t k = 0; k < all_max_ids.size(); ++k) {
        // write it to file
        if (k == 0)
            reg_sizes_file << all_max_ids[k] + window_size << endl;
        else
            reg_sizes_file << all_max_ids[k] - all_max_ids[k-1] << endl;
    }
    reg_sizes_file << ub+window_size - all_max_ids[all_max_ids.size()-1] << endl; // add the last one

    std::cout<<"Segmented region sizes are written to file"<<std::endl;
    return EXIT_SUCCESS;
}