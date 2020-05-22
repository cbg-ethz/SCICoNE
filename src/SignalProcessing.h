//
// Created by Tuncel  Mustafa Anil on 10/21/18.
//

#ifndef SC_DNA_SIGNALPROCESSING_H
#define SC_DNA_SIGNALPROCESSING_H

#include "MathOp.h"
#include <vector>
#include <fstream>

using namespace std;


class SignalProcessing {

public:

    template<class T>
    vector<T> crop(vector<T> &signal, int offset);
    vector<double> diff(vector<double> &signal);
    vector<double> sign(vector<double> &signal);
    void log_transform(vector<double> & signal);
    int evaluate_peak(vector<double> signal, vector<double> sp_cropped_copy, int lb, int ub, double threshold_coefficient,
                          int window_size);
    template<class T>
    int find_highest_peak(vector<T> &signal, int lb, int ub);
    vector<double> make_zero_mean(vector<double>& signal);
    vector<double> subtract_median(vector<double> &signal);
    void median_normalise(vector<double> &signal);
    void attenuate_values_below(vector<double> &signal, double threshold);
    vector<bool> filter_by_val(vector<double> &signal, double val);
    vector<int> create_region_sizes(vector<bool> peaks);
    vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size, int k_star, vector<int> &known_breakpoints, bool compute_lr=true, bool lr_only=false);
    vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size, int k_star, vector<int> &known_breakpoints, vector<vector<double>> &lr_vec, bool compute_lr=true, bool lr_only=false);
};


#endif //SC_DNA_SIGNALPROCESSING_H
