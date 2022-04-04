//
// Created by Tuncel  Mustafa Anil on 10/21/18.
//

#include "SignalProcessing.h"
#include "globals.cpp"

vector<double> SignalProcessing::diff(vector<double> &signal) {
    /*
     * Performs the diff filter on the input signal and returns the value.
     * The returned signal has length 1 less than the input signal.
     * Throws logic_error.
     * */

    if (signal.size() <= 1)
    {
        throw std::logic_error("The vector you perform diff filter must contain at least 2 elements.");
    }

    vector<double> res(signal.size()-1); // -1 because the first one cannot have a diff

    for (int i = 1; i < signal.size()-1; ++i) {
        res[i] = signal[i]-signal[i-1];
    }
    return res;
}

vector<double> SignalProcessing::sign(vector<double> &signal) {
    /*
     * Maps the input signal into a signal of {-1.0,1.0} values by their sign values.
     * Positive values will be mapped to 1.0 while negative will be mapped to -1.0
     * */

    vector<double> res(signal.size());

    for (int i = 0; i < signal.size(); ++i) {
        if (signal[i] > 0)
            res[i] = 1.0;
        else
            res[i] = -1.0;
    }
    return res;
}

vector<bool> SignalProcessing::filter_by_val(vector<double> &signal, double val) {
    /*
     * Creates and returns a vector of boolean containing true for elements in the signal equal to val, false otherwise.
     * */
    vector<bool> res(signal.size());

    for (int i = 0; i < signal.size(); ++i) {
        if (signal[i] == val)
            res[i] = true;
        else
            res[i] = false;
    }

    return res;
}


vector<int> SignalProcessing::create_region_sizes(vector<bool> peaks) {
    /*
     * Segments the bins until a peak is observed.
     * */

    vector<int> region_sizes;

    int size = 0;
    for (int i = 0; i < peaks.size(); ++i) {
        if (!peaks[i])
            size++;
        else
        {
            size++;
            region_sizes.push_back(size);
            size = 0;
        }
    }
    if (size != 0)
        region_sizes.push_back(size);

    return region_sizes;
}

vector<double> SignalProcessing::make_zero_mean(vector<double> &signal) {
    /*
     * Subtracts the mean value from signal and returns a zero mean signal.
     * */
    vector<double> res(signal.size());
    double avg = MathOp::vec_avg(signal);

    for (int i = 0; i < signal.size(); ++i) {
        res[i] = signal[i] - avg;
    }


    return res;
}

void SignalProcessing::attenuate_values_below(vector<double> &signal, double threshold) {
    /*
     * Sets the values below the threshold to zero in the input signal.
     * Modifies the input signal
     * */
    for (int i = 0; i < signal.size(); ++i) {
        if (signal[i] < threshold)
            signal[i] = 0.0;
    }

}

template<class T>
vector<T> SignalProcessing::crop(vector<T> &signal, int offset)
{
    /*
     * Crops the signal by the offset from both ends and returns the cropped signal
     * */

    assert(signal.size() > offset*2); // the signal should be bigger than the offset

    typename std::vector<T>::const_iterator first = signal.begin() + offset;
    typename std::vector<T>::const_iterator last = signal.end()  - offset;
    typename std::vector<T> cropped(first, last);

    return cropped;
}

vector<double> SignalProcessing::subtract_median(vector<double> &signal) {
    /*
     * Substracts the median from the signal and returns the new signal.
     * */

    vector<double> res(signal.size());
    double median_val = MathOp::percentile_val(signal, 0.5);

    for (int i = 0; i < signal.size(); ++i) {
        res[i] = signal[i] - median_val;
    }

    return res;
}

void SignalProcessing::median_normalise(vector<double> &signal) {
    /*
     * Divides each element in the vector by its median value.
     * */

    double med = MathOp::median(signal);

    for (int i = 0; i < signal.size(); ++i) {
        signal[i] /= med;
    }
}

int SignalProcessing::evaluate_peak(vector<double> signal, vector<double> sp_cropped_copy, int lb, int ub, double threshold_coefficient,
                                    int window_size) {
    /*
     * Returns the index of the highest peak above the threshold in a signal between the lb (lower bound) ub (upper bound) intervals.
     * Copy signal contains NaN values and it is passed by value.
     * */

    int max_idx;
    try
    {
        max_idx = this->find_highest_peak(signal, lb, ub);
    }catch (const std::runtime_error& e)
    {
        std::cerr << " a runtime error was caught during the find_highest_peak function, with message '"
                  << e.what() << "'\n";
        return -1; // peak is not selected as a breakpoint
    }

    double max_val = signal[max_idx];

    // remove the values at NaN index from the stdev computation
    size_t n_removed = 0;
    for (int k = 0; k < sp_cropped_copy.size(); ++k) {

        if (std::isnan(sp_cropped_copy[k]))
        {
            signal.erase(signal.begin() + k - n_removed);
            n_removed += 1;
        }
    }
    // take log of the signal
    this->log_transform(signal);

    // double stdev = MathOp::st_deviation(signal);
    // double third_q = MathOp::third_quartile(signal);
    // double median = MathOp::median(signal);
    // double range = third_q - median;
    double range = MathOp::interquartile_range(signal, true); // thirdq - median
    // double range = stdev;
    double threshold = threshold_coefficient * range;

    // use log of max_val
    max_val = log(max_val);

    if (threshold == 0) // reject the breakpoint if stdev is zero
        return -1;

    if (verbosity > 0)
    {
        std::ofstream bp_vals_file("./" + f_name_posfix + "_all_bps_comparison.csv", std::ios_base::app);
        bp_vals_file << max_idx + window_size << ',' << max_val << ',' << range << std::endl; // 10 for window size
    }

    if (max_val > threshold)
        return max_idx;
    else
        return -1;

}

void SignalProcessing::log_transform(vector<double> &signal) {
    /*
     * Takes the input signal and performs log transform on every element of it.
     * Mutates the original matrix
     * */

    for (int i = 0; i < signal.size(); ++i) {
        signal[i] = log(signal[i]);
    }
}

template<class T>
int SignalProcessing::find_highest_peak(vector<T> &signal, int lb, int ub) {
    /*
     * Returns the index of the highest peak above the threshold in a signal between the lb (lower bound) ub (upper bound) intervals.
     * Throws runtime error.
     * */

    assert(lb < ub);
    assert(ub <= signal.size());

    // get the subvector
    vector<double>::const_iterator first = signal.begin() + lb;
    vector<double>::const_iterator last = signal.begin() + ub;
    vector<double> sub(first, last);

    vector<double> peaks = this->diff(sub);
    peaks = this->sign(peaks);
    peaks = this->diff(peaks); // real peaks are peaks -1 because of double differentiation.
    // after 1st diff. peak is the last positive. after 2nd diff. peak is zero but peak+1 is -2

    vector<bool> breakpoints = this->filter_by_val(peaks, -2.0);

    vector<int> bp_indices;

    for (int i = 0; i < breakpoints.size(); ++i) {
        if (breakpoints[i])
            bp_indices.push_back(i-1); // push i-1 because that's the real peak
    }


    if (bp_indices.empty()) // no breakpoints are found within this range
        throw std::runtime_error("no breakpoints are found within the interval: (" + to_string(lb) + ',' + to_string(ub) + ')');

    double max_val = numeric_limits<double>::lowest();
    int max_idx = -1;

    for (int j = 0; j < bp_indices.size(); ++j) {
        if (sub[bp_indices[j]] > max_val)
        {
            max_val = sub[bp_indices[j]];
            max_idx = bp_indices[j];
        }
    }

    return max_idx + lb;
}

vector<double>
SignalProcessing::breakpoint_detection(vector<vector<double>> &mat, int window_size, int k_star, vector<int> &known_breakpoints, bool smoothed, bool compute_lr, bool lr_only) {
  vector<vector<double>> lr_vec;
  return SignalProcessing::breakpoint_detection(mat, window_size, k_star, known_breakpoints, lr_vec, smoothed, compute_lr, lr_only);
}

vector<double>
SignalProcessing::breakpoint_detection(vector<vector<double>> &mat, int window_size, int k_star, vector<int> &known_breakpoints, vector<vector<double>> &lr_vec, bool smoothed, bool compute_lr, bool lr_only) {
    /*
     * Performs the breakpoint detection
     * window_size: there cannot be multiple breakpoints within a window_size
     * k_star: min number of cells to consider for a bin being breakpoint
     * ul: max number of cells to consider for a bin being breakpoint
     * returns the evidence of each bin being a breakpoint
     * */
    // TODO: get rid of unnecessary push_backs

    size_t n_cells = mat.size();

    // compute the LR scores
    if (compute_lr)
      lr_vec = MathOp::likelihood_ratio(mat, window_size, known_breakpoints, smoothed);
    else
      std::cout << "Skipping LR computation" << std::endl;

    if (verbosity > 0)
    {
        std::ofstream lr_vec_file("./" + f_name_posfix + "_lr_vec" + ".csv");
        for (auto const &v1: lr_vec) {
            for (size_t i = 0; i < v1.size(); i++)
            {
                if (i == v1.size()-1) // the last element
                    lr_vec_file << v1[i];
                else // add comma
                    lr_vec_file << v1[i] << ',';
            }
            lr_vec_file << endl;
        }
    }

    if (lr_only) {
      vector<double> sp_vec_dummy;
      return sp_vec_dummy;
    }

    size_t n_breakpoints = lr_vec.size();

    std::cout << "Combining scores..." << std::endl;
    vector<vector<double>> sigma(n_breakpoints,vector<double>(n_cells+1)); // +1 because combine scores considers
    // the breakpoint occurring in zero cells as well

    for (size_t i = 0; i < n_breakpoints; ++i) // compute sigma matrix
        sigma[i] = MathOp::combine_scores(lr_vec[i]);

    vector<double> log_priors;
    log_priors.reserve(n_cells+1);
    for (size_t j = 0; j < n_cells+1; ++j)
        log_priors.push_back(MathOp::breakpoint_log_prior(j, n_cells,0.001));

    vector<vector<double>> log_posterior(n_breakpoints,vector<double>(n_cells+1));
    for (size_t k = 0; k < n_breakpoints; ++k) {
        for (size_t j = 0; j < n_cells+1; ++j) {
            double val = log_priors[j] + sigma[k][j];
            log_posterior[k][j] = val;
        }
    }

    vector<vector<double>> posterior(n_breakpoints,vector<double>(n_cells+1));
    vector<vector<double>> posterior_k(n_breakpoints,vector<double>(n_cells+1));

    vector<double> s_p;
    vector<double> expected_k_vector;

    for (size_t l = 0; l < n_breakpoints; ++l)
    {
        double max_all = *max_element(log_posterior[l].begin(), log_posterior[l].end());

        for (int j = 0; j < log_posterior[l].size(); ++j) {
            log_posterior[l][j] -= max_all;
            posterior[l][j] = exp(log_posterior[l][j]);
            posterior_k[l][j] = j * posterior[l][j];
        }

        double expected_nom = std::accumulate(posterior_k[l].begin(), posterior_k[l].end(), 0.0);
        if (expected_nom != 0.0)
            expected_nom = log(expected_nom);

        double sp_denom = std::accumulate(posterior[l].begin(), posterior[l].end(), 0.0);
        if (sp_denom != 0.0)
            sp_denom = log(sp_denom);

        double log_expected_cells = expected_nom - sp_denom;
        expected_k_vector.push_back(log_expected_cells);

        double max_local = *max_element(log_posterior[l].begin(), log_posterior[l].begin() + k_star - 1);
//        double max_local_ub = *max_element(log_posterior[l].begin() + ul, log_posterior[l].end());

//        double max_local = std::max(max_local, max_local_ub);
        for (int j = 0; j < log_posterior[l].size(); ++j) {
            posterior[l][j] = exp(log_posterior[l][j] - max_local);
            if (posterior[l][j] == 0)
              posterior[l][j] = 1e-8;
        }

        double sp_num_total = std::accumulate(posterior[l].begin(), posterior[l].begin() + k_star - 1, 0.0);
//        double sp_num_ub =  std::accumulate(posterior[l].begin() + ul, posterior[l].end(), 0.0);

        if (sp_num_total != 0.0)
            sp_num_total = log(sp_num_total) + max_local;

        double sp_val = sp_denom - sp_num_total;

        s_p.push_back(sp_val);

    }

    if (verbosity > 0)
    {
        std::ofstream log_posterior_file("./" + f_name_posfix + "_log_posterior_vec.csv");
        std::ofstream expected_k_vec_file("./" + f_name_posfix + "_expected_k_vec.csv");
        std::ofstream sp_file("./" + f_name_posfix + "_sp_vec.csv");

        for (auto const &v1: log_posterior) {
            for (size_t i = 0; i < v1.size(); i++)
            {
                if (i == v1.size()-1) // the last element
                    log_posterior_file << v1[i];
                else // add comma
                    log_posterior_file << v1[i] << ',';
            }
            log_posterior_file << endl;
        }

        for (int i = 0; i < expected_k_vector.size(); ++i) {
            if (i == expected_k_vector.size()-1) // the last element
                expected_k_vec_file << expected_k_vector[i];
            else // add comma
                expected_k_vec_file << expected_k_vector[i] << ',';
        }
        expected_k_vec_file << std::endl;

        for (int i = 0; i < s_p.size(); ++i) {
            if (i == s_p.size()-1) // the last element
                sp_file << s_p[i];
            else // add comma
                sp_file << s_p[i] << ',';
        }
        sp_file << std::endl;

    }
    std::cout << "Done." << std::endl;

    return s_p;

}

template int SignalProcessing::find_highest_peak(vector<double> &signal, int lb, int ub);
template vector<double> SignalProcessing::crop<double>(vector<double>& signal, int offset);
template vector<long double> SignalProcessing::crop<long double>(vector<long double>& signal, int offset);
