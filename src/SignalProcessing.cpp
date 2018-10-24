//
// Created by Tuncel  Mustafa Anil on 10/21/18.
//

#include "SignalProcessing.h"


vector<double> SignalProcessing::diff(vector<double> &signal) {
    /*
     * Performs the diff filter on the input signal and returns the value.
     * The returned signal has length 1 less than the input signal.
     * */

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

template vector<double> SignalProcessing::crop<double>(vector<double>& signal, int offset);
template vector<long double> SignalProcessing::crop<long double>(vector<long double>& signal, int offset);


