//
// Created by Tuncel  Mustafa Anil on 10/21/18.
//

#ifndef SC_DNA_SIGNALPROCESSING_H
#define SC_DNA_SIGNALPROCESSING_H

#include "MathOp.h"
#include <vector>
using namespace std;


class SignalProcessing {

public:


    template<class T>
    vector<T> crop(vector<T> &signal, int offset);
    vector<double> diff(vector<double> &signal);
    vector<double> sign(vector<double> &signal);
    vector<double> make_zero_mean(vector<double>& signal);
    void attenuate_values_below(vector<double> &signal, double threshold);
    vector<bool> filter_by_val(vector<double> &signal, double val);

    vector<int> create_region_sizes(vector<bool> peaks);



};


#endif //SC_DNA_SIGNALPROCESSING_H
