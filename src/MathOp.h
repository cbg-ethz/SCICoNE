//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//

#ifndef SC_DNA_MATHOP_H
#define SC_DNA_MATHOP_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <assert.h>
#include <random>


using namespace std;

class MathOp {

public:
    MathOp() = delete;
    virtual ~MathOp() = delete;

    static double log_likelihood(vector<double>);
    static long double log_add(long double val1, long double val2);
    static vector<double> combine_scores(vector<double> aic_vec);
    static vector<vector<double>> likelihood_ratio(vector<vector<double>> mat, double window_size);

    static int random_uniform(int min, int max);

    template <class T>
    static double vec_avg(vector<T>);

    template <class T>
    static double vec_sum(vector<T>);



};


#endif //SC_DNA_MATHOP_H
