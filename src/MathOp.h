//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//

#ifndef SC_DNA_MATHOP_H
#define SC_DNA_MATHOP_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cassert>
#include <random>
#include <map>
#include "SingletonRandomGenerator.h"

#ifndef SC_DNA_NODE_H
#include "Node.h"
#endif //SC_DNA_NODE_H





using namespace std;

class MathOp {

public:
    MathOp() = delete;
    virtual ~MathOp() = delete;

    static double breakpoint_log_likelihood(vector<double> v, double nu=1.0);
    static long double log_add(long double val1, long double val2);
    static double log_sum(const map<int, double> &map); // map version
    static double log_sum(const vector<double> &vec); // vector version

    static double log_replace_sum(const double& sum, const vector<double>& to_subtract, const vector<double>& to_add);

    static vector<double> combine_scores(vector<double> aic_vec);
    static vector<vector<double>> likelihood_ratio(vector<vector<double>>& mat, double window_size);

    static double breakpoint_log_prior(int k, int m, double mu);
    static double log_n_choose_k(int n, int k);
    static double n_choose_k(int n, int k);

    static int random_uniform(int min, int max);

    template <class T>
    static double vec_avg(vector<T> &);

    template <class T>
    static double vec_sum(vector<T> &);

    static double compute_omega(Node *node, double lambda_r, double lambda_c, double K);


    static double frobenius_avg(vector<vector<int>>& mat, vector<vector<int>>& ground_truth);

    static vector<long double> dirichlet_sample(int len, double alpha = 1.0);


};


#endif //SC_DNA_MATHOP_H
