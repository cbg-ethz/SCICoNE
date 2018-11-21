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
#include <algorithm>
#include "SingletonRandomGenerator.h"

#ifndef SC_DNA_NODE_H
#include "Node.h"
#endif //SC_DNA_NODE_H



using namespace std;

class MathOp {

    /*
     * Provides mathematical operations
     * */

public:

    MathOp() = delete;
    virtual ~MathOp() = delete;
    static double breakpoint_log_likelihood(std::vector<double> v, double lambda, double nu);
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
    static double vec_avg(vector<T> &v);
    template <class T>
    static double vec_sum(vector<T> &);
    template <class T>
    static T percentile_val(vector<T>, double percentile_val);
    template <class T>
    static double median(vector<T> v);
    template <class T>
    static double st_deviation(vector<T>& v);
    static double compute_omega_insert_delete(Node *node, double lambda_r, double lambda_c, double K);
    static double compute_omega_condense_split(Node *node, double lambda_s, int n_regions, bool weighted);
    static double frobenius_avg(vector<vector<int>>& mat, vector<vector<int>>& ground_truth);
    static vector<long double> dirichlet_sample(int len, double alpha = 1.0);

};


#endif //SC_DNA_MATHOP_H
