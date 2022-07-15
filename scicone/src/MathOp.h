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
#include <nlopt.hpp>
#include <cfloat>

#ifndef SC_DNA_NODE_H
#include "Node.h"
#endif //SC_DNA_NODE_H



using namespace std;

typedef struct {
    vector<double> z;
    double nu;
} segment_counts;

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
    static double log_replace_sum(const double &sum, const vector<double> &to_subtract, const vector<double> &to_add,
                                     const map<int, double> &unchanged_vals);
    static vector<double> combine_scores(vector<double> aic_vec);
    static vector<vector<double>> likelihood_ratio(vector<vector<double>> &mat, int window_size, vector<int> &known_breakpoints);
    static double breakpoint_log_prior(int k, int m, double mu);
    static double compute_linear_regression_slope(const vector<double>& x, const vector<double>& y);
    static double log_n_choose_k(unsigned long n, unsigned long k);
    static int random_uniform(int min, int max);
    template <class T>
    static double vec_max(vector<T> v);
    static int median_idx(int l, int r);
    template <class T>
    static double vec_avg(const vector<T> &v);
    template <class T>
    static T percentile_val(vector<T>, double percentile_val);
    template <class T>
    static double median(vector<T> v);
    template <class T>
    static double iqm(vector<T> v);
    template <class T>
    static double st_deviation(vector<T>& v);
    template<class T>
    static double third_quartile(vector<T> &v);
    template<class T>
    static double robust_mean(vector<T> v);
    template <class T>
    static double mat_moment(const vector<vector<T>> &v, int moment);
    template <class T>
    static double interquartile_range(vector<T> &v, bool top);
    static double compute_omega_insert_delete(Node *node, double lambda_r, double lambda_c, unsigned long K);
    static double compute_omega_condense_split(Node *node, double lambda_s, unsigned int n_regions);
    static double frobenius_avg(vector<vector<int>>& mat, vector<vector<int>>& ground_truth);
    static vector<long double> dirichlet_sample(size_t len, double alpha = 1.0);
    static vector<double> dirichlet_sample(vector<double> alphas);
    static double ll_linear_model(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
    static vector<double> compute_linear_regression_parameters(vector<double> &z, int window_size, double nu);
    static double huber_loss(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
    static double huber_mean(vector<double> &z, double delta);
};


#endif //SC_DNA_MATHOP_H
