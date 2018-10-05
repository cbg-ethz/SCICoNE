//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//


#include "MathOp.h"



template<class T>
double MathOp::vec_avg(vector<T> &v) {
    float average = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return average;
}

double MathOp::breakpoint_log_likelihood(std::vector<double> v)
{
    /*
     * Returns the max likelihood value for the Poisson distribution
     *
     * double ln_gamma = 0;
    for (auto const &i : v)
    {
        ln_gamma += log(tgamma(i));
    }*/
    // max likelihood:  std::log(max_ll_val) * sum(v) - (v.size() * max_ll_val)
    double max_ll_val = vec_avg(v);
    double term1,term2;
    // to avoid log(0) * 0
    if (vec_sum(v) == 0 && max_ll_val==0)
        term1 = 0.0;
    else
        term1 = std::log(max_ll_val) * vec_sum(v);

    term2 = (v.size() * max_ll_val);
    double ll =  term1 - term2;

    assert(!std::isnan(ll));
    return ll;
}

template<class T>
double MathOp::vec_sum(vector<T> &v) {
    return accumulate( v.begin(), v.end(), 0.0);
}

vector<vector<double>> MathOp::likelihood_ratio(vector<vector<double>>& mat, double window_size) {
    /*
     *
     * Computes the difference of the AIC_break and AIC_segment cases to tell whether to break or not
     * */

    //MathOp mo = MathOp();
    // the last breakpoint
    size_t bp_size = mat[0].size();

    //cout << bp_size << endl;
    size_t n_cells = mat.size();

    // u_int cell_no = 0;
    vector<vector<double>> aic_vec(bp_size,vector<double>(n_cells)); // fill constructor

    for (int i = 0; i < bp_size; ++i) {

        int start = i - window_size;
        int end = i + window_size;

        // end cannot be equal to bp_size because cellranger DNA is putting 0 to the end
        if (start < 0 || end >= bp_size) {
            //cout << "continue" << endl;
            continue;
        }

        for (int j = 0; j < n_cells; ++j) {

            vector<double> lbins = vector<double>(mat[j].begin() + start, mat[j].begin() + i);
            vector<double> rbins = vector<double>(mat[j].begin() + i, mat[j].begin() + end);
            vector<double> all_bins = vector<double>(mat[j].begin() + start, mat[j].begin() + end);
            /*cout << lbins.size() <<' ' << rbins.size() << ' ' << all_bins.size() << ' ' << i;
            cout << endl;
            cout << avg(lbins) << ' ' << avg(rbins) << ' ' << vec_avg(all_bins) <<endl;
            */
            // k is the degrees of freedom of the segment model
            u_int k_segment = 1;
            double ll_segment = breakpoint_log_likelihood(all_bins);
            double aic_segment = 2 * k_segment - 2 * ll_segment;

            // degrees of freedom is 2, lambda_r and lambda_l
            u_int k_break = 2;
            double ll_break = breakpoint_log_likelihood(lbins) + breakpoint_log_likelihood(rbins);
            double aic_break = 2 * k_break - 2 * ll_break;

            //cout << ll_break << ' ' << ll_segment << ' ' << aic_break << ' ' << aic_segment << endl;
            double aic_p = aic_segment - aic_break;
            //cout << aic_p << endl;

            aic_vec[i][j] = aic_p;
        }
    }
    return aic_vec;
}

long double MathOp::log_add(long double val1, long double val2)
{
    /*
     * Performs addition in the log space
     * uses the following identity:
     * log(a + b) = log(a * (1 + b/a)) = log a + log(1 + b/a)
     * make sure a > b
     * to avoid undefined values
     */

    // the input are in the log space
    long double a = val1;
    long double b = val2;

    if (b>a)
    {
        long double temp = a;
        a = b;
        b = temp;
    }


    // log1p(x) is used instead of log(1+x) to avoid precision loss
    long double res = a + log1p(exp(b-a));

//    if (std::isinf(res))
//        cout << "inf result detected";

    return res;
}



vector<double> MathOp::combine_scores(vector<double> aic_vec)
{

    /*
     * Combines cells.
     * Computes the combined evidence that the breakpoint occurred in any k of m cells.
     * Uses dynamic programming.
     *
     * */
    uint64_t m = aic_vec.size();
    vector<double> row1(m, 0.0);
    vector<double> row2;
    vector<double> res(1,0.0);

    // iterate over n_cells
    for (int j = 0; j < m; ++j) { // j : cells
        row2.clear();
        // inner computation of row2
        for (int k = 0; k < m-j; ++k) {
            double value=0.0;
            if (k==0)
            {
                value = row1[k] + aic_vec[k+j];
            }
            else
            {
                value = log_add(row2[k-1] , row1[k] + aic_vec[k+j]);

            }

            if (std::isinf(value))
                cerr << "inf value detected";

            row2.push_back(value);

        }
        //cout << " row1 size: " << row1.size() << " row2 size: " << row2.size() << " max j: " << aic_vec.size()-j<<endl;
        row1.clear();
        row1 = row2;
        double last = row2.back();
        last -= log_n_choose_k(m,j);

        res.push_back(last);

    }



    return res;
}

int MathOp::random_uniform(int min, int max) {
    int rand_val = 0;
    std::uniform_int_distribution<> dis(min,max);
    std::mt19937 &gen = SingletonRandomGenerator::get_generator();

    rand_val = dis(gen);
    return rand_val;
}

double MathOp::log_sum(const map<int, double> &map) {
/*
 * Performs normal sum in the log space
 * Takes a map as input
 * */

    // init max with the smallest value possible
    double max = numeric_limits<double>::lowest();

    for (auto const &i: map)
        if (i.second > max)
            max = i.second;

    double sum_in_normal_space = 0.0;

    for (auto const &i: map)
        sum_in_normal_space += exp(i.second-max);

    return log(sum_in_normal_space) + max;

}

double MathOp::log_sum(const vector<double> &vec) {
/*
 * Performs normal sum in the log space
 * Takes a vector as input
 * */

    // init max with the smallest value possible
    double max = std::numeric_limits<double>::lowest();

    for (auto const &i: vec)
        if (i > max)
            max = i;

    double sum_in_normal_space = 0.0;

    for (auto const &i: vec)
        sum_in_normal_space += exp(i-max);

    return log(sum_in_normal_space) + max;
}



double MathOp::log_replace_sum(const double &sum, const vector<double> &to_subtract, const vector<double> &to_add) {
/*
 * computes the sum in the log space, replaces the old elements of the vector (subtracts) then adds the new ones
 * */

    double max = sum;
    for (auto i: to_subtract)
        if (i > max)
            max = i;

    double sum_in_normal_space = 0.0;
    sum_in_normal_space += exp(sum-max); // += 1

    for (auto i: to_subtract)
    {
        double temp = exp(i-max);
        sum_in_normal_space -= temp;
    }

    if (sum_in_normal_space<= 1e-10)
        sum_in_normal_space = 1e-10;


    double subtracted_res = log(sum_in_normal_space) + max;

    vector<double> to_log_add = to_add;
    to_log_add.push_back(subtracted_res);

    double res = log_sum(to_log_add);
    assert(!std::isnan(res));

    return res;
}

double MathOp::breakpoint_log_prior(int k, int m, double mu) {
    /*
     * Returns the prior probability of having the breakpoint in k of the m cells.
     * 1-mu is the prob. of being a breakpoint
     * k is the n_cells that the breakpoint is present
     * m is the total n_cells
     * */

    assert((mu >= 0) && (mu <= 1));
    if (k == 0)
        return log(1-mu);
    else if (1 <= k && k <= m)
    {
        // compute the prior
        double res = 0.0;
        res += log(mu);
        res += 2*log_n_choose_k(m,k);
        res -= log(2*k-1);
        res -= log_n_choose_k(2*m,2*k);
        assert(!std::isnan(res));
        return res;
    }
    else
    {
        throw domain_error("k should be less than m and bigger than or equal to zero!");
    }

}

double MathOp::log_n_choose_k(int n, int k) {
    /*
     * Returns the log of n choose k using lgamma function
     * */
    assert(n >= k);
    return(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));

}

double MathOp::n_choose_k(int n, int k) {
    /*
     * Returns the binomial coefficient, n choose k using the gamma function
     * */
    assert (n >= k);
    double res = tgamma(n+1) / (tgamma(k+1) * tgamma(n-k+1));

    return res;
}

double MathOp::compute_omega(Node *node, double lambda_r, double lambda_c, double K) {
    /*
     * Computes the omega value used in assessing the delete node probabilities.
     * Note: omega is named w on the original publication.
     * */

    int r_i = static_cast<int>(node->c_change.size());

    if (r_i == 0)
        return 0.0;
    if (r_i > K)
        throw std::logic_error("There are not r distinct regions to sample");

    vector<int> c;
    for (auto const& x : node->c_change)
        c.push_back(abs(x.second));

    double omega_val = 1.0;
    omega_val *= pow(lambda_r, r_i - 1) * exp(-1 * lambda_r);
    double sum_cj_minus = 0.0;

    for (auto const &elem : c)
        sum_cj_minus += elem - 1;
    omega_val *= pow(lambda_c, sum_cj_minus);
    omega_val *= exp(-1 * r_i * lambda_c);

    omega_val /= pow(2, r_i);
    omega_val /= MathOp::n_choose_k(K, r_i);
    omega_val /= tgamma(r_i);

    double mul_cj_tgamma = 1.0;
    for (auto const &elem : c)
        mul_cj_tgamma *= tgamma(elem);

    omega_val /= mul_cj_tgamma;
    return omega_val;
}

