//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//


#include "MathOp.h"
#include "Utils.h"
#include "Lgamma.h"

template<class T>
double MathOp::vec_avg(const vector<T> &v) {

    double sum = accumulate( v.begin(), v.end(), 0.0);
    int v_size = v.size();
    double average = sum / v_size;
    return average;
}

double MathOp::breakpoint_log_likelihood(std::vector<double> v, double lambda, double nu)
{
    /*
     * Returns the max likelihood value for the Negative binomial distribution with mean: lambda and overdispersion: nu
     * Mean lambda is inferred by maximum likelihood approach.
     *
     */
    // max likelihood:  std::log(lambda) * sum(v) - (v.size() * lambda)

    double term1,term2;
    // to avoid log(0) * 0
    double v_sum = accumulate( v.begin(), v.end(), 0.0);
    if (v_sum == 0 && lambda==0)
        term1 = 0.0;
    else
        term1 = (log(lambda) - log(lambda+nu)) * v_sum;

    term2 = (v.size() * nu * log(lambda+nu));
    double ll =  term1 - term2;

    assert(!std::isnan(ll));
    return ll;
}



vector<vector<double>> MathOp::likelihood_ratio(vector<vector<double>> &mat, int window_size) {
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

    for (unsigned i = 0; i < bp_size; ++i) {

        int start = i - window_size;
        int end = i + window_size;

        // end cannot be equal to bp_size because cellranger DNA is putting 0 to the end
        if (start < 0 || end >= static_cast<int>(bp_size)) {
            //cout << "continue" << endl;
            continue;
        }

        for (unsigned j = 0; j < n_cells; ++j) {

            vector<double> lbins = vector<double>(mat[j].begin() + start, mat[j].begin() + i);
            vector<double> rbins = vector<double>(mat[j].begin() + i, mat[j].begin() + end);
            vector<double> all_bins = vector<double>(mat[j].begin() + start, mat[j].begin() + end);

            size_t n_bins = all_bins.size();

            double lambda_r = median(rbins); // max likelihood value is the avg. (poisson property)
            double lambda_l = median(lbins);
            double lambda_all = vec_avg(all_bins);

            vector<double> bin_positions(n_bins);
            for (size_t l = 0; l < n_bins; ++l) {
                bin_positions[l] = l+1; // 1 indexed
            }


            // predicting the lambdas_regression
            vector<double> lambdas_regression(n_bins);
            double beta_regression = compute_linear_regression_slope(bin_positions, all_bins);
            double mean_x = vec_avg(bin_positions);
            double mean_y = vec_avg(all_bins);

            for (size_t k = 0; k < n_bins; ++k)
            {
                lambdas_regression[k] = mean_y + beta_regression*(bin_positions[k] - mean_x); // prediction
                if (lambdas_regression[k] <= 0)
                    lambdas_regression[k] = 0.0001;
            }


            if (lambda_r == 0)
                lambda_r = 0.0001;
            if (lambda_l == 0)
                lambda_l = 0.0001;

            double aic_p = 0.0;
            // k is the degrees of freedom of the segment model
            u_int k_segment = 1;
            double ll_segment = breakpoint_log_likelihood(all_bins, lambda_all, 1.0);
            double aic_segment = 2 * k_segment - 2 * ll_segment;

            // degrees of freedom is 2, lambda_r and lambda_l
            u_int k_break = 2;
            double ll_break = breakpoint_log_likelihood(lbins, lambda_l, 1.0) +
                              breakpoint_log_likelihood(rbins, lambda_r, 1.0);
            double aic_break = 2 * k_break - 2 * ll_break;

            u_int k_slope = 2;
            double ll_slope = 0.0;

            for (size_t m = 0; m < n_bins; ++m) {
                ll_slope += breakpoint_log_likelihood(vector<double>(all_bins.begin()+m, all_bins.begin()+m+1), lambdas_regression[m], 1.0);
            }
            double aic_slope = 2 * k_slope - 2* ll_slope;


            aic_p = min(aic_segment, aic_slope) - aic_break;

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
    u_int m = aic_vec.size();
    vector<double> row1(m, 0.0);
    vector<double> row2;
    vector<double> res(1,0.0);

    // iterate over n_cells
    for (u_int j = 0; j < m; ++j) { // j : cells
        row2.clear();
        // inner computation of row2
        for (u_int k = 0; k < m-j; ++k) {
            double value=0.0;

            if (k==0)
                value = row1[k] + aic_vec[k+j];
            else
                value = static_cast<double>(log_add(row2[k-1] , row1[k] + aic_vec[k+j]));

            if (std::isinf(value))
                cerr << "inf value detected";

            row2.push_back(value);

        }
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
    boost::random::uniform_int_distribution<> dis(min,max);
    std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;

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



double MathOp::log_replace_sum(const double &sum, const vector<double> &to_subtract, const vector<double> &to_add,
                              const map<int, double> &unchanged_vals) {
    /*
     * computes the sum in the log space, replaces the old elements of the vector (subtracts) then adds the new ones
     * */

    double res = 0.0; // the return value

    if (to_add.size() >= unchanged_vals.size())
    {
        vector<double> to_log_add = to_add;
        for (auto const &u_map : unchanged_vals)
            to_log_add.push_back(u_map.second);
        res = log_sum(to_log_add);
    }
    else
    {

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

        if (sum_in_normal_space<= 1e-6)
        {
            //cout << "POSSIBLE PRECISION LOSS: " << sum_in_normal_space << endl;
            vector<double> to_log_add = to_add;
            for (auto const &u_map : unchanged_vals)
                to_log_add.push_back(u_map.second);
            res = log_sum(to_log_add);

        }
        else
        {
            double subtracted_res = log(sum_in_normal_space) + max;
            vector<double> to_log_add = to_add;
            to_log_add.push_back(subtracted_res);
            res = log_sum(to_log_add);
        }

    }

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
        res -= log_n_choose_k(2*2*m,k);
        assert(!std::isnan(res));
        return res;
    }
    else
    {
        throw domain_error("k should be less than m and bigger than or equal to zero!");
    }

}

double MathOp::log_n_choose_k(unsigned long n, unsigned long k) {
    /*
     * Returns the log of n choose k using lgamma function
     * */
    assert(n >= k);
    return(Lgamma::get_val(n+1) - Lgamma::get_val(k+1) - Lgamma::get_val(n-k+1));

}

double MathOp::compute_omega_insert_delete(Node *node, double lambda_r, double lambda_c, unsigned long K) {
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

    // rewrite in log space
    double omega_val = 0.0;
    omega_val += log(lambda_r) * (r_i - 1) + (-1 * lambda_r);
    double sum_cj_minus = 0.0;

    for (auto const &elem : c)
        sum_cj_minus += elem - 1;
    omega_val += log(lambda_c)* sum_cj_minus;
    omega_val += (-1 * r_i * lambda_c);

    omega_val -= log(2) * r_i;
    omega_val -= MathOp::log_n_choose_k(K, r_i);
    omega_val -= Lgamma::get_val(r_i);

    double mul_cj_lgamma = 0.0;
    for (auto const &elem : c)
        mul_cj_lgamma += Lgamma::get_val(elem);

    omega_val -= mul_cj_lgamma;
    return exp(omega_val);
}

double MathOp::frobenius_avg(vector<vector<int>> &mat, vector<vector<int>> &ground_truth) {
    /*
     * Computes and returns the Frobenius average of 2 matrices
     * */

    assert(mat.size() == ground_truth.size());
    assert(mat[0].size() == ground_truth[0].size());

    double delta = 0.0;

    size_t n = mat.size();
    size_t m = mat[0].size();
    for (unsigned i = 0; i < n; ++i)
    {
        for (unsigned j = 0; j < m; ++j)
        {
            delta += pow(mat[i][j] - ground_truth[i][j] , 2 );
        }
    }
    delta /= (n * m); //take the avg
    delta = sqrt(delta);

    return delta;
}

vector<long double> MathOp::dirichlet_sample(size_t len, double alpha) {
    /*
     * Samples a vector of length len that sum up to 1 using dirichlet with alpha parameter
     * */

    mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
    gamma_distribution<long double> distribution(alpha,1.0);

    vector<long double> y_vals(len);

    // sample from gamma
    long double sum_gamma = 0.0;
    for (unsigned i = 0; i < len; ++i) {
        y_vals[i] = distribution(gen);
        sum_gamma += y_vals[i];
    }

    // normalize them, that's it
    vector<long double> x_vals(len);
    for (unsigned j = 0; j < len; ++j) {
        x_vals[j] = y_vals[j]/sum_gamma;
    }


    // make sure probs sum up to 1
    long double sum = accumulate( x_vals.begin(), x_vals.end(), 0.0l); // l defines long double
    if (sum != 1.0) // add the 0.0000001 to the first element.
    {
        x_vals[0] += 1.0-sum;
    }
    return x_vals;
}

template<class T>
T MathOp::percentile_val(vector<T> vec, double percentile_val) {
    /*
     * Returns the value of the vector at the given percentile.
     * First sorts the vector.
     * The percentile value should be between 0 and 1.
     * */

    assert(percentile_val <=1.0 && percentile_val >= 0.0);

    // sort using the default operator<
    std::sort(vec.begin(), vec.end());

    int percentile_idx = static_cast<int>(percentile_val * vec.size())-1;

    return vec[percentile_idx];
}

double MathOp::compute_omega_condense_split(Node *node, double lambda_s, unsigned int n_regions) {

    /*
     * Computes the omega probability of the condense/split move.
     * */

    double omega_val = 1.0;

    vector<int> res_c_change = {};
    for (u_int i = 0; i < n_regions; ++i) {
        int c_val = 0; // value of the region in the child
        int p_val = 0; // value of the region in the parent

        if (Utils::is_empty_map(node->parent->c_change)) // EXCLUDE: parent is root
            return 0;

        if (node->c_change.count(i))
            c_val = node->c_change[i];
        if (node->parent->c_change.count(i))
            p_val = node->parent->c_change[i];

        if ((p_val == 0) && (c_val == 0)) // this region is not to be considered
            continue;

        if (p_val + c_val == 0) // EXCLUDE: any of the nodes cancel out
            return 0;

        int c = abs(c_val-p_val);
        double f_c = 1.0;

        if (c == 0) // c is zero
        {
            f_c *= exp(-1*lambda_s);
        }
        else if (c % 2 == 0) // c is even
        {
            f_c *= pow(lambda_s, c/2);
            f_c *= exp(-1*lambda_s);
            f_c /= 2*exp(Lgamma::get_val(c/2 +1)); //tgamma +1 = factorial
        }
        else // c is odd
        {
            f_c *= pow(lambda_s, (c-1)/2);
            f_c *= exp(-1*lambda_s);
            f_c /= 2*exp(Lgamma::get_val((c-1)/2 + 1));
        }
        omega_val *= f_c;
        res_c_change.push_back(c_val + p_val);
    }

    if(res_c_change.size() == 1 && res_c_change[0] == 1) // you cannot end up with 1 region and 1 event
        return 0;

    return omega_val;
}

template<class T>
double MathOp::median(vector<T> v) {

    /*
     * Returns the median value of the vector
     * */

    // sort using the default operator<
    std::sort(v.begin(), v.end());

    size_t n = v.size();
    if (n % 2 != 0) // odd
        return (double)v[n/2];
    return (double)(v[(n-1)/2] + v[n/2])/2.0; // even

}

template<class T>
double MathOp::st_deviation(vector<T> &v) {

    /*
     * Returns the standard deviation of a vector of elements.
     * */

    double mean = MathOp::vec_avg(v);
    double var = 0.0;

    for (unsigned i = 0; i < v.size(); ++i) {
        var += pow(v[i] - mean, 2);
    }
    var /= v.size();

    double stdev = sqrt(var);

    return stdev;
}

double MathOp::compute_linear_regression_slope(const vector<double> &x, const vector<double> &y) {
    /*
     * Computes the b1 (slope) of the formula Yi = b0 + b1Xi
     * */

    if(x.size() != y.size()){
        throw length_error("Length of y and x must be equal in simple linear regression");
    }

    double mean_x = vec_avg(x);
    double mean_y = vec_avg(y);

    double numerator = 0.0;
    double denominator = 0.0;

    for(unsigned i=0; i<x.size(); ++i){
        numerator += (x[i] - mean_x) * (y[i] - mean_y);
        denominator += (x[i] - mean_x) * (x[i] - mean_x);
    }

    if(denominator == 0){
        throw domain_error("Denominator cannot be zero");
    }

    return numerator / denominator;
}


template double MathOp::st_deviation(vector<double> &v);
template double MathOp::st_deviation(vector<int> &v);
template double MathOp::median(vector<double> v);
template double MathOp::vec_avg(const vector<double> &v);
template double MathOp::vec_avg(const vector<int> &v);
template double MathOp::percentile_val<double>(vector<double>, double percentile_val);
template long double MathOp::percentile_val<long double>(vector<long double>, double percentile_val);

