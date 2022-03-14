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

template<class T>
double MathOp::mat_moment(const vector<vector<T>> &v, int moment) {

    int rows = v.size();
    int cols = v[0].size();

    double sum = 0.;
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        sum = sum + pow(v[i][j], moment);
      }
    }
    double average = sum / (rows*cols);

    return average;
}

double MathOp::breakpoint_log_likelihood(std::vector<double> v, double lambda, double nu, bool normal)
{
  double ll = 0;
  if (normal)
    ll = normal_log_likelihood(v, lambda, nu);
  else
    ll = nb_log_likelihood(v, lambda, nu);
  return ll;
}

double MathOp::nb_log_likelihood(std::vector<double> v, double lambda, double nu)
{
    /*
     * Returns the log likelihood for the Negative binomial distribution with mean: lambda and overdispersion: nu
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

double MathOp::normal_log_likelihood(std::vector<double> v, double mu, double sigma)
{
    /*
     * Returns the log likelihood for the Normal distribution with mean: mu and stddev: sigma
     */

    double var = pow(sigma,2);
    double term1 = v.size() * 0.5 * log(2*M_PI*var);
    double term2 = 0;
    for (size_t i = 0; i < v.size(); ++i)
      term2 += pow(v[i] - mu, 2);
    term2 *= -1/(2*var);

    double ll =  -term1 + term2;
    assert(!std::isnan(ll));
    return ll;
}


vector<vector<double>> MathOp::likelihood_ratio(vector<vector<double>> &mat, int window_size, vector<int> &known_breakpoints, bool normal) {
    /*
     *
     * Computes the difference of the likelihood_break and likelihood_segment cases to tell whether to break or not
     * */

    // Estimate nu with method of moments, assuming all cells are in the same state
    double global_mean = mat_moment(mat, 1);
    double global_moment_2 = mat_moment(mat, 2);
    double nu = 0.;
    if (not normal)
      nu = pow(global_mean, 2) / global_moment_2;
    else
      nu = global_moment_2 - pow(global_mean, 2);
    std::cout << "Method of moments estimated nu=" << nu << std::endl;

    //MathOp mo = MathOp();
    // the last breakpoint
    size_t n_bins = mat[0].size();

    size_t n_regions = known_breakpoints.size() - 1; // must include first and last bin
    std::cout << "n_regions: " << n_regions << std::endl;
    size_t n_cells = mat.size();

    // u_int cell_no = 0;
    vector<vector<double>> lr_vec(n_bins, vector<double>(n_cells)); // LR of each bin for each cell

    // Parallelize cells and known regions
    #pragma omp parallel for
    for (size_t j = 0; j < n_cells; ++j) {
        for (size_t r = 0; r < n_regions; ++r) {
            for (size_t i = known_breakpoints[r]; i < known_breakpoints[r+1]; ++i) {
                int start = i - window_size;
                int end = i + window_size;

                // end cannot be equal to n_bins because cellranger DNA is puts 0 to the end
                if (start < 0 || end >= static_cast<int>(n_bins)) {
                    continue;
                }

                vector<double> lbins = vector<double>(mat[j].begin() + start, mat[j].begin() + i);
                vector<double> rbins = vector<double>(mat[j].begin() + i, mat[j].begin() + end);
                vector<double> all_bins = vector<double>(mat[j].begin() + start, mat[j].begin() + end);

                size_t n_bins = all_bins.size();

                vector<double> bin_positions(n_bins);
                for (size_t l = 0; l < n_bins; ++l) {
                    bin_positions[l] = l+1; // 1 indexed
                }

                // 1. Likelihood of null model, where mean of all bins in segment follows a linear model.
                //   The mean across the bins in the window can change according to a linear model,
                //   which includes the case where all the bins have the same min (if the slope is zero)
                vector<double> lambdas_segment(n_bins);
                vector<double> regression_parameters = compute_linear_regression_parameters(all_bins, window_size, nu, normal);
                double alpha = regression_parameters[0];
                double beta = regression_parameters[1];

                for (size_t k = 0; k < n_bins; ++k)
                {
                    lambdas_segment[k] = alpha + beta*bin_positions[k]; // prediction
                    if (lambdas_segment[k] <= 0)
                        lambdas_segment[k] = 0.0001;
                }

                double ll_segment = 0;
                for (size_t m = 0; m < n_bins; ++m) {
                    ll_segment += breakpoint_log_likelihood(vector<double>(all_bins.begin()+m, all_bins.begin()+m+1), lambdas_segment[m], nu, normal);
                }

                // 2. Likelihood of breakpoint model, where left and right bin segments have different means
                //    This means that if there is a breakpoint there is a step change between the two semi segments
                double lambda_r = robust_mean(rbins);
                double lambda_l = robust_mean(lbins);
                double lambda_all = vec_avg(all_bins);

                // make sure lambda_all is between the left and right bounds
                lambda_all = std::max(lambda_all, std::min(lambda_r, lambda_l));
                lambda_all = std::min(lambda_all, std::max(lambda_r, lambda_l));

                 // The distance between the left and right segments must be > lambda_all/4, so we update the
                 // segments accordingly
                if (lambda_r > lambda_l) {
                  double gap = lambda_r - lambda_l;
                  double gap_thres = lambda_all/4.0;
                  double scaling = std::max(lambda_all/(4.0*gap), 1.0);
                  if (gap < gap_thres) {
                    lambda_r = lambda_all + scaling*(lambda_r-lambda_all);
                    lambda_l = lambda_all - scaling*(lambda_all-lambda_l);
                  }
                } else if (lambda_r < lambda_l) {
                  double gap = lambda_l - lambda_r;
                  double gap_thres = lambda_all/4.0;
                  double scaling = std::max(lambda_all/(4.0*gap), 1.0);
                  if (gap < gap_thres) {
                    lambda_l = lambda_all + scaling*(lambda_l-lambda_all);
                    lambda_r = lambda_all - scaling*(lambda_all-lambda_r);
                  }
                }
                if (lambda_r == 0)
                    lambda_r = 0.0001;
                if (lambda_l == 0)
                    lambda_l = 0.0001;

                double ll_break = breakpoint_log_likelihood(lbins, lambda_l, nu, normal) +
                                  breakpoint_log_likelihood(rbins, lambda_r, nu, normal);

                // 3. Output the difference between the two models' likelihoods
                lr_vec[i][j] = 2*(ll_break - ll_segment);
            }
        }
    }
    return lr_vec;
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
        res -= log_n_choose_k(2*m,2*k);
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

vector<double> MathOp::dirichlet_sample(vector<double> alphas) {
    /*
     * Samples a vector of length len that sum up to 1 using dirichlet with various alpha parameters
     * */

    mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
    vector<double> y_vals(alphas.size());
    vector<double> x_vals(alphas.size());
    // sample from gamma
    double sum_gamma = 0.0;

    for (int i = 0; i < alphas.size(); ++i) {
        gamma_distribution<double> distribution(alphas[i],1.0);

        y_vals[i] = distribution(gen);
        sum_gamma += y_vals[i];
    }

    // normalize them, that's it
    for (unsigned j = 0; j < alphas.size(); ++j) {
        x_vals[j] = y_vals[j]/sum_gamma;
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
    if (percentile_idx < 0) {
      percentile_idx = 0;
    }

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
double MathOp::vec_max(vector<T> v) {
  double currentMax = -DBL_MAX;

  for (u_int i = 0; i < v.size(); ++i){
    if (v[i] > currentMax)
      currentMax = v[i];
  }

  return currentMax;
}

int MathOp::median_idx(int l, int r) {
    int s = r - l + 1;
    s = (s + 1) / 2 - 1;
    return s + l;
}

template<class T>
double MathOp::iqm(vector<T> v) {
    std::sort(v.begin(), v.end());
    int n = v.size();

    // Index of median of entire data
    int mid_index = median_idx(0, n);

    // Index of median of first half
    int q1_index = median_idx(0, mid_index);

    // Index of median of second half
    int q3_index = median_idx(mid_index + 1, n);

    // Mean within Q1 and Q3
    double sum = 0;
    int counter = 0;
    for (int i = q1_index; i <= q3_index; i++) {
      sum = sum + v[i];
      counter = counter + 1;
    }

    double res = sum / counter;

    return res;
}

template<class T>
double MathOp::robust_mean(vector<T> v) {
    std::sort(v.begin(), v.end());
    int n = v.size();

    // Get data below 5th percentile
    double lower_perc = MathOp::percentile_val(v, 0.05);

    // Get data above 95th percentile
    double upper_perc = MathOp::percentile_val(v, 0.95);

    // Replace with values below and above limits
    vector<double> new_v(n);
    for (int i = 0; i < n; ++i) {
      if (v[i] <= lower_perc)
        new_v[i] = lower_perc;
      else if (v[i] >= upper_perc)
        new_v[i] = upper_perc;
      else
        new_v[i] = v[i];
    }

    // Take mean of new vector
    double res = vec_avg(new_v);

    return res;
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

template<class T>
double MathOp::third_quartile(vector<T> &v) {

    /*
     * Returns the third quartile of a vector of elements.
     * */
     int size = v.size();

     sort(v.begin(), v.end());

     int mid = size/2;
     double median;
     median = size % 2 == 0 ? (v[mid] + v[mid-1])/2 : v[mid];

     vector<double> third;

     for (int i = mid; i!= size; ++i) {
         third.push_back(v[i]);
     }

     double trd;

     int side_length = 0;

     if (size % 2 == 0)
     {
         side_length = size/2;
     }
     else {
         side_length = (size-1)/2;
     }

     trd = (size/2) % 2 == 0 ? (third[side_length/2]/2 + third[(side_length-1)/2])/2 : third[side_length/2];

    return trd;
}

template<class T>
double MathOp::interquartile_range(vector<T> &v, bool top) {
    /*
     * Returns the IQR of a vector
     * */
     int size = v.size();

     auto const Q1 = v.size() / 4;
     auto const Q2 = v.size() / 2;
     auto const Q3 = Q1 + Q2;

     std::nth_element(v.begin(),          v.begin() + Q1, v.end());
     std::nth_element(v.begin() + Q1 + 1, v.begin() + Q2, v.end());
     std::nth_element(v.begin() + Q2 + 1, v.begin() + Q3, v.end());

     int mid = size/2;
     double median;
     median = size % 2 == 0 ? (v[mid] + v[mid-1])/2 : v[mid];

     vector<double> first;
     vector<double> third;

     for (int i = 0; i!=mid; ++i)
        first.push_back(v[i]);

      for (int i = mid; i!= size; ++i)
        third.push_back(v[i]);

     double fst;
     double trd;

     int side_length = 0;

     if (size % 2 == 0)
     {
         side_length = size/2;
     }
     else {
         side_length = (size-1)/2;
     }

     fst = (size/2) % 2 == 0 ? (first[side_length/2]/2 + first[(side_length-1)/2])/2 : first[side_length/2];
     trd = (size/2) % 2 == 0 ? (third[side_length/2]/2 + third[(side_length-1)/2])/2 : third[side_length/2];

     double iqr = trd - fst;

     if (top)
      iqr = trd - median;

    return iqr;
}

double MathOp::huber_loss(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
  segment_counts *d = reinterpret_cast<segment_counts*>(my_func_data);
  vector<double> z = d->z;
  double delta = d->nu;

  int size = z.size();
  double res = 0;
  int a = 0;

  for (size_t l = 0; l < size; ++l) {
    a = z[l] - x[0];
    res = res + pow(delta, 2) * (pow(1 + pow(a/delta, 2), 0.5) - 1);
  }

  return res;
}

double MathOp::huber_mean(vector<double> &z, double nu) {
    /*
     * Computes the robust mean estimate of the bins in z
     * */

    segment_counts func_data = {
      z = z,
      nu = nu,
    };

    nlopt::opt opt(nlopt::LN_BOBYQA, 1);
    std::vector<double> lb(1);
    lb[0] = 0;
    opt.set_lower_bounds(lb);
    std::vector<double> ub(1);
    ub[0] = HUGE_VAL;
    opt.set_upper_bounds(ub);
    opt.set_min_objective(huber_loss, &func_data);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(1);
    x[0] = vec_avg(z);
    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return x[0]; // optimized alpha, beta
}

double MathOp::ll_linear_model_nb(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    segment_counts *d = reinterpret_cast<segment_counts*>(my_func_data);
    vector<double> z = d->z;
    double nu = d->nu; // inverse overdispersion parameter: lower nu gets you more dispersion. Assume low disperson

    int size = z.size();
    double lambda = 0;
    double res = 0;

    for (size_t l = 0; l < size; ++l) {
      lambda = x[0] + (l+1)*x[1];
      // lambda = x[0];
      res = res + z[l]*(log(lambda) - log(lambda + nu)) - nu*log(lambda + nu);
    }

    return res;
}

double MathOp::ll_linear_model_normal(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    segment_counts *d = reinterpret_cast<segment_counts*>(my_func_data);
    vector<double> z = d->z;
    double nu = d->nu;

    double var = pow(nu,2);
    int size = z.size();
    double lambda = 0;
    double res = 0;

    for (size_t l = 0; l < size; ++l) {
      lambda = x[0] + (l+1)*x[1];
      // lambda = x[0];
      res = res - pow(z[l] - lambda, 2) / (2*var);
    }

    return res;
}

vector<double> MathOp::compute_linear_regression_parameters(vector<double> &z, int window_size, double nu, bool normal) {
    /*
     * Computes the alpha and beta of y = alpha + beta * x for z = NB(mean=y, idisp=nu)
     * */

    double lambda_all = vec_avg(z);

    segment_counts func_data = {
      z = z,
      nu = nu,
    };

    nlopt::opt opt(nlopt::LN_BOBYQA, 2);
    std::vector<double> lb(2);
    lb[0] = 0; lb[1] = - 1.0/8.0 * 1.0/window_size; // lower bounds on alpha, beta
    opt.set_lower_bounds(lb);
    std::vector<double> ub(2);
    ub[0] = HUGE_VAL; ub[1] = 1.0/8.0 * 1.0/window_size; // upper bounds on alpha, beta
    opt.set_upper_bounds(ub);
    if (normal)
      opt.set_max_objective(ll_linear_model_normal, &func_data);
    else
      opt.set_max_objective(ll_linear_model_nb, &func_data);
    opt.set_xtol_rel(1e-4);

    std::vector<double> x(2);
    x[0] = vec_avg(z); x[1] = 0.; // initial alpha, beta
    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return x; // optimized alpha, beta
}




template double MathOp::vec_max(vector<double> v);
template double MathOp::st_deviation(vector<double> &v);
template double MathOp::st_deviation(vector<int> &v);
template double MathOp::median(vector<double> v);
template double MathOp::third_quartile(vector<double> &v);
template double MathOp::iqm(vector<double> v);
template double MathOp::vec_avg(const vector<double> &v);
template double MathOp::vec_avg(const vector<int> &v);
template double MathOp::mat_moment(const vector<vector<int>> &v, int moment);
template double MathOp::mat_moment(const vector<vector<double>> &v, int moment);
template double MathOp::percentile_val<double>(vector<double>, double percentile_val);
template long double MathOp::percentile_val<long double>(vector<long double>, double percentile_val);
template double MathOp::interquartile_range(vector<double> &v, bool top);
template double MathOp::robust_mean(vector<double> v);
