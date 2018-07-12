//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//

#include "MathOp.h"

MathOp::MathOp() {}

MathOp::~MathOp() {

}

template<class T>
double MathOp::vec_avg(std::vector<T> v) {
    float average = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return average;
}

double MathOp::log_likelihood(std::vector<double> v)
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

    assert(!isnan(ll));
    return ll;
}

template<class T>
double MathOp::vec_sum(std::vector<T> v) {
    return accumulate( v.begin(), v.end(), 0.0);
}

vector<vector<double>> MathOp::likelihood_ratio(vector<vector<double>> mat, double window_size) {
    /*
     *
     * Computes the difference of the AIC_break and AIC_segment cases to tell whether to break or not
     * */

    MathOp mo = MathOp();
    // the last breakpoint
    u_int bp_size = mat[0].size();

    //cout << bp_size << endl;
    u_int n_cells = mat.size();

    // u_int cell_no = 0;
    vector<vector<double>> aic_vec;
    for (int i = 0; i < bp_size; ++i) {

        int start = i - window_size;
        int end = i + window_size;

        // end cannot be equal to bp_size because cellranger DNA is putting 0 to the end
        if (start < 0 || end >= bp_size) {
            //cout << "continue" << endl;
            continue;
        }

        vector<double> aic_cell = vector<double>();

        for (int j = 0; j < n_cells; ++j) {


            vector<double> vect = mat[j];

            vector<double> lbins = vector<double>(vect.begin() + start, vect.begin() + i);
            vector<double> rbins = vector<double>(vect.begin() + i, vect.begin() + end);
            vector<double> all_bins = vector<double>(vect.begin() + start, vect.begin() + end);
            /*cout << lbins.size() <<' ' << rbins.size() << ' ' << all_bins.size() << ' ' << i;
            cout << endl;
            cout << avg(lbins) << ' ' << avg(rbins) << ' ' << vec_avg(all_bins) <<endl;
            */
            // k is the degrees of freedom of the segment model
            u_int k_segment = 1;
            double ll_segment = mo.log_likelihood(all_bins);
            double aic_segment = 2 * k_segment - 2 * ll_segment;

            // degrees of freedom is 2, lambda_r and lambda_l
            u_int k_break = 2;
            double ll_break = mo.log_likelihood(lbins) + mo.log_likelihood(rbins);
            double aic_break = 2 * k_break - 2 * ll_break;

            //cout << ll_break << ' ' << ll_segment << ' ' << aic_break << ' ' << aic_segment << endl;
            double aic_p = aic_segment - aic_break;
            //cout << aic_p << endl;



            aic_cell.push_back(aic_p);
        }
        aic_vec.push_back(aic_cell);

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

//    if (isinf(res))
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

    vector<double> row1(aic_vec.size(), 0.0);
    vector<double> row2;
    vector<double> res(1,0.0);

    // iterate over n_cells
    for (int i = 0; i < aic_vec.size(); ++i) {
        row2.clear();
        // inner computation of row2
        for (int j = 0; j < aic_vec.size()-i; ++j) {
            double value;
            if (j==0)
            {
                value = row1[j] + aic_vec[j+i];
            }
            else
            {
                value = log_add(row2[j-1] , row1[j] + aic_vec[j+i]);

            }

            if (isinf(value))
                cout << "inf value detected";

            row2.push_back(value);

        }
        //cout << " row1 size: " << row1.size() << " row2 size: " << row2.size() << " max j: " << aic_vec.size()-i<<endl;
        row1.clear();
        row1 = row2;
        double last = row2.back();


        res.push_back(last);

    }



    return res;
}
