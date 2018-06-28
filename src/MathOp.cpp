//
// Created by Tuncel  Mustafa Anil on 6/19/18.
//

#include "MathOp.h"

MathOp::MathOp() {}

MathOp::~MathOp() {

}

template<class T>
double MathOp::avg(std::vector<T> v) {
    float average = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return average;
}

double MathOp::log_likelihood(std::vector<double> v)
{
    /*double ln_gamma = 0;
    for (auto const &i : v)
    {
        ln_gamma += log(tgamma(i));
    }*/
    // max likelihood:  std::log(max_ll_val) * sum(v) - (v.size() * max_ll_val)
    double max_ll_val = avg(v);
    double term1,term2;
    // to avoid log(0) * 0
    if (sum(v) == 0 && max_ll_val==0)
        term1 = 0.0;
    else
        term1 = std::log(max_ll_val) * sum(v);

    term2 = (v.size() * max_ll_val);
    double ll =  term1 - term2;

    assert(!isnan(ll));
    return ll;
}

template<class T>
double MathOp::sum(std::vector<T> v) {
    return accumulate( v.begin(), v.end(), 0.0);
}
