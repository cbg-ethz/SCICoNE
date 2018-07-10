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

using namespace std;

class MathOp {

public:
    MathOp();
    virtual ~MathOp();
    double log_likelihood(vector<double>);

    template <class T>
    double vec_avg(vector<T>);

    template <class T>
    double vec_sum(vector<T>);



};


#endif //SC_DNA_MATHOP_H
