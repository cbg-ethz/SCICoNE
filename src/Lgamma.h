//
// Created by Tuncel  Mustafa Anil on 11/30/18.
//

#ifndef SC_DNA_LGAMMA_H
#define SC_DNA_LGAMMA_H

#include <cmath>
#include <vector>

class Lgamma {
    /*
     * Static class for using the lgamma values.
     * */

public:

    static std::vector<double> lgamma_values;
    static double get_val(int index);

};


#endif //SC_DNA_LGAMMA_H
