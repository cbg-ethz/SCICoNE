//
// Created by Tuncel  Mustafa Anil on 11/30/18.
//

#include "Lgamma.h"

// define the static variables
std::vector<double> Lgamma::lgamma_values;

double Lgamma::get_val(unsigned long index) {
    if (lgamma_values.size() > index)
    {
        return lgamma_values[index];
    }
    else
    {

        if (lgamma_values.empty())
        {
            double inf = std::numeric_limits<double>::infinity();
            lgamma_values = {inf, 0,0};
        }


        while (lgamma_values.size()-1 < index)
        {
            double last_val = lgamma_values[lgamma_values.size()-1];
            lgamma_values.push_back(last_val + log(lgamma_values.size()-1));
        }
        return lgamma_values[index];

    }
}
