//
// Created by Tuncel  Mustafa Anil on 8/30/18.
//

#include "validation.h"
#include <iostream>
#include "SingletonRandomGenerator.h"
#include "globals.cpp"
#include "stats.hpp"


// globals
int print_precision;
int copy_number_limit;
double lambda_s;
double lambda_r;
double lambda_c;

// endof globals


int main()
{

    // test for statslib

    std::mt19937_64 engine(42);
    double ran_val_2 = stats::rnorm(1,2,engine);


    print_precision = 16;
    copy_number_limit = 5;
    lambda_s = 0.5;
    lambda_r = 0.1;
    lambda_c = 0.2;
    std::cout<<"UNIT TESTS" <<std::endl;
//     set a seed number for reproducibility
    SingletonRandomGenerator::get_instance(42);

    test_mathop();
    test_xxhash();
    test_swap_label();
    test_weighted_sample();
    test_prune_reattach();
    test_weighted_prune_reattach();
    test_add_remove_event();
    test_insert_delete_weights();
    test_condense_split_weights(true);
    test_condense_split_weights(false);
    test_reproducibility();
    return EXIT_SUCCESS;
}
