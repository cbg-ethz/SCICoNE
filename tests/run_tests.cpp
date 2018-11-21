//
// Created by Tuncel  Mustafa Anil on 8/30/18.
//

#include "validation.h"
#include <iostream>
#include "SingletonRandomGenerator.h"
#include "globals.cpp"

// globals
int print_precision;
double lambda_s;
double lambda_r;
double lambda_c;

// endof globals


int main()
{
    print_precision = 16;
    lambda_s = 0.5;
    std::cout<<"UNIT TESTS" <<std::endl;
//     set a seed number for reproducibility
    SingletonRandomGenerator::get_generator(42);
    test_mathop();
    test_xxhash();
    test_swap_label();
    test_weighted_sample();
    test_prune_reattach();
    test_weighted_prune_reattach();
    test_add_remove_event();
    test_insert_delete_weights();
    test_condense_split_weights();
//    test_reproducibility();
    return EXIT_SUCCESS;
}
