//
// Created by Tuncel  Mustafa Anil on 8/30/18.
//

#include "validation.h"
#include <iostream>
#include "SingletonRandomGenerator.h"



int main()
{
    std::cout<<"\n UNIT TESTS \n ";
//     set a seed number for reproducibility
    SingletonRandomGenerator::get_generator(42);

    test_xxhash();
    test_swap_label();
    test_weighted_sample();
    test_prune_reattach();
    test_weighted_prune_reattach();
    test_add_remove_event();
    test_reproducibility_five_moves();
    return EXIT_SUCCESS;
}