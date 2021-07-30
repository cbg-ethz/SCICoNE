//
// Created by Tuncel  Mustafa Anil on 8/30/18.
//

#include "validation.h"
#include <iostream>
#include "SingletonRandomGenerator.h"
#include "globals.cpp"

#include <Config.h>

// globals
int print_precision;
int copy_number_limit;
double lambda_s;
double lambda_r;
double lambda_c;
double cf;
double c_penalise;
unsigned is_overdispersed;
string f_name_posfix;
int verbosity;
double eta;
// endof globals


int main()
{
    std::cout << "Hey! " << SOURCE_DIR << std::endl;
    std::string bp_detection_dir = std::string(SOURCE_DIR) + "/tests/bp_detection";
    std::string trees_to_validate_dir = std::string(SOURCE_DIR) + "/tests/trees_to_validate";

    print_precision = 15;
    copy_number_limit = 5;
    lambda_s = 0.5;
    lambda_r = 0.1;
    lambda_c = 0.2;
    cf = 0.0;
    c_penalise = 10.0;
    is_overdispersed = 0;
    verbosity = 0;
    eta = 1e-4;

    std::cout<<"UNIT TESTS" <<std::endl;
//     set a seed number for reproducibility
    SingletonRandomGenerator::get_instance(42);
    test_mathop();
    test_xxhash();
    test_breakpoint_detection(bp_detection_dir);
    test_tree_prior();
    test_event_prior();
    test_n_descendents_computation();
    test_tree_attachment();
    test_prune_reattach();
    test_swap_label();
    test_ploidy_attachment_score();
    test_weighted_sample();
    test_weighted_prune_reattach();
    test_add_remove_event();
    test_insert_delete_weights();
    test_condense_split_weights();
    test_tree_validation(trees_to_validate_dir);
    test_overdispersed_score();
    test_genotype_preserving_move_scores();
    test_apply_multiple_times();
    test_reproducibility();
    test_cluster_scoring();
    test_expand_shrink_block_operation();
    test_expand_shrink_block();
    test_get_event_intersection();
    test_add_common_ancestor_operation();
    test_add_common_ancestor();

    return EXIT_SUCCESS;

}
