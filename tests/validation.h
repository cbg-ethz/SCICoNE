//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_VALIDATION_H
#define SC_DNA_VALIDATION_H


#include "Inference.h"
#include <vector>

using namespace std;

// number of cells
const int m = 5;
// counts per region per cell
const vector<vector<int>> D = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

// region sizes
const vector<int> r = {4,2,3,5,2};

// error tolerance
const float epsilon = 0.001f;

void test_swap_label()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_swap(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] + 552.120f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 413.462f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 670.394f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.325f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.635f)  <= epsilon);

    // compute the log posterior
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, *mcmc.t_prime);
    mcmc.t_prime->score = log_post_t_prime;
    assert(abs(mcmc.t_prime->score + 2625.580f) <= epsilon);

    cout<<"Swap label validation test passed!"<<endl;


}

void test_weighted_sample()
{
    /*
     * Validation test for the weighted sampling
     *
     * */

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();

    // get the subvector
    vector<Node*>::const_iterator first = mcmc.t->all_nodes.begin() + 1;
    vector<Node*>::const_iterator last = mcmc.t->all_nodes.end();
    vector<Node*> nodes_to_sample(first, last);

    vector<float> weights;
    float zeta = 0.0f;
    for (auto const &x : nodes_to_sample)
    {
        float weight = (1.0f / x->n_descendents); // weights are inversely proportional to n_descendents
        weights.push_back(weight);
        zeta += weight;
    }

    assert(abs(weights[0] - 0.2f)  <= epsilon);
    assert(abs(weights[1] - 0.333f)  <= epsilon);
    assert(abs(weights[2] - 1.0f)  <= epsilon);
    assert(abs(weights[3] - 1.0f)  <= epsilon);
    assert(abs(weights[4] - 1.0f)  <= epsilon);
    assert(abs(zeta - 3.533f) <= epsilon);


    cout<<"Weighted sample validation test passed!"<<endl;


}

void test_prune_reattach()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_prune_reattach(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] + 551.804f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 413.292f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 663.804f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.348f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.636f)  <= epsilon);

    assert(abs(mcmc.t_prime->score + 2615.918f) <= epsilon);
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, *mcmc.t_prime);
    mcmc.t_prime->score = log_post_t_prime;
    assert(abs(mcmc.t_prime->score + 2618.528f) <= epsilon);


    cout<<"Prune and reattach validation test passed!"<<endl;
}

void test_weighted_prune_reattach()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_prune_reattach(D,r,false,true);


    // get the subvector
    vector<Node*>::const_iterator first = mcmc.t_prime->all_nodes.begin() + 1;
    vector<Node*>::const_iterator last = mcmc.t_prime->all_nodes.end();
    vector<Node*> nodes_to_sample(first, last);

    vector<float> weights;
    float zeta = 0.0f;
    for (auto const &x : nodes_to_sample)
    {
        float weight = (1.0f / x->n_descendents); // weights are inversely proportional to n_descendents
        weights.push_back(weight);
        zeta += weight;
    }


    assert(abs(weights[0] - 0.2f)  <= epsilon);
    assert(abs(weights[1] - 0.333f)  <= epsilon);
    assert(abs(weights[2] - 1.0f)  <= epsilon);
    assert(abs(weights[3] - 1.0f)  <= epsilon);
    assert(abs(weights[4] - 0.25f)  <= epsilon);
    assert(abs(zeta - 2.783f) <= epsilon);


    cout<<"Weighted prune &reattaach validation test passed!"<<endl;

}

void test_add_remove_event()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_add_remove_events(0.0f, 0.0f, D, r, false, true);

    assert(abs(mcmc.t_prime_sums[0] + 553.134f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 411.720f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 663.925f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.325f)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.654f)  <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, *mcmc.t_prime);
    mcmc.t_prime->score = log_post_t_prime;
    assert(abs(mcmc.t_prime->score + 2620.015f) <= epsilon);

    cout<<"Add / remove event validation test passed!"<<endl;
}


#endif //SC_DNA_VALIDATION_H
