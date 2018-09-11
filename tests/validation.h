//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_VALIDATION_H
#define SC_DNA_VALIDATION_H


#include "Inference.h"
#include <vector>
#include "xxhash.h"

using namespace std;

// number of cells
const int m = 5;
// counts per region per cell
const vector<vector<double>> D = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

// region sizes
const vector<int> r = {4,2,3,5,2};

// error tolerance
const double epsilon = 0.001;
const double epsilon_sens = 1e-06;

u_int ploidy = 2;



void test_xxhash()
{

    vector<int> r1 = {4,2,3,5,2};
    vector<int> r2 = {4,2,3,5};
    vector<int> r3 = {4,2,3,5,1};

    r2.push_back(2);

    uint64_t const r1_hash = XXH64(&r1[0], size(r1) * sizeof(r1[0]) , 0); // seed = 0
    uint64_t const r2_hash = XXH64(&r2[0], size(r2) * sizeof(r2[0]), 0); // seed = 0
    uint64_t const r3_hash = XXH64(&r3[0], size(r3) * sizeof(r3[0]), 0); // seed = 0

    assert(r1_hash == r2_hash);
    assert(r1_hash != r3_hash);

    cout<<"xxhash validation test passed!"<<endl;

}

// TODO write a reprodubilicity test case with 7 moves
// TODO implement the validation tests on the worked examples

void test_reproducibility_five_moves()
{

    // if seed is not set, set it to 42
    SingletonRandomGenerator::get_generator(42);

    Inference mcmc(size(r));

    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    // move probabilities
    const vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f};
    mcmc.infer_mcmc(D, r, move_probs, 500);

    assert(abs(mcmc.best_tree.score + 2615.9176) <= epsilon);
    cout<<"Reproducibility test with 5 moves and 500 iterations is passed!"<<endl;

}

void test_swap_label()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_swap(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] + 552.120)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 413.462)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 670.394)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.325)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.635)  <= epsilon);

    // compute the log posterior
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score + 2625.580) <= epsilon);

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
    vector<Node*>::const_iterator first = mcmc.t.all_nodes_vec.begin() + 1;
    vector<Node*>::const_iterator last = mcmc.t.all_nodes_vec.end();
    vector<Node*> nodes_to_sample(first, last);

    vector<float> weights;
    float zeta = 0.0f;
    for (auto const &x : nodes_to_sample)
    {
        float weight = (1.0f / x->n_descendents); // weights are inversely proportional to n_descendents
        weights.push_back(weight);
        zeta += weight;
    }

    assert(abs(weights[0] - 0.2)  <= epsilon);
    assert(abs(weights[1] - 0.333)  <= epsilon);
    assert(abs(weights[2] - 1.0)  <= epsilon);
    assert(abs(weights[3] - 1.0)  <= epsilon);
    assert(abs(weights[4] - 1.0)  <= epsilon);
    assert(abs(zeta - 3.533) <= epsilon);


    cout<<"Weighted sample validation test passed!"<<endl;


}

void test_insert_delete_weights()
{
    /*
     * Validation test for the weights of the insert_delete node move.
     * */

    Inference mcmc(std::size(r));
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;

    Tree t(ploidy, std::size(r));
    t.random_insert({{0, 1}, {1, 1}}); // 1
    t.insert_at(1,{{3,1}});  // 2
    t.insert_at(2,{{1, 1}, {2, 1}}); // 3
    t.insert_at(3,{{0, -1}}); // 4
    t.insert_at(3,{{3, -1}}); // 5
    t.insert_at(1,{{1, 1}}); // 6

    t.compute_weights();

    int n = static_cast<int>(D.size());
    for (int i = 0; i < n; ++i)
    {
        t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = t.get_children_id_score(t.root);

        t_scores.push_back(scores_vec);
        t_sums.push_back(MathOp::log_sum(scores_vec));
    }

    int m = size(D);
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.score = mcmc.log_posterior(t_sum, m, t);

    // check the score
    assert(abs(t.score + 2632.658) <= epsilon);

    vector<double> omega; // delete weights
    vector<double> upsilon; // cost weighted omega
    vector<double> chi; // add weights;
    vector<double> xi; // cost weighted chi;


    int K = t.n_regions;
    double lambda_r = 2.0;
    double lambda_c = 1.0;

    vector<Node*> all_nodes = t.root->get_descendents(false); // without root
    for(auto &node: all_nodes)
    {
        //int r_i = static_cast<int>(node->c_change.size());
        double omega_val = MathOp::compute_omega(node, lambda_r, lambda_c, K);
        omega.push_back(omega_val);
        upsilon.push_back(omega_val/node->n_descendents);
    }

    for(auto &node: t.all_nodes_vec)
    {
        double chi_val = pow(2, node->get_n_children()); // chi is to be computed for the n_first order children
        chi.push_back(chi_val);
        double xi_val = pow(2, node->get_n_children()+1) / (node->n_descendents+1);
        xi.push_back(xi_val);
    }

    double sum_chi = MathOp::vec_sum(chi);
    double sum_xi = MathOp::vec_sum(xi);


    assert(abs(omega[0] - 0.916e-03) <=epsilon_sens);
    assert(abs(omega.back() - 4.979e-03) <=epsilon_sens);
    assert(abs(sum_chi - 15.0) <= epsilon);
    assert(abs(sum_xi - 7.443) <= epsilon);




    cout<<"Insert and delete node weights validation test passed!"<<endl;

}

void test_prune_reattach()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_prune_reattach(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] + 551.804)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 413.292)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 663.804)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.348)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.636)  <= epsilon);

    assert(abs(mcmc.t_prime.score + 2615.918) <= epsilon);
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score + 2618.528) <= epsilon);


    cout<<"Prune and reattach validation test passed!"<<endl;
}

void test_weighted_prune_reattach()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_prune_reattach(D,r,false,true);


    // get the subvector
    vector<Node*>::const_iterator first = mcmc.t_prime.all_nodes_vec.begin() + 1;
    vector<Node*>::const_iterator last = mcmc.t_prime.all_nodes_vec.end();
    vector<Node*> nodes_to_sample(first, last);

    // re-ordering is needed since prune and reattach does not preserve the order in the all_nodes vector
    std::sort(nodes_to_sample.begin(),nodes_to_sample.end(), [](Node* a, Node* b) { return *a < *b; });



    vector<double> weights;
    double zeta = 0.0;
    for (auto const &x : nodes_to_sample)
    {
        float weight = (1.0 / x->n_descendents); // weights are inversely proportional to n_descendents
        weights.push_back(weight);
        zeta += weight;
    }



    assert(abs(weights[0] - 0.2)  <= epsilon);
    assert(abs(weights[1] - 0.333)  <= epsilon);
    assert(abs(weights[2] - 1.0)  <= epsilon);
    assert(abs(weights[3] - 1.0)  <= epsilon);
    assert(abs(weights[4] - 0.25)  <= epsilon);
    assert(abs(zeta - 2.783) <= epsilon);


    cout<<"Weighted prune &reattaach validation test passed!"<<endl;

}

void test_add_remove_event()
{

    Inference mcmc(std::size(r));
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_add_remove_events(0.0, 0.0, D, r, false, true);

    assert(abs(mcmc.t_prime_sums[0] + 553.134)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] + 411.720)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] + 663.925)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] + 547.325)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] + 406.654)  <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score + 2620.015) <= epsilon);

    cout<<"Add / remove event validation test passed!"<<endl;
}


#endif //SC_DNA_VALIDATION_H
