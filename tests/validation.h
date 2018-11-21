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


void test_mathop()
{
    /*
     * Tests for mathop functions
     * */
    vector<double> vals1 = {12.5, 5.2, 10, 3.8, 15.5};
    vector<double> vals2 = {12.5, 5.2, 10, 44.1, 3.8, 19.5};

    double std = MathOp::st_deviation(vals1);
    double med1 = MathOp::median(vals1);
    double med2 = MathOp::median(vals2);

    assert(abs(std - 4.3858) <= epsilon);
    assert(abs(med1 - 10.0) <= epsilon);
    assert(abs(med2 - 11.25) <= epsilon);

}

void test_xxhash()
{

    vector<int> r1 = {4,2,3,5,2};
    vector<int> r2 = {4,2,3,5};
    vector<int> r3 = {4,2,3,5,1};

    r2.push_back(2);

    uint64_t const r1_hash = XXH64(&r1[0], r1.size() * sizeof(r1[0]) , 0); // seed = 0
    uint64_t const r2_hash = XXH64(&r2[0], r2.size() * sizeof(r2[0]), 0); // seed = 0
    uint64_t const r3_hash = XXH64(&r3[0], r3.size() * sizeof(r3[0]), 0); // seed = 0

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

    Inference mcmc(r.size());

    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    // move probabilities
    const vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f};
    mcmc.infer_mcmc(D, r, move_probs, 500, 0);

    assert(abs(mcmc.best_tree.score + 2615.9176) <= epsilon);
    cout<<"Reproducibility test with 5 moves and 500 iterations is passed!"<<endl;

}

void test_swap_label()
{

    Inference mcmc(r.size());
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r); // t prime is copied from t, order is changed

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });


    mcmc.apply_swap(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] - 2.397)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 2.426)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 22.753)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 7.192)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 9.253)  <= epsilon);

    // compute the log posterior
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score - 11.597) <= epsilon);

    cout<<"Swap label validation test passed!"<<endl;


}

void test_weighted_sample()
{
    /*
     * Validation test for the weighted sampling
     *
     * */

    Inference mcmc(r.size());
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

void test_condense_split_weights()
{
    Inference mcmc(r.size());
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;
    std::vector<std::map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;

    Tree t(ploidy, r.size());
    t.random_insert({{0, 1}, {1, 1}}); // 1
    t.insert_at(1,{{1, 1}, {2, 1}}); // 2
    t.insert_at(2,{{0, -1}}); // 3
    t.insert_at(2,{{3, -1}}); // 4
    t.insert_at(1,{{1, 1}}); // 5

    t.compute_weights();


    Tree t_prime(ploidy, r.size());
    t_prime.random_insert({{0, 1}, {1, 1}}); // 1
    t_prime.insert_at(1,{{2, 1}}); // 2
    t_prime.insert_at(2,{{1, 1}}); // 3
    t_prime.insert_at(3,{{0, -1}}); // 4
    t_prime.insert_at(2,{{3, -1}}); // 5
    t_prime.insert_at(1,{{1, 1}}); // 6

    t_prime.compute_weights();

    int n = static_cast<int>(D.size());
    for (int i = 0; i < n; ++i)
    {
        t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = t.get_children_id_score(t.root);
        t_prime.compute_tree(D[i], r);
        std::map<int, double> scores_vec_prime = t_prime.get_children_id_score(t_prime.root);
        t_scores.push_back(scores_vec);
        t_sums.push_back(MathOp::log_sum(scores_vec));
        t_prime_scores.push_back(scores_vec_prime);
        t_prime_sums.push_back(MathOp::log_sum(scores_vec_prime));
    }

    int m = D.size();
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.score = mcmc.log_posterior(t_sum, m, t);

    double t_prime_sum = accumulate( t_prime_sums.begin(), t_prime_sums.end(), 0.0);
    t_prime.score = mcmc.log_posterior(t_prime_sum, m, t_prime);
    // check the scores
    assert(abs(t.score + 2612.6987) <= epsilon);
    assert(abs(t_prime.score + 2623.16586) <= epsilon);

    bool weighted = true;
    double lambda_s = 0.5;

    vector<double> chi = t.chi_condense_split(weighted);
    double sum_chi = std::accumulate(chi.begin(), chi.end(), 0.0);
    if (weighted)
        assert(abs(sum_chi -1.667) <= epsilon);
    else
        assert(abs(sum_chi - 8) <= epsilon);

    vector<double> omega = t.omega_condense_split(lambda_s, weighted);
    double sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
    if (weighted)
        assert(abs(sum_omega - 0.0878) <= epsilon);
    else
        assert(abs(sum_omega - 0.296) <= epsilon);

    vector<double> chi_prime = t_prime.chi_condense_split(weighted);
    double sum_chi_prime = std::accumulate(chi_prime.begin(), chi_prime.end(), 0.0);
    if (weighted)
        assert(abs(sum_chi_prime - 0.571) <= epsilon);
    else
        assert(abs(sum_chi_prime - 4) <= epsilon);

    vector<double> omega_prime = t_prime.omega_condense_split(lambda_s, weighted);
    double sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
    if (weighted)
        assert(abs(sum_omega_prime - 0.1956) <= epsilon);
    else
        assert(abs(sum_omega_prime - 0.488) <= epsilon);


    double score_diff = t_prime.score - t.score;
    assert(abs(score_diff + 10.467) <= epsilon);

//    double acceptance_prob = exp(score_diff) * (sum_chi+sum_omega) / (sum_chi_prime + sum_omega_prime);
//    if (weighted)
//        assert(abs(acceptance_prob - 0.000244) <= epsilon_sens);
//    else
//        assert(abs(acceptance_prob - 0.000263) <= epsilon_sens);

    cout<<"Condense and split node weights validation test passed!"<<endl;


}

void test_insert_delete_weights()
{
    /*
     * Validation test for the weights of the insert_delete node move.
     * */

    Inference mcmc(r.size());
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;

    Tree t(ploidy, r.size());
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

    int m = D.size();
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.score = mcmc.log_posterior(t_sum, m, t);

    // check the score
    //assert(abs(t.score + 2632.658) <= epsilon);

    int K = t.n_regions;
    double lambda_r = 2.0;
    double lambda_c = 1.0;

    vector<double> omega = t.omega_insert_delete(lambda_r, lambda_c, false); // delete weights
    vector<double> upsilon = t.omega_insert_delete(lambda_r, lambda_c, true); // cost weighted omega
    vector<double> chi = t.chi_insert_delete(false); // weighted = false;
    vector<double> xi = t.chi_insert_delete(true); // cost weighted chi;

    double sum_chi = accumulate(chi.begin(), chi.end(), 0.0);
    double sum_xi = accumulate(xi.begin(), xi.end(), 0.0);

    assert(abs(omega[0] - 0.916e-03) <=epsilon_sens);
    assert(abs(omega.back() - 4.979e-03) <=epsilon_sens);
    assert(abs(sum_chi - 15.0) <= epsilon);
    assert(abs(sum_xi - 7.443) <= epsilon);


    cout<<"Insert and delete node weights validation test passed!"<<endl;

}

void test_prune_reattach()
{

    Inference mcmc(r.size());
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });

    assert(abs(mcmc.t.score - 21.26) <= epsilon);

    mcmc.apply_prune_reattach(D, r, false, false, true);

    assert(abs(mcmc.t_prime_sums[0] - 2.713)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 2.596)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 29.343)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 7.169)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 9.252)  <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score - 18.649) <= epsilon);


    cout<<"Prune and reattach validation test passed!"<<endl;
}

void test_weighted_prune_reattach()
{

    Inference mcmc(r.size());
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });

    mcmc.apply_prune_reattach(D, r, false, false, true);


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

    Inference mcmc(r.size());
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r); // assignment operator changes the order

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });


    mcmc.apply_add_remove_events(0.0, 0.0, D, r, false, true);

    assert(abs(mcmc.t_prime_sums[0] - 1.383)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 4.168)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 29.222)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 7.191)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 9.233)  <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.score = log_post_t_prime;
    assert(abs(mcmc.t_prime.score - 16.47) <= epsilon);

    cout<<"Add / remove event validation test passed!"<<endl;
}


#endif //SC_DNA_VALIDATION_H
