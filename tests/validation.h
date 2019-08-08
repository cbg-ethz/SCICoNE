//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_VALIDATION_H
#define SC_DNA_VALIDATION_H


#include "Inference.h"
#include <vector>
#include <numeric>
#include "xxhash.h"

using namespace std;

// number of cells
const int m = 5;
// counts per region per cell
const vector<vector<double>> D = {
    {39, 37, 45, 49, 30},
    {31, 28, 34, 46, 11},
    {69, 58, 68, 34, 21},
    {72, 30, 31, 46, 21},
    {64, 17, 19, 37, 13}
};

// region sizes
const vector<int> r = {4,2,3,5,2};

// error tolerance
const double epsilon = 0.001;
const double epsilon_sens = 1e-06;

unsigned ploidy = 2;

void test_ploidy_attachment_score()
{
    /*
     * Checks if two nodes on different trees with same ploidy have the same attachment scores
     * */

    vector<int> reg_sizes = {4,2};

    Tree t_two(2, reg_sizes.size());
    t_two.insert_at(0,{{0, 1}, {1,1}});
    t_two.compute_weights();

    Tree t_three(3, reg_sizes.size());
    t_three.compute_weights();

    vector<double> counts = {10, 40, 80, 120, 500};

    t_two.compute_tree(counts, reg_sizes);
    std::map<int, double> scores_vec = t_two.get_children_id_score(t_two.root);
    t_three.compute_tree(counts, reg_sizes);
    std::map<int, double> scores_vec_prime = t_three.get_children_id_score(t_three.root);

    assert(scores_vec[0]==scores_vec_prime[0]);
    cout<<"Cell attachment with different ploidies validation test passed!"<<endl;

}

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

    /*
     * Tests the functionality of the xxhash module
     * */

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

void test_reproducibility()
{
    /*
     * Tests the reproducibility of the markov chain
     * */

    // if seed is not set, set it to 42
    SingletonRandomGenerator::get_instance(42);

    Inference mcmc(r.size(), ploidy, verbosity);

    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);
    mcmc.update_t_prime(); // set t_prime to t

    // move probabilities
    vector<float> move_probs = {0.0f,1.0f,0.0f,1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f};
    //---------------------------pr--w-pr--sw--w-sw---ar----w-ar--id---w-id---cs---w-cs--geno--

    unsigned size_limit = std::numeric_limits<unsigned>::max();
    mcmc.infer_mcmc(D, r, move_probs, 2500, size_limit);

    cout<<"Reproducibility score: " << mcmc.best_tree.posterior_score << endl;
    cout<<"Epsilon: " << epsilon << endl;
    assert(abs(mcmc.best_tree.posterior_score - 30.653) <= epsilon);
    cout<<"Reproducibility test is passed!"<<endl;

}

void test_swap_label()
{
    /*
     * Tests the swap label move
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);
    mcmc.update_t_prime(); // set t_prime to t

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });


    mcmc.apply_swap(D,r,false,true);

    assert(abs(mcmc.t_prime_sums[0] - 3.027)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 2.709)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 29.222)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 7.445)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 7.415)  <= epsilon);

    // compute the root score
    size_t n_cells = D.size();
    double sum_root_score = 0.0;
    for (u_int i = 0; i < n_cells; ++i) {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        double root_score = mcmc.t_prime.get_od_root_score(r,sum_d,D[i]);
        sum_root_score += root_score;
    }

    // compute the log posterior
    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_tree_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.posterior_score = log_post_t_prime;

    double total_score = mcmc.t_prime.posterior_score + sum_root_score;
    double total_score_gt = -1503.33;
    assert(abs(total_score - total_score_gt) <= epsilon);

    std::cout<<"Swap label validation test passed!"<<std::endl;


}

void test_weighted_sample()
{
    /*
     * Validation test for the weighted sampling
     *
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
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
    /*
     * Tests the condense split weights
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    std::vector<std::map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;

    Tree t(ploidy, r.size());
    t.random_insert({{0, 1}, {1, 1}}); // 1
    t.insert_at(1,{{1, 1}, {2, 1}}); // 2
    t.insert_at(2,{{0, -1}}); // 3
    t.insert_at(2,{{3, -1}}); // 4
    t.insert_at(1,{{0, 1}}); // 5

    t.compute_weights();


    Tree t_prime(ploidy, r.size());
    t_prime.random_insert({{0, 1}, {1, 1}}); // 1
    t_prime.insert_at(1,{{2, 1}}); // 2
    t_prime.insert_at(2,{{1, 1}}); // 3
    t_prime.insert_at(3,{{0, -1}}); // 4
    t_prime.insert_at(2,{{3, -1}}); // 5
    t_prime.insert_at(1,{{0, 1}}); // 6

    t_prime.compute_weights();

    size_t n_cells = static_cast<int>(D.size());
    for (size_t i = 0; i < n_cells; ++i)
    {
        t_prime.compute_tree(D[i], r);
        std::map<int, double> scores_vec_prime = t_prime.get_children_id_score(t_prime.root);
        t_prime_scores.push_back(scores_vec_prime);
        t_prime_sums.push_back(MathOp::log_sum(scores_vec_prime));
    }
    double sum_root_score = 0.0;
    for (size_t i = 0; i < n_cells; ++i)
    {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        double root_score = mcmc.t_prime.get_od_root_score(r,sum_d,D[i]);
        sum_root_score += root_score;
    }

    double event_prior = t_prime.event_prior();
    double event_prior_tp_gt = -30.395;
    assert(abs(event_prior - event_prior_tp_gt) <= epsilon);


    double t_prime_sum = accumulate( t_prime_sums.begin(), t_prime_sums.end(), 0.0);
    t_prime.posterior_score = mcmc.log_tree_posterior(t_prime_sum, m, t_prime);

    double total_score_tp = sum_root_score + t_prime.posterior_score;
    double total_score_tp_gt = -1507.806;
    assert(abs(total_score_tp - total_score_tp_gt) <= epsilon);

    // intentional override
    double lambda_s = 0.5;

    double sum_chi = t.chi_condense_split_reweighted(false); // weighted = false;
    double sum_chi_gt = 0.0595;
    assert(abs(sum_chi - sum_chi_gt) <= epsilon);

    vector<double> omega = t.omega_condense_split(lambda_s, false);
    double sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
    double sum_omega_gt = 0.296;
    assert(abs(sum_omega - sum_omega_gt) <= epsilon);

    double sum_chi_prime = t_prime.chi_condense_split_reweighted(false);
    double sum_chi_tp_gt = 0.0198;
    assert(abs(sum_chi_prime - sum_chi_tp_gt) <= epsilon);

    vector<double> omega_prime = t_prime.omega_condense_split(lambda_s, false);
    double sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
    double sum_omega_prime_gt = 0.488;
    assert(abs(sum_omega_prime - sum_omega_prime_gt) <= epsilon);

    // weighted versions
    double sum_xi = t.chi_condense_split_reweighted(true);
    double sum_xi_gt = 0.0124;
    assert(abs(sum_xi - sum_xi_gt) <= epsilon);

    vector<double> upsilon = t.omega_condense_split(lambda_s, true);
    double sum_upsilon = accumulate(upsilon.begin(), upsilon.end(), 0.0);
    double sum_upsilon_gt = 0.0878;
    assert(abs(sum_upsilon - sum_upsilon_gt) <= epsilon);

    double sum_xi_tp = t_prime.chi_condense_split_reweighted(true);
    double sum_xi_tp_gt = 0.002825;
    assert(abs(sum_xi_tp - sum_xi_tp_gt) <= epsilon);

    vector<double> upsilon_tp = t_prime.omega_condense_split(lambda_s, true);
    double sum_upsilon_tp = accumulate(upsilon_tp.begin(), upsilon_tp.end(), 0.0);
    double sum_upsilon_tp_gt = 0.1956;
    assert(abs(sum_upsilon_tp - sum_upsilon_tp_gt) <= epsilon);

    std::cout<<"Condense and split node weights validation test passed! " << std::endl;

}

void test_insert_delete_weights()
{
    /*
     * Validation test for the weights of the insert_delete node move.
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;

    Tree t(ploidy, r.size());
    t.random_insert({{0, 1}, {1, 1}}); // 1
    t.insert_at(1,{{3,1}});  // 2
    t.insert_at(2,{{1, 1}, {2, 1}}); // 3
    t.insert_at(3,{{0, -1}}); // 4
    t.insert_at(3,{{3, -1}}); // 5
    t.insert_at(1,{{0, 1}}); // 6

    t.compute_weights();

    size_t n_cells = static_cast<int>(D.size());
    for (size_t i = 0; i < n_cells; ++i)
    {
        t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = t.get_children_id_score(t.root);

        t_scores.push_back(scores_vec);
        t_sums.push_back(MathOp::log_sum(scores_vec));
    }
    double sum_root_score = 0.0;
    for (u_int i = 0; i < n_cells; ++i) {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        double root_score = mcmc.t_prime.get_od_root_score(r,sum_d,D[i]);
        sum_root_score += root_score;
    }

    int m = D.size();
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.posterior_score = mcmc.log_tree_posterior(t_sum, m, t);

    double event_prior = t.event_prior();
    double event_prior_gt = -40.395;
    assert(abs(event_prior - event_prior_gt) <= epsilon);

    double total_score = sum_root_score + t.posterior_score;
    double total_score_gt = -1524.053;
    assert(abs(total_score - total_score_gt) <= epsilon);

    // intentional overriding of parameters
    double lambda_r = 2.0;
    double lambda_c = 1.0;

    vector<double> omega = t.omega_insert_delete(lambda_r, lambda_c, false); // delete weights
    assert(abs(omega[0] - 0.916e-03) <=epsilon_sens);
    assert(abs(omega.back() - 4.979e-03) <=epsilon_sens);

    double sum_omega = accumulate(omega.begin(), omega.end(), 0.0);
    double sum_omega_gt = 0.0217;
    assert(abs(sum_omega - sum_omega_gt) <= epsilon);


    double sum_chi = t.chi_insert_delete_reweighted(false); // weighted = false;
    double sum_chi_gt = 0.0742;
    assert(abs(sum_chi - sum_chi_gt) <= epsilon);

    vector<double> upsilon = t.omega_insert_delete(lambda_r, lambda_c, true); // cost weighted omega
    vector<double> xi = t.chi_insert_delete(true); // cost weighted chi;

    double sum_upsilon = accumulate(upsilon.begin(), upsilon.end(), 0.0);
    double sum_upsilon_gt = 0.0166;
    assert(abs(sum_upsilon - sum_upsilon_gt) <= epsilon);

    double sum_xi = t.chi_insert_delete_reweighted(true);
    double sum_xi_gt = 0.0368;

    assert(abs(sum_xi - sum_xi_gt) <= epsilon);


    cout<<"Insert and delete node weights validation test passed!"<<endl;

}

void test_tree_prior()
{
    /*
     * Validation test for the tree prior computation
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    u_int n = mcmc.t.get_n_nodes();


    double tree_prior = mcmc.log_tree_prior(m,n);
    double tree_prior_gt = -16.12583;
    assert(abs(tree_prior - tree_prior_gt) <= epsilon);

    std::cout<<"Tree prior validation test passed!"<<std::endl;

}

void test_event_prior()
{
    /*
     * Validation test for the event prior computation
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    double event_prior = mcmc.t.event_prior();
    double event_prior_gt = -26.3;
    assert(abs(event_prior - event_prior_gt)  <= epsilon);

    std::cout<<"Event prior computation test passed!"<<std::endl;

}

void test_tree_attachment()
{
    /*
     * Tests the correctness of tree score computation methods on worked example
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    std::vector<std::vector<double>> t_scores_gt = {
        {0.0,-3.555, 0.613, 4.327, -8.629, -12.352},
        {0.0, -1.855, 1.444, 3.520, -11.899, -7.949},
        {0.0, 8.532, 21.885, 18.318, 29.222, 3.361},
        {0.0, 6.987, 3.464, -6.201, -3.698, 7.684},
        {0.0, 7.065, 1.117, -10.187, -5.987, 10.464}
    };

    // the natural sum on the log space values
    std::vector<double> t_sums_gt = {4.364, 3.668, 29.222, 8.098, 10.497};

    for (size_t i = 0; i < t_scores_gt.size(); ++i)
        for (size_t j = 0; j < t_scores_gt[0].size(); ++ j)
            assert(abs(mcmc.t_scores[i][j] - t_scores_gt[i][j])  <= epsilon);

    for (size_t i = 0; i < t_sums_gt.size(); ++i)
        assert(abs(mcmc.t_sums[i] - t_sums_gt[i])  <= epsilon);

    std::cout<<"Tree attachment computation test passed!"<<std::endl;
}

void test_prune_reattach()
{
    /*
     * Tests the prune and reattach move
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);
    mcmc.update_t_prime(); // set t_prime to t

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });

    assert(abs(mcmc.t.posterior_score - 13.424) <= epsilon);

    mcmc.apply_prune_reattach(D, r, false, true);

    assert(abs(mcmc.t_prime_sums[0] - 1.057)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 1.695)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 24.169)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 8.264)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 10.511)  <= epsilon);

    // compute the root score
    size_t n_cells = D.size();
    double sum_root_score = 0.0;
    for (u_int i = 0; i < n_cells; ++i) {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        double root_score = mcmc.t_prime.get_od_root_score(r,sum_d,D[i]);
        sum_root_score += root_score;
    }

    double sum_root_score_gt = -1510.724;
    assert(abs(sum_root_score - sum_root_score_gt) <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_tree_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.posterior_score = log_post_t_prime;
    assert(abs(mcmc.t_prime.posterior_score - 3.269) <= epsilon);

    // total score
    double total_score_gt = -1507.455;
    double total_score = mcmc.t_prime.posterior_score + sum_root_score;
    assert(abs(total_score - total_score_gt) <= epsilon);


    cout<<"Prune and reattach validation test passed!"<<std::endl;
}

void test_weighted_prune_reattach()
{

    /*
     * Tests the weighted prune reattach move
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);
    mcmc.update_t_prime(); // set t_prime to t

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });

    mcmc.apply_prune_reattach(D, r, false, true);


    double zeta = mcmc.t.cost();
    double zeta_prime = mcmc.t_prime.cost();

    assert(abs(zeta - 3.533) <= epsilon);
    assert(abs(zeta_prime - 2.783) <= epsilon);

    cout<<"Weighted prune &reattach validation test passed!"<<endl;

}

void test_add_remove_event()
{
    /*
     * Tests the add remove event move
     * */

    Inference mcmc(r.size(), ploidy, verbosity);
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r); // assignment operator changes the order
    mcmc.update_t_prime(); // set t_prime to t

    // re-ordering is needed since the copy_tree method does not preserve the order in the all_nodes vector
    std::sort(mcmc.t_prime.all_nodes_vec.begin(),mcmc.t_prime.all_nodes_vec.end(), [](Node* a, Node* b) { return *a < *b; });

    mcmc.apply_add_remove_events(lambda_r, lambda_c, D, r, false, true);

    assert(abs(mcmc.t_prime_sums[0] - 1.383)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[1] - 4.168)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[2] - 29.222)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[3] - 7.191)  <= epsilon);
    assert(abs(mcmc.t_prime_sums[4] - 9.233)  <= epsilon);

    double t_prime_sum = accumulate( mcmc.t_prime_sums.begin(), mcmc.t_prime_sums.end(), 0.0);
    double log_post_t_prime = mcmc.log_tree_posterior(t_prime_sum, m, mcmc.t_prime);
    mcmc.t_prime.posterior_score = log_post_t_prime;
    assert(abs(mcmc.t_prime.posterior_score - 16.47 + 1*c_penalise) <= epsilon);

    cout<<"Add / remove event validation test passed!"<<endl;
}

void test_children_repeat_genotype()
{
    /*
     * Tests if two siblings repeat the same region sign
     * */

    Tree t(ploidy, r.size());
    t.random_insert({{0, 1}, {1, 1}}); // 1
    t.insert_at(0,{{3,1}});  // 2
    t.insert_at(0,{{1, -2}, {2, 1}}); // 3
    t.insert_at(3,{{0, -1}}); // 4

    assert(!t.root->children_repeat_genotype());

    t.insert_at(0,{{3,-1}, {2,2}}); // 5
    assert(t.root->children_repeat_genotype());

    cout<<"Children repeat genotype validation test passed!"<<endl;
}

void test_tree_validation()
{
    /*
     * Makes sure the valid trees are considered valid while invalid trees are detected
     * */

    u_int n_regions = 0;
    Tree* t = new Tree(ploidy, n_regions);

    t->load_from_file("../tests/trees_to_validate/valid_tree_1.txt");
    assert(t->is_valid_subtree(t->root));
    assert(not t->is_redundant());
    delete t;

    t = new Tree(ploidy, n_regions);
    t->load_from_file("../tests/trees_to_validate/invalid_tree_1.txt");
    assert(not t->is_valid_subtree(t->root)); // NOT VALID, first order children (siblings) repeat genotype
    assert(not t->is_redundant());
    delete t;

    t = new Tree(ploidy, n_regions);
    t->load_from_file("../tests/trees_to_validate/invalid_tree_2.txt");
    assert(t->is_valid_subtree(t->root));
    assert(t->is_redundant()); // IS REDUNDANT, two nodes carry the same genotype
    delete t;

    t = new Tree(ploidy, n_regions);
    t->load_from_file("../tests/trees_to_validate/invalid_tree_3.txt");
    assert(t->zero_ploidy_changes(t->root)); // NOT VALID, zero ploidy changes back!
    assert(not t->is_redundant());
    delete t;

    t = new Tree(ploidy, n_regions);
    t->load_from_file("../tests/trees_to_validate/invalid_tree_4.txt");
    assert(t->subtree_out_of_bound(t->root)); // NOT VALID, tree is out of bounds!
    assert(t->zero_ploidy_changes(t->root)); // NOT VALID, zero ploidy changes back!
    assert(not t->is_redundant());
    delete t;

    std::cout<<"Is valid tree? validation test passed!"<<std::endl;
}

#endif //SC_DNA_VALIDATION_H
