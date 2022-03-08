//
// Created by Tuncel  Mustafa Anil on 7/17/18.
//

#ifndef SC_DNA_INFERENCE_H
#define SC_DNA_INFERENCE_H

#include "Tree.h"
#include "SingletonRandomGenerator.h"
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <array>
#include <functional>
#include <cfloat>
#include "globals.cpp"
#include "Lgamma.h"

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

class Inference {
/*
 * Contains functionality to perform Markov chain Monte Carlo (mcmc) inference
 *
 * */
public:
    // TODO: store n_cells as well and use whenever required
    Tree t;
    Tree t_prime;
    Tree best_tree;
    u_int n_regions;
    int ploidy;
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;
    std::vector<std::pair<int, double>> t_maxs;
    std::vector<std::map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;
    std::vector<std::pair<int, double>> t_prime_maxs;
    int verbosity;
    bool max_scoring;
    bool smoothed;
    int n_cells;
    std::vector<int> region_neutral_states;

public:
    Inference(u_int n_regions, int n_cells, vector<int> &region_neutral_states, int ploidy=2, int verbosity=2, bool max_scoring=false, bool smoothed=false);
    ~Inference();
    void destroy();
    void compute_t_root_scores(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
    void compute_t_prime_root_scores(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
    std::vector<double>
    get_tree_root_scores(const vector<vector<double>> &D, const vector<int> &r, const Tree &tree);
    void compute_t_table(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
    void compute_t_prime_scores(Node *attached_node, const vector<vector<double>> &D, const vector<int> &r);
    void compute_t_prime_sums(const vector<vector<double>> &D);
    void update_t_prime();

    double log_tree_prior(int m, int n);
    double log_tree_posterior(double tree_sum, int m, Tree &tree);
    template<class ...Ts, class AnyFunction>
    bool apply_multiple_times(unsigned n, AnyFunction func, Ts &...args);
    bool apply_prune_reattach(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                              bool validation_test_mode);
    bool apply_genotype_preserving_pr(double gamma);
    bool apply_delete_leaf(const vector<vector<double>> &D, const vector<int> &r);
    bool apply_add_common_ancestor(const vector<vector<double>> &D, const vector<int> &r, bool validation_test_mode);
    bool apply_add_remove_events(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                                 bool validation_test_mode);
    bool apply_expand_shrink_blocks(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                                bool validation_test_mode);
    bool apply_insert_delete_node(const vector<vector<double>> &D, const vector<int> &r, unsigned int size_limit,
                                  bool weighted);
    bool apply_condense_split(const vector<vector<double>> &D, const vector<int> &r, unsigned int size_limit,
                              bool weighted);
    bool apply_swap(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                    bool test_mode);
    bool apply_overdispersion_change(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
    bool apply_sigma_change(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
    Tree *comparison(int m, double gamma, unsigned move_id, const vector<int> &cluster_sizes);
    void infer_mcmc(const vector<vector<double>> &D, const vector<int> &r, const vector<float> &move_probs,
                    int n_iters, unsigned int size_limit, double alpha, double gamma, const vector<int> &cluster_sizes);

    void update_t_scores();
    void random_initialize(u_int n_nodes, u_int n_regions, int max_iters, int max_regions_per_node=1); // randomly initializes a tree and copies it into the other
    void initialize_worked_example(int mode=1); // initializes the trees based on the test example
    void initialize_from_file(string path);
    vector<vector<int>> assign_cells_to_nodes(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes);
private:
    int deleted_node_idx();
};




void Inference::random_initialize(u_int n_nodes, u_int n_regions, int max_iters, int max_regions_per_node) {

    Tree *random_tree;
    int i = 0;

    if (n_nodes == 0)
        return;

    while(true)
    {
        i++;
        random_tree = new Tree(ploidy, n_regions, this->region_neutral_states);
        for (unsigned j = 0; j < n_nodes; ++j)
        {
            /*
             * Create a c_change hashmap using poisson and bernoulli
             *
             * */
            // create the map
            // create a map, fill it properly with r amount of labels
            map<u_int, int> distinct_regions;
            try {
                Utils::random_initialize_labels_map(distinct_regions, n_regions, max_regions_per_node, lambda_r); // modifies the distinct_regions
            }catch (const std::out_of_range& e)
            {
                if (verbosity > 2)
                    std::cout << " an out of range error was caught during the initialize labels map method, with message '"
                              << e.what() << "'\n";
                delete random_tree; // delete the tree
                random_tree = new Tree(ploidy, n_regions, this->region_neutral_states);
                break;
            }
            random_tree->random_insert(static_cast<map<u_int, int> &&>(distinct_regions));
        }

        if (i > max_iters)
        {
            throw runtime_error("a valid tree cannot be found after " + to_string(max_iters)  + " iterations. Please re-set the lambda_r, lambda_c and n_nodes variables.");
        }

        if (random_tree->get_n_nodes() == 0)
            continue;

        bool is_valid_tree = random_tree->is_valid_subtree(random_tree->root);
        bool is_redundant = random_tree->is_redundant();
        if (!is_valid_tree || is_redundant)
            delete random_tree;
        else
            break;
    }

    t = *random_tree;
    delete random_tree; // deallocate
    t.compute_weights();



}

void Inference::initialize_worked_example(int mode) {

    // build tree
    // tree that generated the data
    if (mode == 1) {
      t.random_insert({{0, 1}, {1, 1}});
      t.insert_at(1,{{1, 1}, {2, 1}});
      t.insert_at(2,{{0, -1}});
      t.insert_at(2,{{3, -1}});
      t.insert_at(1,{{0, 1}});
    } else {
      t.random_insert({{0, 1}}); // 1
      t.insert_at(1,{{1, 1}, {2, 1}, {3, -2}, {4, -2}}); // 2
      t.insert_at(1,{{1, 1}, {2, -1}, {3, -1}}); // 3
    }

    // Tree score: -2605.9655
//    t.insert_at(0,{{0,1},{1,1}}); // 1
//    t.insert_at(1,{{0,-1},{2,1}}); // 2
//    t.insert_at(2,{{3,-1}}); // 3

    t.compute_weights();

}

Inference::Inference(u_int n_regions, int n_cells, vector<int> &region_neutral_states, int ploidy, int verbosity, bool max_scoring, bool smoothed) : t(ploidy, n_regions, region_neutral_states),
                                                                   t_prime(ploidy, n_regions, region_neutral_states),
                                                                   best_tree(ploidy, n_regions, region_neutral_states),
                                                                   t_prime_sums(n_cells),
                                                                   t_sums(n_cells),
                                                                   t_maxs(n_cells),
                                                                   t_prime_maxs(n_cells)
                                                                   {

    this->n_regions = n_regions;
    this->ploidy = ploidy;
    this->verbosity = verbosity;
    this->max_scoring = max_scoring;
    this->n_cells = n_cells;
    this->region_neutral_states = region_neutral_states;
    this->smoothed = smoothed;
}

Inference::~Inference() {
    destroy();
}

bool Inference::apply_prune_reattach(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                                     bool validation_test_mode) {
    /*
     * Applies prune and reattach to t_prime
     * Updates the sums and scores tables partially
     * */

    Node* attached_node;
    attached_node = t_prime.prune_reattach(weighted, validation_test_mode);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);
        return true;
    }
    else
        return false;
}

void Inference::compute_t_table(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {
    // If D is cluster by region, we iterate over each cluster and compute the sum of each cluster's
    // score weighted by its size
    this->t_sums.clear();

    double currentMax = -DBL_MAX;
    std::pair<int, double> currentMax_node;

    int j = static_cast<int>(D.size());
    double t_sum = 0.0;
    for (int i = 0; i < j; ++i)
    {

        this->t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = this->t.root->get_children_id_score();

        this->t_scores.push_back(scores_vec);

        double aux = 0.0;
        if (max_scoring) {
          currentMax = -DBL_MAX;
          for (auto it = scores_vec.cbegin(); it != scores_vec.cend(); ++it ) {
              if (it->second > currentMax) {
                  currentMax_node.first = it->first;
                  currentMax_node.second = it->second;
                  currentMax = currentMax_node.second;
              }
          }
          aux = currentMax;
          this->t_maxs[i].first = currentMax_node.first;
          this->t_maxs[i].second = currentMax;
        }
        // if (max_scoring) {
        //   // Get maximum score for cell
        //   double currentMax = -DBL_MAX;
        //   unsigned arg_max = 0;
        //   for (auto it = scores_vec.cbegin(); it != scores_vec.cend(); ++it ) {
        //       if (it->second > currentMax)
        //           currentMax = it->second;
        //   }
        //   t_sum = currentMax;
        // }
        else {
            aux = MathOp::log_sum(scores_vec);
          }
        t_sum = t_sum + aux * cluster_sizes[i];
        this->t_sums.push_back(aux);
    }

    int m = D.size();
    int n = t.get_n_nodes();
    t.total_attachment_score = t_sum;
    t.prior_score = log_tree_prior(m, n);
    t.posterior_score = log_tree_posterior(t_sum, m, t);

}

void Inference::destroy() {
    // nothing to deallocate
}

Tree * Inference::comparison(int m, double gamma, unsigned move_id, const vector<int> &cluster_sizes) {
    /*
     * Returns the pointer to the accepted tree
     * m is size(D), i.e. number of cells
     * Throws std::out_of_range exception
     * */

    double t_prime_sum = 0.;
    for (int i=0; i<m;i++)
      t_prime_sum = t_prime_sum + t_prime_sums[i]*cluster_sizes[i];

    u_int t_n_nodes = t.get_n_nodes();
    u_int t_prime_n_nodes = t_prime.get_n_nodes();

    // compare the posteriors
    t_prime.posterior_score = log_tree_posterior(t_prime_sum, m, t_prime);
    // update the tree prior as well
    t_prime.prior_score = log_tree_prior(m, t_prime_n_nodes);
    // update the total attachment score
    t_prime.total_attachment_score = t_prime_sum;

    // check if move is weighted
    bool weighted;
    if ((move_id %2 == 1) && (move_id <=9)) // weighted move ids
        weighted = true;
    else
        weighted = false;

    // acceptance probability computations
    double score_diff = 0.0;
    if (move_id == 13) // dispersion change
    {
        double root_score_diff = t_prime.root_score - t.root_score;
        double posterior_score_diff = t_prime.posterior_score - t.posterior_score;

        score_diff = root_score_diff + posterior_score_diff;
    }
    else
        score_diff = t_prime.posterior_score - t.posterior_score;

    double log_acceptance_prob = 0.0; // later gets modified


    // compute nbd correction if posterior sampling
    double total_nbd_corr = 0.0;
    if (not max_scoring)
    {
        total_nbd_corr = 1.0;
        double nbd_corr= 1.0;
        double sum_chi=0.0, sum_chi_prime=0.0, sum_omega=0.0, sum_omega_prime=0.0;

        if (move_id == 1) // weighted prune-reattach
        {
            nbd_corr = t.cost() / t_prime.cost();
            assert(!std::isinf(nbd_corr));
            if (std::isinf(nbd_corr))
                return &t; // reject

            // ro variable
            total_nbd_corr *= nbd_corr;
        }
        else if (move_id == 6 || move_id == 7) // insert/delete move or weighted insert/delete move
        {

            vector<double> chi = t.chi_insert_delete(weighted);
            sum_chi = std::accumulate(chi.begin(), chi.end(), 0.0);
            vector<double> chi_prime = t_prime.chi_insert_delete(weighted);
            sum_chi_prime = std::accumulate(chi_prime.begin(), chi_prime.end(), 0.0);

            if (std::isinf(sum_chi))
                throw std::out_of_range("sum_chi is infinity, insert/delete move will be rejected");

            vector<double> omega = t.omega_insert_delete(lambda_r, lambda_c, weighted, max_scoring);
            sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
            vector<double> omega_prime = t_prime.omega_insert_delete(lambda_r, lambda_c, weighted, max_scoring);
            sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
        }
        else if (move_id == 8 || move_id == 9) // condense/split move or weighted cs
        {

            vector<double> chi = t.chi_condense_split(weighted);
            sum_chi = std::accumulate(chi.begin(), chi.end(), 0.0);
            vector<double> chi_prime = t_prime.chi_condense_split(weighted);
            sum_chi_prime = std::accumulate(chi_prime.begin(), chi_prime.end(), 0.0);

            vector<double> omega = t.omega_condense_split(lambda_s, weighted, max_scoring);
            sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
            vector<double> omega_prime = t_prime.omega_condense_split(lambda_s, weighted, max_scoring);
            sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
        }

        if (move_id == 6 || move_id == 7 || move_id == 8 || move_id == 9) // moves that require nbd correction
        {
            double n = static_cast<double>(t_n_nodes);
            if (t_n_nodes < t_prime_n_nodes) // insert, split
            {
                double weight = sum_chi/sum_omega_prime;
                total_nbd_corr *= weight;
            }
            else // delete, condense
            {
                double weight = sum_omega/sum_chi_prime;
                total_nbd_corr *= weight;
            }
        }

        if (move_id == 7 || move_id == 9) // weighted insert-delete or weighted condense-split
        {
            if (t_n_nodes > t_prime_n_nodes) // delete
            {
                // find the node that is deleted
                // use the get_id_score function since it returns a map having node id as a key
                map<int,double> t_prime_scores = t_prime.root->get_children_id_score();
                Node* deleted = nullptr;
                for (auto const &t_node : t.root->get_descendents(false)) // root won't be contained!
                    if (!t_prime_scores.count(t_node->id))
                    {
                        deleted = t_node;
                        break;
                    }
                if (deleted == nullptr)
                    throw std::logic_error("The deleted node could not be found in t!");

                unsigned d_i_T = deleted->n_descendents;
                unsigned d_i_T_prime = 0;

                // find it's parent in t_prime
                int parent_id = deleted->parent->id;
                for (auto const &t_prime_node : t_prime.root->get_descendents(true))
                    if (parent_id == t_prime_node->id)
                    {
                        d_i_T_prime = t_prime_node->n_descendents;
                        break;
                    }
                if (d_i_T == 0)
                    throw std::logic_error("The deleted node's parent could not be found in t_prime!");

                double p_add_corr_num = d_i_T_prime + 1;
                double p_add_corr_denom = 2 * d_i_T;
                double ratio = p_add_corr_num / p_add_corr_denom;
                total_nbd_corr *= ratio;
            }
            else // insert
            {
                // find the node that is inserted in t_prime
                map<int,double> t_scores = t.root->get_children_id_score();

                Node* added = nullptr;
                for (auto const &t_prime_node : t_prime.root->get_descendents(false)) // without root
                    if (!t_scores.count(t_prime_node->id))
                    {
                        added = t_prime_node;
                        break;
                    }
                if (added == nullptr)
                    throw std::logic_error("The inserted node could not be found in t_prime!");

                unsigned d_i_T_prime = added->n_descendents;
                unsigned d_i_T = 0;
                // find it's parent in t
                int parent_id = added->parent->id;
                for (auto const &t_node : t.root->get_descendents(true))
                    if (parent_id == t_node->id)
                    {
                        d_i_T = t_node->n_descendents;
                        break;
                    }
                if (d_i_T == 0)
                    throw std::logic_error("The inserted node's parent could not be found in t!");

                double p_add_corr_num = 2 * d_i_T_prime;
                double p_add_corr_denom = d_i_T + 1;
                double ratio = p_add_corr_num / p_add_corr_denom;
                total_nbd_corr *= ratio;
            }
        }
    }

    double log_total_nbd_corr;
    if (total_nbd_corr == 0.0)
        log_total_nbd_corr = 0.0;
    else
        log_total_nbd_corr = std::log(total_nbd_corr);

    log_acceptance_prob = log_total_nbd_corr + gamma * score_diff;

    // reject the move if nan or inf occurs
    if (std::isnan(log_acceptance_prob) || std::isinf(log_acceptance_prob))
        return &t;

    if (verbosity > 2)
        cout << "log acceptance prob: " << log_acceptance_prob << endl;

    if (log_acceptance_prob > 0)
    {
        if (verbosity > 2)
            std::cout << "Move accepted." << std::endl;
        return &t_prime;
    }

    else
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        boost::random::uniform_real_distribution<double> distribution(0.0,1.0);
        double rand_val = distribution(gen);
        rand_val = std::log(rand_val); // take the log

        if (verbosity > 2)
            cout<<"log rand_val: "<<rand_val<<endl;

        if (log_acceptance_prob > rand_val)
        {
            if (verbosity > 2)
                std::cout << "Move accepted." << std::endl;
            return &t_prime;
        }
        else
        {
            if (verbosity > 2)
                std::cout << "Move rejected." << std::endl;
            return &t;
        }

    }
}

void Inference::infer_mcmc(const vector<vector<double>> &D, const vector<int> &r, const vector<float> &move_probs,
                           int n_iters, unsigned int size_limit, double alpha, double gamma, const vector<int> &cluster_sizes) {

    int m = static_cast<int>(D.size());

    u_int n_rejected = 0;
    u_int n_accepted = 0;

    // for writing the posteriors on file
    std::ofstream mcmc_scores_file;
    std::ofstream rel_mcmc_scores_file;
    std::ofstream acceptance_ratio_file;
    std::ofstream gamma_file;

    if (verbosity > 1)
    {
        mcmc_scores_file.open(f_name_posfix + "_markov_chain.csv", std::ios_base::out);
        rel_mcmc_scores_file.open(f_name_posfix + "_rel_markov_chain.csv", std::ios_base::out);
        acceptance_ratio_file.open(f_name_posfix + "_acceptance_ratio.csv", std::ios_base::out);
        gamma_file.open(f_name_posfix + "_gamma_values.csv", std::ios_base::out);
    }

    best_tree = t; //start with the t

    using namespace std::placeholders;
    unsigned n_apply_move = 50; // number of times move is applied
    for (int i = 0; i < n_iters; ++i) {

        if (verbosity > 2)
            std::cout << "MCMC iteration " << i << ":" <<std::endl;


        if (verbosity > 1 && (i % 10000 == 0))
        {
            std::cout << "iteration" << i <<  "tree score" << t.posterior_score + t.root_score  << std::endl;
            if (smoothed)
              std::cout << "sigma: " << t.sigma << "root score: " << t.root_score << std::endl;
            else
              std::cout << "nu: " << t.nu << "root score: " << t.root_score << std::endl;
        }


        bool rejected_before_comparison = false;

        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        boost::random::discrete_distribution<> d(move_probs.begin(), move_probs.end());

        unsigned move_id = d(gen);

        switch (move_id)
        {
            case 0:
            {
                // prune & reattach
                if (verbosity > 2)
                    cout << "Prune and reattach" << endl;

                auto func = std::bind(&Inference::apply_prune_reattach, this, _1, _2, false, false);
                bool prune_reattach_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not prune_reattach_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Prune and reattach rejected before comparison"<<endl;
                }
                break;
            }
            case 1:
            {
                // weighted prune & reattach
                if (verbosity > 2)
                    cout<<"Weighted prune and reattach"<<endl;

                auto func = std::bind(&Inference::apply_prune_reattach, this, _1, _2, true, false);
                bool weighted_prune_reattach_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not weighted_prune_reattach_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Weighted prune and reattach rejected before comparison"<<endl;
                }
                break;
            }
            case 2:
            {
                // swap labels
                if (verbosity > 2)
                    cout << "swap labels" << endl;

                auto func = std::bind(&Inference::apply_swap, this, _1, _2, false, false);
                bool swap_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not swap_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Swap labels rejected before comparison" << endl;
                }
                break;
            }
            case 3:
            {
                // weighted swap labels
                if (verbosity > 2)
                    cout << "weighted swap labels" << endl;

                auto func = std::bind(&Inference::apply_swap, this, _1, _2, true, false); // weighted: true
                bool weighted_swap_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not weighted_swap_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Weighted swap labels rejected before comparison" << endl;
                }
                break;
            }
            case 4:
            {
                // add or remove event
                if (verbosity > 2)
                    cout << "add or remove event" << endl;

                auto func = std::bind(&Inference::apply_add_remove_events, this, _1, _2, false, false);
                bool add_remove_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not add_remove_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Add or remove event rejected before comparison"<<endl;
                }
                break;
            }
            case 5:
            {
                // weighted add or remove event
                if (verbosity > 2)
                    cout << "weighted add or remove event" << endl;

                auto func = std::bind(&Inference::apply_add_remove_events, this, _1, _2, true, false); // weighted=true
                bool add_remove_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not add_remove_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Weighted add or remove event rejected before comparison"<<endl;
                }
                break;
            }
            case 6:
            {
                // insert delete node
                if (verbosity > 2)
                    cout << "insert/delete node" << endl;

                auto func = std::bind(&Inference::apply_insert_delete_node, this, _1, _2, _3, false); // weighted=false
                bool insert_delete_success = apply_multiple_times(n_apply_move, func, D, r, size_limit);

                if (not insert_delete_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "insert/delete rejected before comparison" << endl;
                }
                break;
            }
            case 7:
            {
                // weighted insert delete node
                if (verbosity > 2)
                    cout << "weighted insert/delete node" << endl;

                auto func = std::bind(&Inference::apply_insert_delete_node, this, _1, _2, _3, true); // weighted=true
                bool insert_delete_success = apply_multiple_times(n_apply_move, func, D, r, size_limit);

                if (not insert_delete_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "weighted insert/delete rejected before comparison" << endl;
                }
                break;
            }
            case 8:
            {
                // condense split move
                if (verbosity > 2)
                    cout << "condense split move " <<endl;

                auto func = std::bind(&Inference::apply_condense_split, this, _1, _2, _3, false); // weighted=false
                bool condense_split_success = apply_multiple_times(n_apply_move, func, D, r, size_limit);

                if (not condense_split_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "condense/split move rejected before comparison"<<endl;
                }
                break;
            }
            case 9:
            {
                // weighted condense split move
                if (verbosity > 2)
                    cout << "weighted condense split move " <<endl;

                auto func = std::bind(&Inference::apply_condense_split, this, _1, _2, _3, true); // weighted=true
                bool condense_split_success = apply_multiple_times(n_apply_move, func, D, r, size_limit);

                if (not condense_split_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "weighted condense/split move rejected before comparison"<<endl;
                }
                break;
            }
            case 10:
            {
                // genotype_preserving prune & reattach
                if (verbosity > 2)
                    cout<<"Genotype preserving prune and reattach"<<endl;
                bool genotype_prune_reattach_success = apply_genotype_preserving_pr(gamma);
                if (not genotype_prune_reattach_success)
                {
                    n_rejected++;
                    acceptance_ratio_file << std::setprecision(print_precision) << -1 << ((i==n_iters-1) ? "" : ",");
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Genotype preserving prune/reattach rejected"<<endl;
                }
                else // accepted
                {
                    n_accepted++;
                    acceptance_ratio_file << std::setprecision(print_precision) << static_cast<int>(move_id) << ((i==n_iters-1) ? "" : ",");
                    // update t score
                    double t_sum = 0.0;
                    for (int i=0; i<m;i++)
                      t_sum = t_sum + t_sums[i]*cluster_sizes[i];
                    t.posterior_score = log_tree_posterior(t_sum, m, t); // the prior score will change
                    t_prime = t; // update t_prime
                }
                break;
            }
            case 11:
            {
              // expand or shrink block
              if (verbosity > 2)
                  cout << "expand or shrink block" << endl;

              auto func = std::bind(&Inference::apply_expand_shrink_blocks, this, _1, _2, false, false);
              bool add_remove_success = apply_multiple_times(n_apply_move, func, D, r);

              if (not add_remove_success) {
                  rejected_before_comparison = true;
                  if (verbosity > 2)
                      cout << "Expand or shrink block rejected before comparison"<<endl;
              }
              break;
            }
            case 12:
            {
                // add common ancestor move
                if (verbosity > 2)
                    cout << "Add common ancestor" << endl;

                auto func = std::bind(&Inference::apply_add_common_ancestor, this, _1, _2, false);
                bool add_common_ancestor_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not add_common_ancestor_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Add common ancestor move rejected before comparison" << endl;
                }
                break;
            }
            case 13:
            {
                // changing dispersion (overdispersion or standard deviation) move
                if (verbosity > 2)
                    cout<<"Update dispersion move"<<endl;

                auto func = std::bind(&Inference::apply_overdispersion_change, this, _1, _2, _3);
                if (smoothed)
                    auto func = std::bind(&Inference::apply_sigma_change, this, _1, _2, _3);
                bool disp_change_success = apply_multiple_times(n_apply_move, func, D, r, cluster_sizes);

                if (not disp_change_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Update dispersion move rejected before comparison"<<endl;
                }
                break;
            }
            case 14:
            {
                // delete leaf move
                if (verbosity > 2)
                    cout << "Delete leaf move" << endl;

                auto func = std::bind(&Inference::apply_delete_leaf, this, _1, _2);
                bool delete_leaf_success = apply_multiple_times(n_apply_move, func, D, r);

                if (not delete_leaf_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 2)
                        cout << "Delete leaf move rejected before comparison" << endl;
                }
                break;
            }
            default:
                throw std::logic_error("undefined move index");
        } // end of switch

        Tree* accepted;

        if (move_id == 10) //genotype preserving prune reattach, gibbs sampling
            accepted = &t;
        else // metropolis hasting sampling
        {
            // compare the trees
            if (rejected_before_comparison)
                accepted = &t;
            else
                try
                {
                    accepted = comparison(m, gamma, move_id, cluster_sizes);
                } catch (const std::out_of_range& e)
                {
                    if (verbosity > 2)
                    {
                        std::cout << " an out of range error was caught during the comparison function, with message '"
                                  << e.what() << '\'' << std::endl;
                        accepted = &t;
                    }
                }
                catch (const std::exception &e)
                {
                    if (verbosity > 2)
                    {
                        std::cout << " an exception was caught during the comparison function, with message '"
                                  << e.what() << '\'' << std::endl;
                        accepted = &t;
                    }
                }

            double score_diff = t_prime.posterior_score - t.posterior_score;
            // update trees and the matrices
            if (accepted == &t_prime)
            {
                if (move_id != 13 && std::abs(score_diff) > 0.1)
                    gamma *= exp((1.0-alpha)*alpha);
                acceptance_ratio_file << std::setprecision(print_precision) << static_cast<int>(move_id) << ((i==n_iters-1) ? "" : ",");
                n_accepted++;
                t_sums = t_prime_sums;
                t_maxs = t_prime_maxs;
                update_t_scores(); // this should be called before t=tprime, because it checks the tree sizes in both.
                t = t_prime;
                if ((t_prime.posterior_score + t_prime.root_score)  > (best_tree.posterior_score + best_tree.root_score))
                    best_tree = t_prime;
            }
            else
            {
                // print acceptance ratio
                if (move_id != 13 && std::abs(score_diff) > 0.1)
                    gamma *= exp((0.0-alpha)*alpha);
                acceptance_ratio_file << std::setprecision(print_precision) << -1 << ((i==n_iters-1) ? "" : ",");
                n_rejected++;
                t_prime = t;
            }

            if (gamma > 100.0)
                gamma = 100.0;
            if (gamma < 1.0/100.0)
                gamma = 1.0/100.0;

            // t_prime_sums.clear();
            t_prime_scores.clear();
        }

        if (verbosity > 1)
        {
            static double first_score = accepted->posterior_score + accepted->root_score; // the first value will be kept in whole program
            // print accepted log_posterior

            if (i == n_iters - 1)
            {
                mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->root_score;
                rel_mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->root_score - first_score;
                gamma_file << gamma;
            }
            else
            {
                mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->root_score << ',';
                rel_mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->root_score - first_score << ',';
                gamma_file << gamma << ',';
            }


        }
    }
}

bool Inference::apply_overdispersion_change(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {
    /*
     * Changes the current overdispersion parameter nu using gaussian random walk.
     * */

    try
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        double rand_val = 0.0;
        boost::random::normal_distribution<double> distribution(0.0,0.02);
        rand_val = distribution(gen);

        double log_t_prime_nu = std::log(t_prime.nu) + rand_val;

        if (log_t_prime_nu > 10.0)
            return false; // reject the move

        if (verbosity > 2)
          std::cout<<"Old nu value: " << t_prime.nu << ",\t";

        t_prime.nu = std::exp(log_t_prime_nu); // real space

        if (verbosity > 2)
          std::cout<<"new nu value: " << t_prime.nu << std::endl;


        this->compute_t_prime_root_scores(D,r, cluster_sizes);
        compute_t_prime_scores(t_prime.root, D, r);
        compute_t_prime_sums(D);
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 2)
            std::cout << " a standard exception was caught during the apply overdispersion change move, with message '"
                      << e.what() << "'\n";
        return false;
    }

    return true;
}


bool Inference::apply_sigma_change(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {
    /*
     * Changes the current standard deviation parameter sigma using gaussian random walk.
     * */

    try
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        double rand_val = 0.0;
        boost::random::normal_distribution<double> distribution(0.0,0.02);
        rand_val = distribution(gen);

        double log_t_prime_sigma = std::log(t_prime.sigma) + rand_val;

        if (log_t_prime_sigma > 10.0)
            return false; // reject the move

        if (verbosity > 2)
          std::cout<<"Old sigma value: " << t_prime.sigma << ",\t";

        t_prime.sigma = std::exp(log_t_prime_sigma); // real space

        if (verbosity > 2)
          std::cout<<"new sigma value: " << t_prime.sigma << std::endl;


        this->compute_t_prime_root_scores(D,r, cluster_sizes);
        compute_t_prime_scores(t_prime.root, D, r);
        compute_t_prime_sums(D);
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 2)
            std::cout << " a standard exception was caught during the apply sigma change move, with message '"
                      << e.what() << "'\n";
        return false;
    }

    return true;
}

double Inference::log_tree_prior(int m, int n) {
    /*
     * Computes and returns the tree prior.
     * m: n_cells
     * n: n_nodes
     * tree_sum: sum of all scores across all cells and attachment points
     * */

//    double log_prior = - (n -1 + m) * log(n+1) -m * n * log(2); // tree prior
    double combinatorial_penalization = 0.0;
    combinatorial_penalization = - cf*m*n*log(2);

    double tree_prior = 0.0;
    if (max_scoring) {
      m = 0.;
    }
    tree_prior = -(n-1+m)*log(n+1);

    double log_prior = tree_prior + combinatorial_penalization;

    return log_prior;
}

double Inference::log_tree_posterior(double tree_sum, int m, Tree &tree) {
    /*
     * Computes and returns the posterior log score of the tree
     * m: n_cells
     * */

    // n: n_nodes
    int n = tree.get_n_nodes();
    double log_posterior = 0.0;
    log_posterior = tree_sum + this->log_tree_prior(m, n); // initialise posterior with prior then add posterior

    double PV = tree.event_prior();
    log_posterior += PV;

    return log_posterior;
}

void Inference::update_t_scores() {

    // iterate over t_prime_scores
    // if index exists t_scores, update, else insert

    int deleted_index = deleted_node_idx(); // if -1 then not deleted, otherwise the index of the deleted


    for (unsigned k=0; k < t_prime_scores.size(); k++) // iterates over n_cells
        for (auto const& x : t_prime_scores[k])
        {
            t_scores[k][x.first] = t_prime_scores[k][x.first]; // if found in t_scores[k] map, then updates. Else inserts.
            if (deleted_index != -1)
                t_scores[k].erase(deleted_index);
        }

}


void Inference::compute_t_prime_scores(Node *attached_node, const vector<vector<double>> &D, const vector<int> &r) {

    // if the t_prime_scores is empty then fill it with size(D) elements
    // otherwise append the results to the first size(D) places in t_prime_scores
    // size of t_prime_scores should always be size(D)

    bool is_empty_table = t_prime_scores.empty();

    int j = 0;

    // I was doing this to make sure t_prime_scores would have every node when I compute the maximum, which makes this slower
    // if (max_scoring)
    //   attached_node = t_prime.root;

    for (auto const &d: D)
    {
        double sum_d = accumulate( d.begin(), d.end(), 0.0);

        if (attached_node != t_prime.root)
            attached_node->parent->attachment_score = t_scores[j][attached_node->parent->id]; // the indices must match
        // attached node->parent->id must match the all_nodes_vec index
        t_prime.compute_stack(attached_node, d, sum_d,r);

        if (is_empty_table)
            t_prime_scores.push_back(attached_node->get_children_id_score());
        else
        {
            // append the contents of second hashmap into the first
            // note that they cannot have overlapping keys
            for (auto const& map_item : attached_node->get_children_id_score()) // go to each child of the changed node and get its score
            {
                t_prime_scores[j][map_item.first] = map_item.second;
            }
        }

        j++;
    }
}

bool Inference::apply_swap(const vector<vector<double>> &D, const vector<int> &r, bool weighted, bool test_mode) {

    vector<Node*> swapped_nodes;
    swapped_nodes = t_prime.swap_labels(weighted, test_mode);

    for (auto const &node : swapped_nodes)
    {
        compute_t_prime_scores(node, D, r);
    }
    compute_t_prime_sums(D);

    return true;

}

void Inference::compute_t_prime_sums(const vector<vector<double>> &D) {

    /*
     * Computes the t_prime sums that represent the partial computed sub-tree
     * Takes the structural changes and tree size changes into account.
     * */


    // In case of a delete node the removed node is also added to the old_vals
    // in delete case: add all t_scores to old_vals and remove the ones not found in t_prime_scores
    // TODO: same lines of code is also used in update_t_scores method, make it reusable

    int deleted_index = deleted_node_idx(); // if -1 then not deleted, otherwise the index of the deleted

    if (max_scoring) {
      // Parallelize cells
      #pragma omp parallel for
        for (unsigned i = 0; i < D.size(); i++) { // each cell i
            vector<double> old_vals; // old values of kept nodes
            old_vals.reserve(t_scores[i].size()); // the max possible size

            map<int,double> unchanged_vals = t_scores[i]; // initiate with all values and remove the changed ones

            // Instead of going through all the nodes, go through the updated ones to get their best
            // and compare with previous best
            double currentMax_newval = -DBL_MAX;
            pair<int, double> currentMax_newnode;
            bool prev_max_in_new = false;
            for (auto &u_map : t_prime_scores[i]) { // go through each updated node in T'
                if (t_scores[i].count(u_map.first)) // if this node was already present
                {
                    unchanged_vals.erase(u_map.first); // erase it from the unchanged vals
                }

                // Check if an updated node was the previous maximum
                if (t_maxs[i].first == u_map.first) {
                  prev_max_in_new = true;
                }

                // Best new node, including the inserted ones
                if (u_map.second > currentMax_newval) {
                  currentMax_newval = u_map.second;
                  currentMax_newnode.first = u_map.first;
                  currentMax_newnode.second = currentMax_newval;
                }
            }

            if (deleted_index != -1)
            {
                unchanged_vals.erase(deleted_index); // erase it from the unchanged vals
            }

            pair<int, double> prev_best;
            prev_best.first = t_maxs[i].first; // If the previous max was not changed, take it
            prev_best.second = t_maxs[i].second;
            // If the previous max was changed (or if it was deleted), we need find the maximum among the unchanged nodes.
            if (prev_max_in_new || (t_maxs[i].first == deleted_index)) { // find best unchanged and compare against it
              double currentMax = -DBL_MAX;
              pair<int, double> currentMax_unchangednode;
              for (auto &u_map : unchanged_vals) {
                  if (u_map.second > currentMax) {
                      currentMax_unchangednode.first = u_map.first;
                      currentMax_unchangednode.second = u_map.second;
                      currentMax = u_map.second;
                  }
              }
              prev_best.first = currentMax_unchangednode.first;
              prev_best.second = currentMax_unchangednode.second;
            }

            // Update cell score
            pair<int, double> newMax_node;
            newMax_node.first = currentMax_newnode.first;
            newMax_node.second = currentMax_newval;
            if (prev_best.second > currentMax_newval) {
              newMax_node.first = prev_best.first;
              newMax_node.second = prev_best.second;
            }
            t_prime_maxs[i].first = newMax_node.first;
            t_prime_maxs[i].second = newMax_node.second;

            assert(!std::isnan(newMax_node.second));
            t_prime_sums[i] = newMax_node.second;
        }
    }
    else {
      // Parallelize cells
      #pragma omp parallel for
        for (unsigned i = 0; i < D.size(); ++i) {
            vector<double> old_vals; // old values of kept nodes
            old_vals.reserve(t_scores[i].size()); // the max possible size
            vector<double> new_vals; // new values of kept nodes
            new_vals.reserve(t_scores[i].size());

            map<int,double> unchanged_vals = t_scores[i]; // initiate with all values and remove the changed ones


            for (auto &u_map : t_prime_scores[i]) {
                if (t_scores[i].count(u_map.first)) // add only if it is existing in the old vals // for the insertion case
                {
                    old_vals.push_back(t_scores[i][u_map.first]); // again the indices should match
                    unchanged_vals.erase(u_map.first); // erase it from the unchanged vals
                }


                new_vals.push_back(u_map.second);
            }

            if (deleted_index != -1)
            {
                old_vals.push_back(t_scores[i][deleted_index]);
                unchanged_vals.erase(deleted_index); // erase it from the unchanged vals
            }


            // in case of delete, the deleted val is in old_vals not in unchanged

            double res = MathOp::log_replace_sum(t_sums[i], old_vals, new_vals, unchanged_vals); // it takes t_sums[i]
            // subtracts the olds and adds the news
            // in case of delete, subtract an extra value
            // in case of insert, add an extra value
            // if the tree size changes, update it (tree.n_nodes). Posterior takes that into account
            assert(!std::isnan(res));

            t_prime_sums[i] = res;
        }
    }
}

bool Inference::apply_add_remove_events(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                                        bool validation_test_mode)
{
    /*
     * Applies add/remove event to t_prime
     * Updates the sums and scores tables partially
     * */

    // weighted = false
    Node* attached_node;

    attached_node = t_prime.add_remove_events(weighted, validation_test_mode);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}


bool Inference::apply_expand_shrink_blocks(const vector<vector<double>> &D, const vector<int> &r, bool weighted,
                                        bool validation_test_mode)
{
    /*
     * Applies expand/shrink blocks to t_prime
     * Updates the sums and scores tables partially
     * */

    // weighted = false
    Node* attached_node;

    attached_node = t_prime.expand_shrink_blocks(weighted, validation_test_mode);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}

bool Inference::apply_insert_delete_node(const vector<vector<double>> &D, const vector<int> &r, unsigned int size_limit,
                                         bool weighted) {
    /*
     * Applies the insert/delete move on t_prime
     * Updates the sums and scores tables partially
     * */

    Node* tobe_computed;

    tobe_computed = t_prime.insert_delete_node(size_limit, weighted, max_scoring);


    if (tobe_computed != nullptr)
    {
        compute_t_prime_scores(tobe_computed, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}

int Inference::deleted_node_idx() {
    /*
     * Finds and returns the index of the deleted node, should a node be deleted.
     * Return -1 if not found.
     * */

    int deleted_index = -1;
    if (t_prime.all_nodes_vec.size() < t.all_nodes_vec.size()) // a node is deleted
    {
        // find the index of the deleted
        set<int> set_tprime_nodes;
        for (auto item : t_prime.all_nodes_vec)
            set_tprime_nodes.insert(item->id);
        for (auto &item : t.all_nodes_vec)
            if (set_tprime_nodes.count(item->id) == 0)
                deleted_index = item->id;
    }

    return deleted_index;
}

bool Inference::apply_condense_split(const vector<vector<double>> &D, const vector<int> &r, unsigned int size_limit,
                                     bool weighted) {
    /*
     * Applies the condense/delete move on t_prime
     * Updates the sums and scores tables partially
     * */

    Node* tobe_computed;
    tobe_computed = t_prime.condense_split_node(size_limit, weighted, max_scoring);

    if (tobe_computed != nullptr)
    {
        compute_t_prime_scores(tobe_computed, D, r);
        compute_t_prime_sums(D);
        return true;
    }
    else
        return false;
}

vector<vector<int>> Inference::assign_cells_to_nodes(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {


    /*
     * re-computes the best tree to assigns cells to nodes by the maximum score
     * Writes the cell-node_ids, region sizes, and the inferred CNVs files.
     * Returns the inferred CNVs matrix (with the ploidy added)
     * */

    t_scores.clear();
    t_sums.clear();
    t = best_tree;
    this->compute_t_table(D,r,cluster_sizes);

    std::ofstream cell_node_ids_file;
    std::ofstream cell_region_cnvs_file;
    if (verbosity > 1)
    {
        cell_node_ids_file.open(f_name_posfix + "_cell_node_ids.tsv");
        cell_region_cnvs_file.open(f_name_posfix + "_cell_region_cnvs.csv");
    }
    // create a hashmap of nodes for constant access by id
    unordered_map<uint64_t , Node*> hash_map;
    for (unsigned i=0; i < t.all_nodes_vec.size(); i++)
    {
        hash_map[t.all_nodes_vec[i]->id] = t.all_nodes_vec[i];
    }

    size_t n_cells = t_scores.size();

    vector<vector<int>> cell_regions(n_cells, vector<int>(this->n_regions, ploidy)); //fill vector, initialise with ploidy
    for (u_int i = 0; i < n_cells; ++i) {
      cell_regions[i] = this->region_neutral_states;
    }

    for (size_t j = 0; j < n_cells; ++j) {
        // t_scores[i] is the map
        pair<const int, double> max_pair = *max_element(t_scores[j].begin(), t_scores[j].end(), [] (const pair<const int, double>& p1, const pair<const int, double>& p2)
        {
            return p1.second < p2.second;
        }) ;
        if (verbosity > 1)
            cell_node_ids_file << j << '\t' << max_pair.first << '\n';

        Node* max_node = hash_map[max_pair.first];

        for (auto const& x : max_node->c) // iterate over map
        {
            cell_regions[j][x.first] = x.second + this->region_neutral_states[x.first];
        }

    }
    if (verbosity > 1)
    {
        for (size_t k = 0; k < n_cells; ++k) {
            for (u_int i = 0; i < n_regions; ++i) {
                if (i == n_regions-1) // the last element
                    cell_region_cnvs_file << cell_regions[k][i];
                else // add comma
                    cell_region_cnvs_file << cell_regions[k][i] << ',';
            }
            cell_region_cnvs_file << '\n';
        }
    }

    return cell_regions;
}


void Inference::initialize_from_file(string path) {
    /*
     * Initializes the tree t from the file
     * */
    t.load_from_file(path);


}

std::vector<double>
Inference::get_tree_root_scores(const vector<vector<double>> &D, const vector<int> &r, const Tree &tree) {
    /*
     * Returns the vector of root scores for every cell
     * */

    u_int n_cells = D.size();

    vector<double> scores(n_cells);
    for (u_int i = 0; i < n_cells; ++i) {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        scores[i] = tree.get_root_score(r, sum_d, D[i]);
    }

    return scores;
}

void Inference::compute_t_root_scores(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {
/*
 * Computes and stores the overdispersed root sum for tree t across all cells
 * */
    vector<double> t_od_root_scores = get_tree_root_scores(D,r,this->t);

    double sum_t_od_root_scores = 0.;
    int m = D.size();
    for (int i=0; i < m; i++)
      sum_t_od_root_scores = sum_t_od_root_scores + t_od_root_scores[i]*cluster_sizes[i];

    this->t.root_score = sum_t_od_root_scores;

}

void Inference::compute_t_prime_root_scores(const vector<vector<double>> &D, const vector<int> &r, const vector<int> &cluster_sizes) {
/*
 * Computes and stores the overdispersed root sum for tree t prime across all cells
 * */
    vector<double> t_prime_od_root_scores = get_tree_root_scores(D,r,this->t_prime);

    double sum_t_prime_od_root_scores = 0.;
    int m = D.size();
    for (int i=0; i < m; i++)
      sum_t_prime_od_root_scores = sum_t_prime_od_root_scores + t_prime_od_root_scores[i]*cluster_sizes[i];

    this->t_prime.root_score = sum_t_prime_od_root_scores;
}

void Inference::update_t_prime() {
    /*
     * Sets t_prime to t
     * */

    this->t_prime = this->t;
}

bool Inference::apply_genotype_preserving_pr(double gamma) {
    /*
     * Applies the genotype preserving prune and reattach move in a gibbs sample setting.
     * */

    try
    {
        this->t.genotype_preserving_prune_reattach(gamma);
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 2)
            std::cout << " a standard exception was caught during the genotype preserving prune and reattach move,"
                         " with message '" << e.what() << '\'' << std::endl;
        return false;
    }

    return true;
}

template<class... Ts, class AnyFunction>
bool Inference::apply_multiple_times(unsigned n, AnyFunction func, Ts &...args) {
    /*
     * Applies the boolean function func multiple times with the variable arguments
     * n: the number of tries
     * func: function to be called
     * args: the variable size of arguments
     * */

    for (unsigned i = 0; i < n; ++i) {

        bool is_successful = false;
        try {
            is_successful = func(args...);
        }
        catch (const InvalidMove& e) // don't try again if the move is invalid
        {
            if (verbosity > 2)
                std::cout << " an invalid move has occurred, with message '"
                          << e.what() << "'\n";
            break;
        }
        catch (const std::exception& e)
        {
            if (verbosity > 2)
                std::cout << " a standard exception was caught during apply_multiple_times, with message '"
                          << e.what() << "'\n";
            t_prime = t; // undo the changes made on t_prime
            continue; // try again
        }
        if (is_successful)
            return true;
        else // reset the changes on t_prime
            t_prime = t;
    }
    return false;
}

bool Inference::apply_delete_leaf(const vector<vector<double>> &D, const vector<int> &r) {
    /*
     * Applies the delete leaf move on t_prime.
     * */

    Node* tobe_computed;

    tobe_computed = t_prime.delete_leaf();

    compute_t_prime_scores(tobe_computed, D, r);
    compute_t_prime_sums(D);

    return true;

}

bool Inference::apply_add_common_ancestor(const vector<vector<double>> &D, const vector<int> &r, bool validation_test_mode) {
    /*
     * Applies the delete leaf move on t_prime.
     * */

    Node* tobe_computed;
    tobe_computed = t_prime.add_common_ancestor(validation_test_mode);

    compute_t_prime_scores(tobe_computed, D, r);
    compute_t_prime_sums(D);

    return true;

}

#endif //SC_DNA_INFERENCE_H
