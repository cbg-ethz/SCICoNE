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
#include <math.h>
#include <array>
#include "globals.cpp"
#include "Lgamma.h"

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

class Inference {
/*
 * Contains functionality to perform monte carlo markov chains (mcmc) inference
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
    std::vector<std::map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;
    std::string f_name;
    int verbosity;

public:
    Inference(u_int n_regions, int ploidy=2, int verbosity=2);
    ~Inference();
    void destroy();
    void compute_t_od_scores(const vector<vector<double>> &D, const vector<int> &r);
    void compute_t_prime_od_scores(const vector<vector<double>> &D, const vector<int> &r);
    std::vector<double>
    get_tree_od_root_scores(const vector<vector<double>> &D, const vector<int> &r, const Tree &tree);
    void compute_t_table(const vector<vector<double>> &D, const vector<int> &r);
    void compute_t_prime_scores(Node *attached_node, const vector<vector<double>> &D, const vector<int> &r);
    void compute_t_prime_sums(const vector<vector<double>> &D);
    double log_tree_prior(int m, int n);
    double log_posterior(double tree_sum, int m, Tree &tree);
    bool apply_prune_reattach(const vector<vector<double>> &D, const vector<int> &r, bool genotype_preserving,
                                  bool weighted, bool validation_test_mode);
    bool apply_add_remove_events(double lambda_r, double lambda_c, const vector<vector<double>> &D,
                                 const vector<int> &r,
                                 bool weighted = false,
                                 bool validation_test_mode = false);
    bool apply_insert_delete_node(double lambda_r, double lambda_c, const vector<vector<double>> &D,
                                      const vector<int> &r, unsigned int size_limit, bool weighted);
    bool apply_condense_split(double lambda_s, const vector<vector<double>> &D, const vector<int> &r,
                                  unsigned int size_limit, bool weighted);
    bool apply_swap(const vector<vector<double>> &D, const vector<int> &r, bool weighted = false,
                    bool test_mode = false);
    bool apply_overdispersion_change(const vector<vector<double>> &D, const vector<int> &r);
    Tree *comparison(int m, double gamma, unsigned move_id);
    void infer_mcmc(const vector<vector<double>> &D, const vector<int> &r, const vector<float> &move_probs, int n_iters,
                    unsigned int size_limit);
    void write_best_tree();
    void update_t_scores();
    void random_initialize(u_int n_nodes, u_int n_regions, int max_iters); // randomly initializes a tree and copies it into the other
    void initialize_worked_example(); // initializes the trees based on the test example
    void initialize_from_file(string path);
    vector<vector<int>> assign_cells_to_nodes(const vector<vector<double>> &D, const vector<int> &r);
private:
    int deleted_node_idx();
};




void Inference::random_initialize(u_int n_nodes, u_int n_regions, int max_iters) {

    Tree *random_tree;
    int i = 0;
    while(true)
    {
        i++;
        random_tree = new Tree(ploidy, n_regions);
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
                Utils::random_initialize_labels_map(distinct_regions, n_regions); // modifies the distinct_regions
            }catch (const std::out_of_range& e)
            {
                if (verbosity > 0)
                    std::cout << " an out of range error was caught during the initialize labels map method, with message '"
                              << e.what() << "'\n";
                delete random_tree; // delete the tree
                random_tree = new Tree(ploidy, n_regions);
                break;
            }
            random_tree->random_insert(static_cast<map<u_int, int> &&>(distinct_regions));
        }
        if (random_tree->get_n_nodes() == 0)
            continue;

        if (i > max_iters)
        {
            throw runtime_error("a valid tree cannot be found after " + to_string(max_iters)  + " iterations. Please re-set the lambda_r, lambda_c and n_nodes variables.");
        }

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

void Inference::initialize_worked_example() {

    // build tree
    // tree that generated the data
    t.random_insert({{0, 1}, {1, 1}});
    t.insert_at(1,{{1, 1}, {2, 1}});
    t.insert_at(2,{{0, -1}});
    t.insert_at(2,{{3, -1}});
    t.insert_at(1,{{1, 1}});

    // Tree score: -2605.9655
//    t.insert_at(0,{{0,1},{1,1}}); // 1
//    t.insert_at(1,{{0,-1},{2,1}}); // 2
//    t.insert_at(2,{{3,-1}}); // 3

    t.compute_weights();

}

Inference::Inference(u_int n_regions, int ploidy, int verbosity) : t(ploidy, n_regions),
                                                                   t_prime(ploidy, n_regions),
                                                                   best_tree(ploidy, n_regions) {

    this->n_regions = n_regions;
    this->ploidy = ploidy;
    this->verbosity = verbosity;
    std::ofstream outfile;
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count(); // get a seed from time
    f_name = std::to_string(seed);
}

Inference::~Inference() {
    destroy();
}

bool Inference::apply_prune_reattach(const vector<vector<double>> &D, const vector<int> &r, bool genotype_preserving,
                                     bool weighted, bool validation_test_mode) {
    /*
     * Applies prune and reattach to t_prime
     * Updates the sums and scores tables partially
     * */

    Node* attached_node;
    try {
        attached_node = t_prime.prune_reattach(genotype_preserving, weighted, validation_test_mode);
    }catch (const std::logic_error& e)
    {
        if (verbosity > 0)
            std::cout << " a logic error was caught during the prune and reattach move, with message '"
                      << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the prune and reattach move, with message '"
                      << e.what() << "'\n";
        return false;
    }

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);
        return true;
    }
    else
        return false;
}

void Inference::compute_t_table(const vector<vector<double>> &D, const vector<int> &r) {

    int j = static_cast<int>(D.size());
    for (int i = 0; i < j; ++i)
    {
        this->t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = this->t.get_children_id_score(this->t.root);

        this->t_scores.push_back(scores_vec);
        this->t_sums.push_back(MathOp::log_sum(scores_vec));
    }

    int m = D.size();
    int n = t.get_n_nodes();
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.prior_score = log_tree_prior(m, n);
    t.posterior_score = log_posterior(t_sum, m, t);

    // update t_prime
    // calls the copy constructor
    t_prime = t;

}

void Inference::destroy() {
    // nothing to deallocate
}

Tree * Inference::comparison(int m, double gamma, unsigned move_id) {
    /*
     * Returns the pointer to the accepted tree
     * m is size(D), i.e. number of cells
     * */

    double t_prime_sum = accumulate(t_prime_sums.begin(), t_prime_sums.end(), 0.0);

    u_int t_n_nodes = t.get_n_nodes();
    u_int t_prime_n_nodes = t_prime.get_n_nodes();

    // compare the posteriors
    t_prime.posterior_score = log_posterior(t_prime_sum, m, t_prime);
    // update the tree prior as well
    t_prime.prior_score = log_tree_prior(m, t_prime_n_nodes);

    // check if move is weighted
    bool weighted;
    if ((move_id %2 == 1) && (move_id <=9)) // weighted move ids
        weighted = true;
    else
        weighted = false;

    // acceptance probability computations
    double score_diff = 0.0;
    if (move_id == 11) // overdispersion change
    {
        double od_score_diff = t_prime.od_score - t.od_score;
        double posterior_score_diff = t_prime.posterior_score - t.posterior_score;

        score_diff = od_score_diff + posterior_score_diff;
    }
    else
        score_diff = t_prime.posterior_score - t.posterior_score;

    double acceptance_prob = exp(gamma*score_diff); // later gets modified

    // compute nbd correction
    double nbd_corr= 1.0;
    double sum_chi=0.0, sum_chi_prime=0.0, sum_omega=0.0, sum_omega_prime=0.0;

    if (move_id == 1) // weighted prune-reattach
    {
        nbd_corr = t.cost() / t_prime.cost();
        if (std::isinf(nbd_corr))
            return &t;

        // ro variable
        acceptance_prob *= nbd_corr;
    }
    else if (move_id == 6 || move_id == 7) // insert/delete move or weighted insert/delete move
    {
        if (tree_prior_in_chi)
        {
            sum_chi = t.chi_insert_delete_reweighted(weighted);
            sum_chi_prime = t_prime.chi_insert_delete_reweighted(weighted);
        }
        else
        {
            vector<double> chi = t.chi_insert_delete(weighted);
            sum_chi = std::accumulate(chi.begin(), chi.end(), 0.0);
            vector<double> chi_prime = t_prime.chi_insert_delete(weighted);
            sum_chi_prime = std::accumulate(chi_prime.begin(), chi_prime.end(), 0.0);
        }

        vector<double> omega = t.omega_insert_delete(lambda_r, lambda_c, weighted);
        sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
        vector<double> omega_prime = t_prime.omega_insert_delete(lambda_r, lambda_c, weighted);
        sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
    }
    else if (move_id == 8 || move_id == 9) // condense/split move or weighted cs
    {

        if (tree_prior_in_chi)
        {
            // new chi weighting
            sum_chi = t.chi_condense_split_reweighted(weighted);
            sum_chi_prime = t_prime.chi_condense_split_reweighted(weighted);
        }
        else
        {
            vector<double> chi = t.chi_condense_split(weighted);
            sum_chi = std::accumulate(chi.begin(), chi.end(), 0.0);
            vector<double> chi_prime = t_prime.chi_condense_split(weighted);
            sum_chi_prime = std::accumulate(chi_prime.begin(), chi_prime.end(), 0.0);
        }

        vector<double> omega = t.omega_condense_split(lambda_s, weighted);
        sum_omega = std::accumulate(omega.begin(), omega.end(), 0.0);
        vector<double> omega_prime = t_prime.omega_condense_split(lambda_s, weighted);
        sum_omega_prime = std::accumulate(omega_prime.begin(), omega_prime.end(), 0.0);
    }

    if (move_id == 6 || move_id == 7 || move_id == 8 || move_id == 9) // moves that require nbd correction
    {
        double v_prime = 0.0;
        if (std::isnan(v)) //v is global
        {
            v = sum_chi / (sum_chi + sum_omega);
            v_prime = sum_chi_prime / (sum_chi_prime + sum_omega_prime);
        }
        else
            v_prime = v;

        if (t_n_nodes < t_prime_n_nodes) // insert, split
        {
            double n = static_cast<double>(t_n_nodes);
            double weight = std::pow(n,2) * std::pow((n+2)/n, n);
            acceptance_prob *= weight;

            // add v' X' v W terms
            double correction = (1.0 - v_prime)*sum_chi/(v*sum_omega_prime);
            acceptance_prob *= correction;
        }
        else // delete, condense
        {
            double n = static_cast<double>(t_n_nodes);
            double weight = std::pow(n,-2) * std::pow(n/(n+2), n);
            acceptance_prob *= weight;

            // add v' X' v W terms
            double correction = v_prime*sum_omega/((1.0-v)*sum_chi_prime);
            acceptance_prob *= correction;
        }
    }

    if (move_id == 7 || move_id == 9) // weighted insert-delete or weighted condense-split
    {
        if (t_n_nodes > t_prime_n_nodes) // delete
        {
            // find the node that is deleted
            // use the get_id_score function since it returns a map having node id as a key
            map<int,double> t_prime_scores = t_prime.get_children_id_score(t_prime.root);
            Node* deleted = nullptr;
            for (auto const &t_node : t.root->get_descendents(false)) // root won't be contained!
            {
                if (!t_prime_scores.count(t_node->id))
                {
                    deleted = t_node;
                    break;
                }
            }
            unsigned d_i_T = deleted->n_descendents;
            unsigned d_i_T_prime = 0;

            // find it's parent in t_prime
            int parent_id = deleted->parent->id;
            for (auto const &t_prime_node : t_prime.root->get_descendents(true))
            {
                if (parent_id == t_prime_node->id)
                {
                    d_i_T_prime = t_prime_node->n_descendents;
                    break;
                }
            }

            double p_add_corr_num = d_i_T_prime + 1;
            double p_add_corr_denom = 2 * d_i_T;
            double ratio = p_add_corr_num / p_add_corr_denom;
            acceptance_prob *= ratio;
        }
        else // insert
        {
            // find the node that is inserted in t_prime
            map<int,double> t_scores = t.get_children_id_score(t.root);

            Node* added = nullptr;
            for (auto const &t_prime_node : t_prime.root->get_descendents(false)) // without root
            {
                if (!t_scores.count(t_prime_node->id))
                {
                    added = t_prime_node;
                    break;
                }
            }

            unsigned d_i_T_prime = added->n_descendents;
            unsigned d_i_T = 0;
            // find it's parent in t
            int parent_id = added->parent->id;
            for (auto const &t_node : t.root->get_descendents(true))
            {
                if (parent_id == t_node->id)
                {
                    d_i_T = t_node->n_descendents;
                    break;
                }
            }

            double p_add_corr_num = 2 * d_i_T_prime;
            double p_add_corr_denom = d_i_T + 1;
            double ratio = p_add_corr_num / p_add_corr_denom;
            acceptance_prob *= ratio;
        }
    }

    assert(!std::isnan(acceptance_prob));

    if (verbosity > 0)
        cout<<"acceptance prob: "<<acceptance_prob<<endl;

    if (acceptance_prob > 1)
        return &t_prime;
    else
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        boost::random::uniform_real_distribution<double> distribution(0.0,1.0);
        double rand_val = distribution(gen);

        if (verbosity > 0)
            cout<<"rand_val: "<<rand_val<<endl;

        if (acceptance_prob > rand_val)
            return &t_prime;
        else
            return &t;
    }
}

void Inference::infer_mcmc(const vector<vector<double>> &D, const vector<int> &r, const vector<float> &move_probs,
                           int n_iters,
                           unsigned int size_limit) {

    int m = static_cast<int>(D.size());

    u_int n_rejected = 0;
    u_int n_accepted = 0;

    double gamma = 1.0; // gamma param to amplify the difference in log likelihoods

    // for writing the posteriors on file
    std::ofstream mcmc_scores_file;
    std::ofstream rel_mcmc_scores_file;
    if (verbosity > 1)
    {
        mcmc_scores_file.open(f_name + "_markov_chain.txt", std::ios_base::app);
        rel_mcmc_scores_file.open(f_name + "_rel_markov_chain.txt", std::ios_base::app);
    }

    best_tree = t; //start with the t

    for (int i = 0; i < n_iters; ++i) {


        bool rejected_before_comparison = false;

        std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
        boost::random::discrete_distribution<> d(move_probs.begin(), move_probs.end());

        unsigned move_id = d(gen);

        switch (move_id)
        {
            case 0:
            {
                // prune & reattach
                if (verbosity > 0)
                    cout << "Prune and reattach" << endl;
                bool prune_reattach_success = apply_prune_reattach(D, r, false, false, false);
                if (not prune_reattach_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Prune and reattach is rejected before comparison"<<endl;
                }
                break;
            }
            case 1:
            {
                // weighted prune & reattach
                if (verbosity > 0)
                    cout<<"Weighted prune and reattach"<<endl;
                bool weighted_prune_reattach_success = apply_prune_reattach(D, r, false, true, false); // weighted=true
                if (not weighted_prune_reattach_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Weighted prune and reattach is rejected before comparison"<<endl;
                }
                break;
            }
            case 2:
            {
                // swap labels
                if (verbosity > 0)
                    cout << "swap labels" << endl;
                bool swap_success = apply_swap(D, r, false); // weighted=false
                if (not swap_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Swap labels is rejected before comparison" << endl;
                }
                break;
            }
            case 3:
                {
                // weighted swap labels
                if (verbosity > 0)
                    cout << "weighted swap labels" << endl;
                bool weighted_swap_success = apply_swap(D, r, true); // weighted=true
                if (not weighted_swap_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Weighted swap labels is rejected before comparison" << endl;
                }
                break;
            }
            case 4:
            {
                // add or remove event
                if (verbosity > 0)
                    cout << "add or remove event" << endl;
                // pass 0.0 to the poisson distributions to have 1 event added/removed
                bool add_remove_success = apply_add_remove_events(lambda_r, lambda_c, D, r, false); // weighted=false
                if (not add_remove_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Add or remove event is rejected before comparison"<<endl;
                }
                break;
            }
            case 5:
            {
                // weighted add or remove event
                if (verbosity > 0)
                    cout << "weighted add or remove event" << endl;
                // pass 0.0 to the poisson distributions to have 1 event added/removed
                bool add_remove_success = apply_add_remove_events(lambda_r, lambda_c, D, r, true); // weighted=true
                if (not add_remove_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Weighted add or remove event is rejected before comparison"<<endl;
                }
                break;
            }
            case 6:
            {
                // insert delete node
                if (verbosity > 0)
                    cout << "insert/delete node" << endl;
                bool insert_delete_success = apply_insert_delete_node(lambda_r, lambda_c, D, r, size_limit, false); // weighted=false
                if (not insert_delete_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "insert/delete rejected before comparison" << endl;
                }
                break;
            }
            case 7:
            {
                // weighted insert delete node
                if (verbosity > 0)
                    cout << "weighted insert/delete node" << endl;
                bool insert_delete_success = apply_insert_delete_node(lambda_r, lambda_c, D, r, size_limit, true); // weighted=true
                if (not insert_delete_success) {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "weighted insert/delete rejected before comparison" << endl;
                }
                break;
            }
            case 8:
            {
                // condense split move
                if (verbosity > 0)
                    cout << "condense split move " <<endl;
                bool condense_split_success = apply_condense_split(lambda_s, D, r, size_limit, false);
                if (not condense_split_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "condense/split move is rejected before comparison"<<endl;
                }
                break;
            }
            case 9:
            {
                // weighted condense split move
                if (verbosity > 0)
                    cout << "weighted condense split move " <<endl;
                bool condense_split_success = apply_condense_split(lambda_s, D, r, size_limit, true); //weighted=true
                if (not condense_split_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "weighted condense/split move is rejected before comparison"<<endl;
                }
                break;
            }
            case 10:
            {
                // genotype_preserving prune & reattach
                if (verbosity > 0)
                    cout<<"Genotype preserving prune and reattach"<<endl;
                bool genotype_prune_reattach_success = apply_prune_reattach(D, r, true, false, false); // weighted=false
                if (not genotype_prune_reattach_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Genotype preserving prune/reattach is rejected before comparison"<<endl;
                }
                break;
            }
            case 11:
            {
                // changing overdispersion move
                if (verbosity > 0)
                    cout<<"Overdispersion changing move"<<endl;
                bool od_change_success = apply_overdispersion_change(D, r);
                if (not od_change_success)
                {
                    rejected_before_comparison = true;
                    if (verbosity > 0)
                        cout << "Overdispersion changing move is rejected before comparison"<<endl;
                }
                break;
            }
            default:
                throw std::logic_error("undefined move index");
        }

        // compare the trees

        Tree* accepted;
        if (rejected_before_comparison)
            accepted = &t;
        else
            accepted = comparison(m, gamma, move_id);

        static double first_score = accepted->posterior_score + accepted->od_score; // the first value will be kept in whole program
        // print accepted log_posterior
        mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->od_score << ',';
        rel_mcmc_scores_file << std::setprecision(print_precision) << accepted->posterior_score + accepted->od_score - first_score << ',';

        // update trees and the matrices
        if (accepted == &t_prime)
        {
            n_accepted++;
            t_sums = t_prime_sums;
            update_t_scores(); // this should be called before t=tprime, because it checks the tree sizes in both.
            t = t_prime;
            if ((t_prime.posterior_score + t_prime.od_score)  > (best_tree.posterior_score + best_tree.od_score))
                best_tree = t_prime;
        }
        else
        {
            n_rejected++;
            t_prime = t;
        }
        t_prime_sums.clear();
        t_prime_scores.clear();

        // N: n_nodes of tree t
        if ((i > 10000) && (i % 10000 == 0))
        {
            cout << "iteration" << i <<  "tree score" << t.posterior_score + t.od_score  << endl;
        }
    }

}

bool Inference::apply_overdispersion_change(const vector<vector<double>> &D, const vector<int> &r) {
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

        t_prime.nu = std::exp(log_t_prime_nu); // real space

        this->compute_t_prime_od_scores(D,r);
        compute_t_prime_scores(t_prime.root, D, r);
        compute_t_prime_sums(D);
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the apply overdispersion change move, with message '"
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

    double log_prior = - (n -1 + m) * log(n+1); // tree prior
    return log_prior;
}

double Inference::log_posterior(double tree_sum, int m, Tree &tree) {
    // TODO: move to the mathop
    // m: n_cells, n: n_nodes

    int n = tree.get_n_nodes();
    double log_posterior = 0.0;
    log_posterior = tree_sum + this->log_tree_prior(m, n); // initialise posterior with prior then add posterior

    int repetition_count = 0; // the repetition count to be used in the penalisation
    double c_penalisation = c_penalise; // the penalisation coefficient, global var

    /*
     * compute penalization term
     * K: max region index
     *
     * */

    vector<double> p_v;
    for (auto it = tree.all_nodes_vec.begin()+1; it != tree.all_nodes_vec.end(); ++it) // without the root
    {
        Node* node = *it;
        map<u_int,int>& c_change = node->c_change;
        int v = 0;
        int v_prev = 0; // the first region is zero
        int i_prev = -1; // the initial index is -1, it'll be updated later

        auto last_elem_id = c_change.rbegin()->first;

        for (auto const &event_it : c_change)
        {

            // penalisation for repetition
            int parent_state = 0;

            try
            {
                parent_state = node->parent->c.at(event_it.first);
                int c_change_val = event_it.second;

                if (signbit(c_change_val) != signbit(parent_state))
                    repetition_count++;
            }
            catch (const std::out_of_range& e)
            {
                // pass
            }

            int diff;
            if (static_cast<int>(event_it.first) - 1 != i_prev) // if the region is adjacent to its previous
            {
                int diff_right = 0 - v_prev; // the right hand side change at the end of the last consecutive region
                if (diff_right > 0)
                    v += diff_right;
                v_prev = 0;
            }
            diff = event_it.second - v_prev;
            if (diff > 0)
                v += diff;
            v_prev = event_it.second;
            i_prev = event_it.first;

            if (event_it.first == last_elem_id)
            {
                int diff_last = 0 - v_prev;

                if (diff_last > 0)
                {
                    v += diff_last;
                    assert(v>0);
                }
            }
        }

        double pv_i = 0.0;

        int K = this->n_regions;
        pv_i -= v*log(2*K); // the event prior
        p_v.push_back(pv_i);

    }

    assert(n==static_cast<int>(p_v.size()));
    double PV = 0.0;
    PV += std::accumulate(p_v.begin(), p_v.end(), 0.0);
    PV -= Lgamma::get_val(n+1);

    log_posterior += PV;

    log_posterior -= c_penalisation*repetition_count; // penalise the repetitions

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


void Inference::write_best_tree() {
    std::ofstream outfile;
    outfile.open(f_name+"_tree.txt", std::ios_base::app);
    outfile << "The resulting tree is: "<<std::endl;
    outfile << std::setprecision(print_precision) << best_tree;
    std::cout << std::setprecision(print_precision) << best_tree;
}

void Inference::compute_t_prime_scores(Node *attached_node, const vector<vector<double>> &D, const vector<int> &r) {

    // if the t_prime_scores is empty then fill it with size(D) elements
    // otherwise append the results to the first size(D) places in t_prime_scores
    // size of t_prime_scores should always be size(D)

    bool is_empty_table = t_prime_scores.empty();

    int j = 0;
    for (auto const &d: D)
    {
        double sum_d = accumulate( d.begin(), d.end(), 0.0);

        if (attached_node != t_prime.root)
            attached_node->parent->log_score = t_scores[j][attached_node->parent->id]; // the indices must match
        // attached node->parent->id must match the all_nodes_vec index
        t_prime.compute_stack(attached_node, d, sum_d,r);

        if (is_empty_table)
            t_prime_scores.push_back(t_prime.get_children_id_score(attached_node));
        else
        {
            // append the contents of second hashmap into the first
            // note that they cannot have overlapping keys
            for (auto const& map_item : t_prime.get_children_id_score(attached_node))
            {
                t_prime_scores[j][map_item.first] = map_item.second;
            }
        }

        j++;
    }
}

bool Inference::apply_swap(const vector<vector<double>> &D, const vector<int> &r, bool weighted, bool test_mode) {

    vector<Node*> swapped_nodes;
    try {
        swapped_nodes = t_prime.swap_labels(weighted, test_mode);
    }catch (const std::logic_error& e)
    {
        if (verbosity > 0)
            std::cout << " a logic error was caught during the swap labels move, with message '"
                      << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the swap labels move, with message '"
                      << e.what() << "'\n";
        return false;
    }


    if (swapped_nodes.empty()) // it can be empty if an exception is thrown or the move is rejected
        return false;

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

    int i = 0;
    for (auto const &d: D) {
        vector<double> old_vals;
        old_vals.reserve(t_scores[i].size()); // the max possible size
        vector<double> new_vals;
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

        t_prime_sums.push_back(res);
        i++;
    }
}

bool Inference::apply_add_remove_events(double lambda_r, double lambda_c, const vector<vector<double>> &D,
                                        const vector<int> &r, bool weighted,
                                        bool validation_test_mode)
{
    /*
     * Applies add/remove event to t_prime
     * Updates the sums and scores tables partially
     * */

    // weighted = false
    Node* attached_node;

    try {
        attached_node = t_prime.add_remove_events(lambda_r,lambda_c,weighted, validation_test_mode);
    }catch (const std::logic_error& e)
    {
        if (verbosity > 0)
            std::cout << " a logic error was caught during the add remove events move, with message '"
                      << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the add remove events move, with message '"
                      << e.what() << "'\n";
        return false;
    }

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}

bool Inference::apply_insert_delete_node(double lambda_r, double lambda_c, const vector<vector<double>> &D,
                                         const vector<int> &r, unsigned int size_limit, bool weighted) {
    /*
     * Applies the insert/delete move on t_prime
     * Updates the sums and scores tables partially
     * */

    Node* tobe_computed;
    try {
        tobe_computed = t_prime.insert_delete_node(lambda_r, lambda_c, size_limit, weighted);
    }catch (const std::out_of_range& e)
    {
        if (verbosity > 0)
            std::cout << " an out of range error was caught during the insert/delete node move, with message '"
                      << e.what() << "'\n";
        return false; // reject the move
    }catch (const std::exception& e) {
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the insert/delete node move, with message '"
                      << e.what() << "'\n";
        return false;
    }

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

bool Inference::apply_condense_split(double lambda_s, const vector<vector<double>> &D, const vector<int> &r,
                                     unsigned int size_limit, bool weighted) {
    /*
     * Applies the condense/delete move on t_prime
     * Updates the sums and scores tables partially
     * */

    Node* tobe_computed;
    try
    {
        tobe_computed = t_prime.condense_split_node(lambda_s, size_limit, weighted);
    }catch (const std::exception& e) {
        if (verbosity > 0)
            std::cout << " a standard exception was caught during the split/condense node move, with message '"
                      << e.what() << "'\n";
        return false;
    }

    if (tobe_computed != nullptr)
    {
        compute_t_prime_scores(tobe_computed, D, r);
        compute_t_prime_sums(D);
        return true;
    }
    else
        return false;
}

vector<vector<int>> Inference::assign_cells_to_nodes(const vector<vector<double>> &D, const vector<int> &r) {


    /*
     * re-computes the best tree to assigns cells to nodes by the maximum score
     * Writes the cell-node_ids, region sizes, and the inferred CNVs files.
     * Returns the inferred CNVs matrix (with the ploidy added)
     * */

    t_scores.clear();
    t_sums.clear();
    t = best_tree;
    this->compute_t_table(D,r);

    std::ofstream cell_node_ids_file;
    std::ofstream cell_node_cnvs_file;
    std::ofstream region_sizes_file;

    if (verbosity > 1)
    {
        cell_node_ids_file.open(f_name + "_cell_node_ids.txt");
        cell_node_cnvs_file.open(f_name + "_cell_node_cnvs.txt");
        region_sizes_file.open(f_name + "_region_sizes.txt");

        for (const auto &r_it : r) region_sizes_file << r_it << "\n";
    }

    // create a hashmap of nodes for constant access by id
    unordered_map<uint64_t , Node*> hash_map;
    for (unsigned i=0; i < t.all_nodes_vec.size(); i++)
    {
        hash_map[t.all_nodes_vec[i]->id] = t.all_nodes_vec[i];
    }

    size_t n_cells = t_scores.size();

    vector<vector<int>> cell_regions(n_cells, vector<int>(this->n_regions, ploidy)); //fill ctor, initialise with ploidy

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
            cell_regions[j][x.first] = x.second + ploidy;
        }

    }
    for (size_t k = 0; k < n_cells; ++k) {
        for (u_int i = 0; i < n_regions; ++i) {
            if (verbosity > 1)
                cell_node_cnvs_file << cell_regions[k][i] << '\t';
        }
        if (verbosity > 1)
            cell_node_cnvs_file << '\n';
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
Inference::get_tree_od_root_scores(const vector<vector<double>> &D, const vector<int> &r, const Tree &tree) {
    /*
     * Returns the vector of overdispersed root scores for every cell
     * */

    u_int n_cells = D.size();

    vector<double> scores(n_cells);
    for (u_int i = 0; i < n_cells; ++i) {
        double sum_d = std::accumulate(D[i].begin(), D[i].end(), 0.0);
        scores[i] = tree.get_od_root_score(r, sum_d, D[i]);
    }

    return scores;
}

void Inference::compute_t_od_scores(const vector<vector<double>> &D, const vector<int> &r) {
/*
 * Computes and stores the overdispersed root sum for tree t across all cells
 * */
    vector<double> t_od_root_scores = get_tree_od_root_scores(D,r,this->t);
    double sum_t_od_root_scores = std::accumulate(t_od_root_scores.begin(), t_od_root_scores.end(), 0.0);
    this->t.od_score = sum_t_od_root_scores;

}

void Inference::compute_t_prime_od_scores(const vector<vector<double>> &D, const vector<int> &r) {
/*
 * Computes and stores the overdispersed root sum for tree t prime across all cells
 * */
    vector<double> t_prime_od_root_scores = get_tree_od_root_scores(D,r,this->t_prime);
    double sum_t_prime_od_root_scores = std::accumulate(t_prime_od_root_scores.begin(), t_prime_od_root_scores.end(), 0.0);
    this->t_prime.od_score = sum_t_prime_od_root_scores;
}


#endif //SC_DNA_INFERENCE_H
