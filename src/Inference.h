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

class Inference {
/*
 * Contains functionality to perform monte carlo markov chains (mcmc) inference
 * //TODO use stack dynamic Trees for t & t_prime
 * */
public:
    Tree *t;
    Tree *t_prime;
    Tree *best_tree;
    std::vector<std::vector<double>> t_scores;
    std::vector<double> t_sums;
    std::vector<unordered_map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;
    std::string f_name;

public:
    Inference(u_int n_regions, u_int ploidy=2);
    ~Inference();
    void destroy();
    void compute_t_table(const vector<vector<int>> &D, const vector<int>& r);
    void compute_t_prime_scores(Node *attached_node, const vector<vector<int>> &D, const vector<int> &r);
    void compute_t_prime_sums(const vector<vector<int>> &D);
    double log_posterior(double tree_sum, int m, int n);

    bool apply_prune_reattach(const vector<vector<int>> &D, const vector<int> &r, bool weighted = false,
                              bool validation_test_mode = false);
    bool apply_add_remove_events(float lambda_r, float lambda_c, const vector<vector<int>> &D, const vector<int> &r,
                                 bool weighted = false,
                                 bool validation_test_mode = false);


    bool apply_swap(const vector<vector<int>> &D, const vector<int> &r, bool weighted = false, bool test_mode = false);
    Tree * comparison(int m);

    void infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs);
    void write_best_tree();
    void update_t_scores();

    Tree *get_t() const;

    Tree *get_t_prime() const;
    void random_initialize(); // randomly initializes a tree and copies it into the other
    void initialize_worked_example(); // initializes the trees based on the test example
};



void Inference::random_initialize() {

    t->random_insert({{0, 1}, {1, 1}});
    t->random_insert({{1, 1}, {2, 1}});
    t->random_insert({{0, -1}});
    t->random_insert({{3, -1}});
    t->random_insert({{1, 1}});

    t->compute_weights();

}

void Inference::initialize_worked_example() {

    // build tree
    // tree that generated the data
    t->random_insert({{0, 1}, {1, 1}});
    t->insert_at(1,{{1, 1}, {2, 1}});
    t->insert_at(2,{{0, -1}});
    t->insert_at(2,{{3, -1}});
    t->insert_at(1,{{1, 1}});

    // score the -2596.33 tree
//    t->random_insert({{0, 1}, {1, 1}}); //1
//    t->insert_at(1,{{1, 1}}); // 2
//    t->insert_at(1,{{0, -1}}); // 3
//    t->insert_at(3,{{1, 1}, {2, 1}}); // 4
//    t->insert_at(4,{{3, -1}}); // 5

    t->compute_weights();

}

Inference::Inference(u_int n_regions, u_int ploidy) {
    t = new Tree(ploidy, n_regions);
    t_prime =new Tree(ploidy, n_regions);
    best_tree = new Tree(ploidy, n_regions);

    std::ofstream outfile;
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count(); // get a seed from time
    f_name = std::to_string(seed) + ".txt";


}

Inference::~Inference() {
    destroy();
}

bool Inference::apply_prune_reattach(const vector<vector<int>> &D, const vector<int> &r, bool weighted,
                                     bool validation_test_mode) {
    /*
     * Applies prune and reattach to t_prime
     * Updates the sums and scores tables partially
     * */

    using namespace std;
    // weighted = false
    Node* attached_node = t_prime->prune_reattach(weighted, validation_test_mode);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}

void Inference::compute_t_table(const vector<vector<int>> &D, const vector<int>& r) {

    int n = static_cast<int>(D.size());
    for (int i = 0; i < n; ++i)
    {
        this->t->compute_tree(D[i], r);
        vector<double> scores_vec = this->t->get_scores();

        this->t_scores.push_back(scores_vec);
        this->t_sums.push_back(MathOp::log_sum(scores_vec));
    }

    int m = size(D);
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    int t_n = t->get_n_nodes();
    t->score = log_posterior(t_sum, m, t_n);

    // update t_prime
    // calls the copy constructor
    *t_prime = *t;

}

void Inference::destroy() {
    delete t;
    t = nullptr;
    delete t_prime;
    t_prime = nullptr;
}

Tree * Inference::comparison(int m) {
    /*
     * Returns the pointer to the accepted tree
     * m is size(D)
     * */


    double log_post_t = 0.0;
    double log_post_t_prime = 0.0;

    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    int t_n = t->get_n_nodes();
    log_post_t = log_posterior(t_sum, m, t_n);

    // assign the tree score
    t->score = log_post_t;

    double t_prime_sum = accumulate( t_prime_sums.begin(), t_prime_sums.end(), 0.0);
    int tp_n = t_prime->get_n_nodes();
    log_post_t_prime = log_posterior(t_prime_sum, m, tp_n);

    t_prime->score = log_post_t_prime;


    double acceptance_prob = exp(log_post_t_prime - log_post_t);

    cout<<"acceptance prob: "<<acceptance_prob<<endl;

    if (acceptance_prob > 1)
        return t_prime;

    else
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_generator();
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        double rand_val = distribution(gen);

        cout<<"rand_val: "<<rand_val<<endl;

        if (acceptance_prob > rand_val)
            return t_prime;
        else
            return t;

    }
}

void Inference::infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs) {


    int m = static_cast<int>(D.size());
    int n_accepted = 0;
    int n_rejected = 0;
    int n_attached_to_the_same_pos = 0;
    int add_remove_move_rejected = 0;

    // for writing the posteriors on file
    std::ofstream outfile;
    outfile.open(f_name, std::ios_base::app);

     *best_tree = *t; //start with the t

    for (int i = 0; i < 1000; ++i) {


        bool rejected_before_comparison = false;

        std::mt19937 &gen = SingletonRandomGenerator::get_generator();
        std::discrete_distribution<> d(move_probs.begin(), move_probs.end());

        unsigned move_id = d(gen);

        switch (move_id)
        {
            case 0:
            {
                // prune & reattach
                cout << "Prune and reattach" << endl;
                bool prune_reattach_success = apply_prune_reattach(D, r, false);
                if (not prune_reattach_success) {
                    n_attached_to_the_same_pos++;
                    rejected_before_comparison = true;
                }
                break;
            }
            case 1:
            {
                // weighted prune & reattach
                cout<<"Weighted prune and reattach"<<endl;
                bool weighted_prune_reattach_success = apply_prune_reattach(D, r, true); // weighted=true
                if (not weighted_prune_reattach_success)
                {
                    n_attached_to_the_same_pos++;
                    rejected_before_comparison = true;
                }
                break;
            }
            case 2:
            {
                // swap labels
                cout << "swap labels" << endl;
                bool swap_success = apply_swap(D, r, false); // weighted=false
                if (not swap_success)
                    rejected_before_comparison = true;
                break;
            }
            case 3:
                {
                // weighted swap labels
                cout << "weighted swap labels" << endl;
                bool weighted_swap_success = apply_swap(D, r, true); // weighted=true
                if (not weighted_swap_success)
                    rejected_before_comparison = true;
                break;
            }
            case 4:
            {
                // TODO: start allowing this move after considering the forbidden states
                // add or remove event
                cout << "add or remove event" << endl;
                // pass 0.0f to the poisson distributions to have 1 event added/removed
                bool add_remove_success = apply_add_remove_events(0.0f, 0.0f, D, r, true); // weighted=true
                if (not add_remove_success) {
                    add_remove_move_rejected++;
                    rejected_before_comparison = true;
                }
                break;
            }
            default:
                throw std::logic_error("undefined move index");
        }

        // compare the trees
        Tree* accepted;
        if (rejected_before_comparison)
            accepted = t;
        else
            accepted = comparison(m);

        // print accepted log_posterior
        outfile << std::setprecision(8) << accepted->score << ',';

        // update trees and the matrices
        if (accepted == t_prime)
        {
            n_accepted++;

            t_sums = t_prime_sums;
            update_t_scores();
            *t = *t_prime;
            if (t_prime->score > best_tree->score)
                *best_tree = *t_prime;
        }
        else
        {
            n_rejected++;
            *t_prime = *t;
        }
        t_prime_sums.clear();
        t_prime_scores.clear();

    }
    cout<<"n_accepted: "<<n_accepted<<endl;
    cout<<"n_rejected: "<<n_rejected<<endl;
    cout<<"n_attached_to_the_same_pos: "<<n_attached_to_the_same_pos<<endl;
    cout<<"add_remove_move_rejected: "<<add_remove_move_rejected<<endl;
}

double Inference::log_posterior(double tree_sum, int m, int n) {

    // TODO: implement this then call this to compute log posteriors
    // TODO: set the log posterior of the best tree to the initial one first
    // later updated it upon acceptance of t_primes
    double log_posterior = 0.0;
    log_posterior = tree_sum - (n -1 + m ) * log(n+1);

    return log_posterior;
}

void Inference::update_t_scores() {

    for (unsigned k=0; k < t_scores.size(); k++)
        for (unsigned i=0; i < t_scores[k].size(); i++)
            if(t_prime_scores[k].find(i) != t_prime_scores[k].end()) // if found in hashmap
                t_scores[k][i] = t_prime_scores[k][i];

}

Tree *Inference::get_t() const {
    return t;
}

Tree *Inference::get_t_prime() const {
    return t_prime;
}

void Inference::write_best_tree() {
    std::ofstream outfile;
    outfile.open(f_name, std::ios_base::app);
    outfile << "The resulting tree is: "<<std::endl;
    outfile << std::setprecision(8) << *best_tree;
    std::cout << std::setprecision(8) << *best_tree;
}

void Inference::compute_t_prime_scores(Node *attached_node, const vector<vector<int>> &D, const vector<int> &r) {

    // if the t_prime_scores is empty then fill it with size(D) elements
    // otherwise append the results to the first size(D) places in t_prime_scores
    // size of t_prime_scores should always be size(D)

    bool is_empty_table = t_prime_scores.empty();

    int j = 0;
    for (auto const &d: D)
    {
        int sum_d = accumulate( d.begin(), d.end(), 0);
        attached_node->parent->log_score = t_scores[j][attached_node->parent->id];
        t_prime->compute_stack(attached_node, d, sum_d,r);

        if (is_empty_table)
            t_prime_scores.push_back(t_prime->get_children_id_score(attached_node));
        else
        {
            // append the contents of second hashmap into the first
            // note that they cannot have overlapping keys
            for (auto const& map_item : t_prime->get_children_id_score(attached_node))
            {
                t_prime_scores[j][map_item.first] = map_item.second;
            }
        }

        j++;
    }
}

bool Inference::apply_swap(const vector<vector<int>> &D, const vector<int> &r, bool weighted, bool test_mode) {

    vector<Node*> swapped_nodes = t_prime->swap_labels(weighted, test_mode);

    if (swapped_nodes.empty())
        return false;

    for (auto const &node : swapped_nodes)
    {
        compute_t_prime_scores(node, D, r);
    }
    compute_t_prime_sums(D);

    return true;

}

void Inference::compute_t_prime_sums(const vector<vector<int>> &D) {

    int i = 0;
    for (auto const &d: D)
    {
        vector<double> old_vals;
        vector<double> new_vals;

        for (auto &u_map : t_prime_scores[i])
        {
            old_vals.push_back(t_scores[i][u_map.first]);
            new_vals.push_back(u_map.second);
        }

        double res =MathOp::log_replace_sum(t_sums[i],old_vals,new_vals);
        assert(!isnan(res));

        t_prime_sums.push_back(res);
        i++;
    }
}

bool Inference::apply_add_remove_events(float lambda_r, float lambda_c, const vector<vector<int>> &D,
                                        const vector<int> &r, bool weighted,
                                        bool validation_test_mode) {
    /*
     * Applies add/remove event to t_prime
     * Updates the sums and scores tables partially
     * */

    // weighted = false
    Node* attached_node = t_prime->add_remove_events(lambda_r,lambda_c,weighted, validation_test_mode);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);
        compute_t_prime_sums(D);

        return true;
    }
    else
        return false;
}


#endif //SC_DNA_INFERENCE_H
