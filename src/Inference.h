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

class Inference {
/*
 * Contains functionality to perform monte carlo markov chains (mcmc) inference
 * //TODO use stack dynamic Trees for t & t_prime
 * */
public:
    Tree t;
    Tree t_prime;
    Tree best_tree;
    u_int n_regions;
    std::vector<std::map<int, double>> t_scores;
    std::vector<double> t_sums;
    std::vector<std::map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;
    std::string f_name;

public:
    Inference(u_int n_regions, u_int ploidy=2);
    ~Inference();
    void destroy();
    void compute_t_table(const vector<vector<int>> &D, const vector<int>& r);
    void compute_t_prime_scores(Node *attached_node, const vector<vector<int>> &D, const vector<int> &r);
    void compute_t_prime_sums(const vector<vector<int>> &D);
    double log_posterior(double tree_sum, int m, Tree &tree);
    bool apply_prune_reattach(const vector<vector<int>> &D, const vector<int> &r, bool weighted = false,
                              bool validation_test_mode = false);
    bool apply_add_remove_events(double lambda_r, double lambda_c, const vector<vector<int>> &D, const vector<int> &r,
                                 bool weighted = false,
                                 bool validation_test_mode = false);
    bool apply_insert_delete_node(double lambda_r, double lambda_c, const vector<vector<int>> &D, const vector<int> &r, bool weighted=false, bool validation_test_mode=false);
    bool apply_swap(const vector<vector<int>> &D, const vector<int> &r, bool weighted = false, bool test_mode = false);
    Tree * comparison(int m);
    void infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs, int n_iters);
    void write_best_tree();
    void update_t_scores();
    void random_initialize(); // randomly initializes a tree and copies it into the other
    void initialize_worked_example(); // initializes the trees based on the test example
private:
    int deleted_node_idx();
};




void Inference::random_initialize() {

    t.random_insert({{0, 1}, {1, 1}});
    t.random_insert({{1, 1}, {2, 1}});
    t.random_insert({{0, -1}});
    t.random_insert({{3, -1}});
    t.random_insert({{1, 1}});

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

Inference::Inference(u_int n_regions, u_int ploidy): t(ploidy, n_regions), t_prime(ploidy, n_regions), best_tree(ploidy, n_regions)  {
    this->n_regions = n_regions;



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

    Node* attached_node;
    try {
        attached_node = t_prime.prune_reattach(weighted, validation_test_mode);
    }catch (const std::logic_error& e)
    {
        std::cout << " a logic error was caught during the prune and reattach move, with message '"
                  << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
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

void Inference::compute_t_table(const vector<vector<int>> &D, const vector<int>& r) {

    int n = static_cast<int>(D.size());
    for (int i = 0; i < n; ++i)
    {
        this->t.compute_tree(D[i], r);
        std::map<int, double> scores_vec = this->t.get_children_id_score(this->t.root);

        this->t_scores.push_back(scores_vec);
        this->t_sums.push_back(MathOp::log_sum(scores_vec));
    }

    int m = size(D);
    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    t.score = log_posterior(t_sum, m, t);

    // update t_prime
    // calls the copy constructor
    t_prime = t;

}

void Inference::destroy() {

}

Tree* Inference::comparison(int m) {
    /*
     * Returns the pointer to the accepted tree
     * m is size(D)
     * */


    double log_post_t = 0.0;
    double log_post_t_prime = 0.0;

    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    log_post_t = log_posterior(t_sum, m, t);

    // assign the tree score
    t.score = log_post_t;

    double t_prime_sum = accumulate( t_prime_sums.begin(), t_prime_sums.end(), 0.0);
    log_post_t_prime = log_posterior(t_prime_sum, m, t_prime);

    t_prime.score = log_post_t_prime;

    double acceptance_prob = exp(log_post_t_prime - log_post_t);

    if (log_post_t_prime > -2000.0)
        cout<<"debug"; // remove it afterwards

    cout<<"acceptance prob: "<<acceptance_prob<<endl;

    if (acceptance_prob > 1)
        return &t_prime;

    else
    {
        std::mt19937 &gen = SingletonRandomGenerator::get_generator();
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        double rand_val = distribution(gen);

        cout<<"rand_val: "<<rand_val<<endl;

        if (acceptance_prob > rand_val)
            return &t_prime;
        else
            return &t;

    }
}

void Inference::infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs, int n_iters) {


    int m = static_cast<int>(D.size());
    int n_accepted = 0;
    int n_rejected = 0;
    int n_attached_to_the_same_pos = 0;
    int add_remove_move_rejected = 0;
    int insert_delete_move_rejection = 0;

    // for writing the posteriors on file
    std::ofstream outfile;
    outfile.open(f_name, std::ios_base::app);

    best_tree = t; //start with the t

    for (int i = 0; i < n_iters; ++i) {


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
                // add or remove event
                cout << "add or remove event" << endl;
                // pass 0.0 to the poisson distributions to have 1 event added/removed
                bool add_remove_success = apply_add_remove_events(0.0, 0.0, D, r, true); // weighted=true
                if (not add_remove_success) {
                    add_remove_move_rejected++;
                    rejected_before_comparison = true;
                }
                break;
            }
            case 5:
            {
                // insert delete node
                cout << "insert/delete node" << endl;
                bool insert_delete_success = apply_insert_delete_node(1.0, 1.0, D, r, true, false); // weighted=false
                if (not insert_delete_success) {
                    insert_delete_move_rejection++;
                    rejected_before_comparison = true;
                    cout << "insert/delete rejected before comparison" << endl;
                }
                else
                    cout<< "insert/delete accepted!" << endl;
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
            accepted = comparison(m);

        // print accepted log_posterior
        outfile << std::setprecision(8) << accepted->score << ',';

        // update trees and the matrices
        if (accepted == &t_prime)
        {
            n_accepted++;

            t_sums = t_prime_sums;
            update_t_scores(); // this should be called before t=tprime, because it checks the tree sizes in both.
            t = t_prime;
            if (t_prime.score > best_tree.score)
                best_tree = t_prime;
        }
        else
        {
            n_rejected++;
            t_prime = t;
        }
        t_prime_sums.clear();
        t_prime_scores.clear();

    }
    cout<<"n_accepted: "<<n_accepted<<endl;
    cout<<"n_rejected: "<<n_rejected<<endl;
    cout<<"n_attached_to_the_same_pos: "<<n_attached_to_the_same_pos<<endl;
    cout<<"add_remove_move_rejected: "<<add_remove_move_rejected<<endl;
}

double Inference::log_posterior(double tree_sum, int m, Tree &tree) {
    // TODO: move to the mathop
    // m: n_cells, n: n_nodes

    int n = tree.get_n_nodes();

    double log_posterior = 0.0;
    log_posterior = tree_sum - (n -1 + m ) * log(n+1);

    /*
     * compute penalization term
     * K: max region index
     * V: the event vector (vector of labels) list of hashmap c_change
     *
     *
     *
     * compute K from v, v is c_change
     * 1) compute all v! and put in hashmap
     *
     * */




    map<int, double> vfact_hash;

    int K = this->n_regions;
    for (auto it = tree.all_nodes_vec.begin()+1; it != tree.all_nodes_vec.end(); ++it)
    {
        Node* node = *it;
        map<u_int,int>& c_change = node->c_change;
        int v = 0;
        for (auto const &it : c_change)
            v += it.second;

        vfact_hash[v] = log(tgamma(abs(v)+1)); // log of factorial
    }

    vector<double> p_v;
    for (auto it = tree.all_nodes_vec.begin()+1; it != tree.all_nodes_vec.end(); ++it)
    {
        Node* node = *it;
        map<u_int,int>& c_change = node->c_change;
        int v = 0;
        for (auto const &it : c_change)
            v += abs(it.second);
        double pv_i = 0.0;

        pv_i += vfact_hash[v];
        pv_i -= v*log(2*K);

        for (auto const &it : c_change)
            pv_i -= log(tgamma(abs(it.second) + 1)); // +1 because we are using gamma func for factorial

        p_v.push_back(pv_i);

    }

    assert(n==p_v.size());
    double PV = 0.0;
    PV += std::accumulate(p_v.begin(), p_v.end(), 0.0);
    PV -= log(tgamma(n+1));


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


void Inference::write_best_tree() {
    std::ofstream outfile;
    outfile.open(f_name, std::ios_base::app);
    outfile << "The resulting tree is: "<<std::endl;
    outfile << std::setprecision(8) << best_tree;
    std::cout << std::setprecision(8) << best_tree;
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

bool Inference::apply_swap(const vector<vector<int>> &D, const vector<int> &r, bool weighted, bool test_mode) {

    vector<Node*> swapped_nodes;
    try {
        swapped_nodes = t_prime.swap_labels(weighted, test_mode);
    }catch (const std::logic_error& e)
    {
        std::cout << " a logic error was caught during the swap labels move, with message '"
                  << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
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

void Inference::compute_t_prime_sums(const vector<vector<int>> &D) {

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
        vector<double> new_vals;

        for (auto &u_map : t_prime_scores[i]) {
            if (t_scores[i].count(u_map.first)) // add only if it is existing in the old vals // for the insertion case
                old_vals.push_back(t_scores[i][u_map.first]); // again the indices should match
            new_vals.push_back(u_map.second);
        }

        if (deleted_index != -1)
            old_vals.push_back(t_scores[i][deleted_index]);

        double res = MathOp::log_replace_sum(t_sums[i], old_vals, new_vals); // it takes t_sums[i]
        // subtracts the olds and adds the news
        // in case of delete, subtract an extra value
        // in case of insert, add an extra value
        // if the tree size changes, update it (tree.n_nodes). Posterior takes that into account
        assert(!isnan(res));

        t_prime_sums.push_back(res);
        i++;
    }
}

bool Inference::apply_add_remove_events(double lambda_r, double lambda_c, const vector<vector<int>> &D,
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
        std::cout << " a logic error was caught during the add remove events move, with message '"
                  << e.what() << "'\n";
        return false; // reject the move
    }
    catch (const std::exception& e) { // caught by reference to base
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

bool Inference::apply_insert_delete_node(double lambda_r, double lambda_c, const vector<vector<int>> &D,
                                         const vector<int> &r, bool weighted, bool validation_test_mode) {
    /*
     * Applies the insert/delete move on t_prime
     * Updates the sums and scores tables partially
     * */

    Node* tobe_computed;
    try {
        tobe_computed = t_prime.insert_delete_node(lambda_r, lambda_c, weighted, validation_test_mode);
    }catch (const std::out_of_range& e)
    {
        std::cout << " an out of range error was caught during the insert/delete node move, with message '"
                  << e.what() << "'\n";
        return false; // reject the move
    }catch (const std::exception& e) {
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


#endif //SC_DNA_INFERENCE_H
