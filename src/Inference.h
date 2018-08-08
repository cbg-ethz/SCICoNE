//
// Created by Tuncel  Mustafa Anil on 7/17/18.
//

#ifndef SC_DNA_INFERENCE_H
#define SC_DNA_INFERENCE_H

#include "Tree.h"
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <string>

class Inference {
/*
 * Contains functionality to perform monte carlo markov chains (mcmc) inference
 * //TODO use stack dynamic Trees for t & t_prime
 * */
public:
    Tree *t;
    Tree *t_prime;
    std::vector<std::vector<double>> t_scores;
    std::vector<double> t_sums;
    std::vector<unordered_map<int, double>> t_prime_scores;
    std::vector<double> t_prime_sums;
    std::string f_name;

public:
    Inference(u_int ploidy=2);
    ~Inference();
    void destroy();
    void compute_t_table(const vector<vector<int>> &D, const vector<int>& r);
    void compute_t_prime_scores(Node *attached_node, const vector<vector<int>> &D, const vector<int> &r);
    Node * apply_prune_reattach(const vector<vector<int>> &D, const vector<int> &r, bool weighted=false);
    void apply_swap(bool weighted=false);
    bool comparison(int m);

    void infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs);
    void write_best_tree();
    void update_t_scores();





    Tree *get_t() const;

    Tree *get_t_prime() const;
    void random_initialize(); // randomly initializes a tree and copies it into the other
    void initialize_example(); // initializes the trees based on the test example
};



void Inference::random_initialize() {

    t->random_insert({{0, 1}, {1, 1}});
    t->random_insert({{1, 1}, {2, 1}});
    t->random_insert({{0, -1}});
    t->random_insert({{3, -1}});
    t->random_insert({{1, 1}});

}

void Inference::initialize_example() {

    // build tree
    // tree that generated the data
    t->random_insert({{0, 1}, {1, 1}});
    t->insert_at(1,{{1, 1}, {2, 1}});
    t->insert_at(2,{{0, -1}});
    t->insert_at(2,{{3, -1}});
    t->insert_at(1,{{1, 1}});


    // 1532098496221375.txt
//    t->random_insert({{0, 1}, {1,1}});
//    t->insert_at(0,{{1,1},{2,1}});
//    t->insert_at(1,{{3,-1}});
//    t->insert_at(3,{{1,1}});
//    t->insert_at(0,{{0,-1}});


}

Inference::Inference(u_int ploidy) {
    t = new Tree(ploidy);
    t_prime =new Tree(ploidy);

    std::ofstream outfile;
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count(); // get a seed from time
    f_name = std::to_string(seed) + ".txt";


}

Inference::~Inference() {
    destroy();
}

Node* Inference::apply_prune_reattach(const vector<vector<int>> &D, const vector<int> &r, bool weighted) {
    /*
     * Applies prune and reattach to t_prime
     * Updates the sums and scores tables partially
     * */

    using namespace std;
    // weighted = false
    Node* attached_node = t_prime->prune_reattach(weighted);

    if (attached_node != nullptr)
    {
        compute_t_prime_scores(attached_node, D, r);

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
        return attached_node;
    }
    else
        return nullptr;
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

bool Inference::comparison(int m) {
    // m is

    double log_post_t = 0.0;
    double log_post_t_prime = 0.0;

    double t_sum = accumulate( t_sums.begin(), t_sums.end(), 0.0);
    int t_n = t->get_n_nodes();
    log_post_t = t_sum - (t_n -1 + m ) * log(t_n+1);

    double t_prime_sum = accumulate( t_prime_sums.begin(), t_prime_sums.end(), 0.0);
    int tp_n = t_prime->get_n_nodes();
    log_post_t_prime = t_prime_sum - (tp_n -1 + m ) * log(tp_n+1);

    double acceptance_prob = exp(log_post_t_prime - log_post_t);

    cout<<"acceptance prob: "<<acceptance_prob<<endl;
    std::ofstream outfile;
    outfile.open(f_name, std::ios_base::app);

    if (acceptance_prob > 1)
    {
        outfile << log_post_t_prime << ',';
        return true;
    }

    else
    {
        long long int seed = std::chrono::system_clock::now().time_since_epoch().count(); // get a seed from time
        static std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        double rand_val = distribution(generator);

        cout<<"rand_val: "<<rand_val<<endl;

        if (acceptance_prob > rand_val)
        {
            outfile << log_post_t_prime << ',';
            return true;
        }

        else
        {
            outfile << log_post_t << ',';
            return false;
        }

    }
}

void Inference::infer_mcmc(const vector<vector<int>> &D, const vector<int> &r, const vector<float> &move_probs) {

    int m = static_cast<int>(D.size());
    int n_accepted = 0;
    int n_rejected = 0;
    int n_attached_to_the_same_pos = 0;
    for (int i = 0; i < 500; ++i) {



        long long int seed = std::chrono::system_clock::now().time_since_epoch().count(); // get a seed from time
        std::default_random_engine generator(seed);
        std::discrete_distribution<> d(move_probs.begin(), move_probs.end());

        unsigned move_id = d(generator);

        switch (move_id)
        {
            case 0:
            {
                // prune & reattach
                cout << "Prune and reattach" << endl;
                Node *node = apply_prune_reattach(D, r, false);
                if (node == nullptr) {
                    n_attached_to_the_same_pos++;
                    continue;
                }
                break;
            }
            case 1:
            {
                // weighted prune & reattach
                cout<<"Weighted prune and reattach"<<endl;
                Node* node = apply_prune_reattach(D, r, true); // weighted=true
                if (node == nullptr)
                {
                    n_attached_to_the_same_pos++;
                    continue;
                }
                break;
            }
            case 2:
                // swap labels
                cout<<"swap labels"<<endl;
            case 3:
                // weighted swap labels
                cout<<"weighted swap labels"<<endl;
                break;
            default:
                throw std::logic_error("undefined move index");
        }

        // compare the trees
        bool accepted = comparison(m);
        // update trees and the matrices
        if (accepted)
        {
            n_accepted++;

            t_sums = t_prime_sums;
            update_t_scores();
            *t = *t_prime;
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
    outfile << *t;
}

void Inference::compute_t_prime_scores(Node *attached_node, const vector<vector<int>> &D, const vector<int> &r) {

    int j = 0;
    for (auto const &d: D)
    {
        int sum_d = accumulate( d.begin(), d.end(), 0.0);
        attached_node->parent->log_score = t_scores[j][attached_node->parent->id];
        t_prime->compute_stack(attached_node, d, sum_d,r);
        t_prime_scores.push_back(t_prime->get_children_id_score(attached_node));
        j++;
    }
}

void Inference::apply_swap(bool weighted) {

    vector<Node*> swapped_nodes = t_prime->swap_labels(weighted);



}


#endif //SC_DNA_INFERENCE_H
