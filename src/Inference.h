//
// Created by Tuncel  Mustafa Anil on 7/17/18.
//

#ifndef SC_DNA_INFERENCE_H
#define SC_DNA_INFERENCE_H

#include "Tree.h"
#include <vector>

class Inference {
/*
 * Contains functionality to perform monte carlo markov chains (mcmc) inference
 * */
private:
    Tree *t;
    Tree *t_prime;
    std::vector<std::vector<double>> t_scores;

public:
    Inference(u_int ploidy=2);

    virtual ~Inference();

    template<int N>
    void compute_t_table(int (&D)[N], int (&r)[N]);

    void prune_reattach();
    void w_prune_reattach();

    void swap();
    void w_swap();


    void random_initialize(); // randomly initializes a tree and copies it into the other
    void test_initialize(); // initializes the trees based on the test example
};



void Inference::random_initialize() {

}

void Inference::test_initialize() {

    // build tree
    t->random_insert({{0, 1}, {1, 1}});
    t->insert_at(1,{{1, 1}, {2, 1}});
    t->insert_at(2,{{0, -1}});
    t->insert_at(2,{{3, -1}});
    t->insert_at(1,{{1, 1}});
}

Inference::Inference(u_int ploidy) {
    t = new Tree(ploidy);
    t_prime = new Tree(ploidy);


}

Inference::~Inference() {
    t->destroy();
    t_prime->destroy();
}

void Inference::prune_reattach() {
    auto attached_node = t_prime->prune_reattach();


}

template<int N>
void Inference::compute_t_table(int (&D)[N], int (&r)[N]) {

//    int n = std::size(D);
//    for (int i = 0; i < n; ++i) {
//        this->t->compute_tree(D[i], r);
//        sum_scores.push_back(t.sum_score());
//    }

};

#endif //SC_DNA_INFERENCE_H
