//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <numeric>

struct Node{
    std::unordered_map<u_int,int> c;
    std::unordered_map<u_int,int> c_change;
    double log_score = 0.0;
    Node* first_child = nullptr;
    Node* next = nullptr;

    ~Node()
    {
        std::cout<<"node struct destructor"<<std::endl;
    };
};

class Tree {

public:
    Node* root;
    std::vector<Node*> all_nodes; // for uniform selection
    u_int ploidy; // to be added to values of the unordered map for each node

    bool is_leaf(Node*) const;
    Node* uniform_select();
    void random_insert(std::unordered_map<u_int, int>&&);
    void insert_at(u_int pos, std::unordered_map<u_int, int>&&);
    void insert_child(Node*, std::unordered_map<u_int, int>&&);
    Tree(u_int ploidy);



    virtual ~Tree();
    void traverse_tree();
    void destroy();
    void print_node(Node&);

    template<int N>
    void compute_root_score(int (&D)[N], int (&r)[N]);

    template<int N>
    void compute_score(Node* node, int (&D)[N], int& sum_D, int (&r)[N], std::unordered_map<u_int,int>& c_prev, double& log_parent);
private:
    void traverse(Node*);
    void update_label(std::unordered_map<u_int,int>& c_parent, Node* node);

    template<int N>
    int compute_z(int (&r)[N], std::unordered_map<u_int, int> &c);



};


template<int N>
void Tree::compute_root_score(int (&D)[N], int (&r)[N]) {

    int sum_d;
    sum_d = std::accumulate(D, D + std::size(D), 0);
    int z = 0;

    for (auto const &x : r)
        z += x * this->ploidy;

    root->log_score = sum_d * log(this->ploidy) - sum_d * log(z);
}

template<int N>
void Tree::compute_score(Node* node, int (&D)[N], int& sum_D, int (&r)[N], std::unordered_map<u_int,int>& c_prev, double& log_parent) {
    double val= 0.0;

    val += log_parent;

    for (auto const &x : node->c_change)
    {
        int cf = node->c[x.first];
        val += D[x.first] * (log(cf+ploidy));

        int cp_f = c_prev[x.first];
        val -= D[x.first] * (log(cp_f+ploidy));
    }
    int z = compute_z(r, node->c);
    int z_prev = compute_z(r, c_prev);

    val -= sum_D*log(z);
    val += sum_D*log(z_prev);

    node->log_score = val;

}

template<int N>
int Tree::compute_z(int (&r)[N], std::unordered_map<u_int, int> &c)
{
    int z = 0;

    for (int i = 0; i < std::size(r); ++i) {
        z += r[i] * (c[i]+ploidy);
    }

    return z;
}


//void Tree::traverse(Node* node) {
//    if (is_leaf(node))
//        std::cout<<"leaf: ";
//    else
//        std::cout<<"internal: ";
//    print_node(*node);
//    //std::cout<<node->log_score<<std::endl;
//    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
//        traverse(temp);
//    }
//}


#endif //SC_DNA_TREE_H
