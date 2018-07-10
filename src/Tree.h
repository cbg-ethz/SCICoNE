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
#include <limits>

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

    double average_score();


    virtual ~Tree();
    void traverse_tree();
    void destroy();
    void print_node(Node&);

    template<int N>
    void compute_tree(int (&D)[N], int (&r)[N]);



    template<int N>
    void compute_root_score(int (&D)[N], int (&r)[N]);

    template<int N>
    void compute_score(Node* node, int (&D)[N], int& sum_D, int (&r)[N], std::unordered_map<u_int,int>& c_prev, double& log_parent);

private:

    void traverse(Node*);
    void update_label(std::unordered_map<u_int,int>& c_parent, Node* node);

    template<int N>
    int compute_z(int (&r)[N], std::unordered_map<u_int, int> &c);

    template<int N>
    void tail_compute(Node *node, int (&D)[N], int &sum_D, int (&r)[N], std::unordered_map<u_int, int> &c_prev,
                      double &log_parent);


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

template<int N>
void Tree::compute_tree(int (&D)[N], int (&r)[N]) {

    // for the root
    compute_root_score(D,r);
    double root_score = root->log_score;
    //reuse the computed sum in each node
    int sum_d = std::accumulate(D, D + std::size(D), 0);
    tail_compute(root->first_child, D, sum_d, r, root->c, root_score);

}

template<int N>
void Tree::tail_compute(Node *node, int (&D)[N], int &sum_D, int (&r)[N], std::unordered_map<u_int, int> &c_prev,
                        double &log_parent)
{
    compute_score(node, D, sum_D, r, c_prev, log_parent);
    //std::cout<<node->log_score<<std::endl;
    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
        tail_compute(temp, D, sum_D, r, node->c, node->log_score);
    }
}


Tree::Tree(u_int ploidy)
{
    root = new Node();
    this->ploidy = ploidy;
    // creates a copy of the root ptr and stores it in the vector

    all_nodes.push_back(root);


}

Tree::~Tree() {
    destroy();
}

void Tree::traverse_tree() {
    traverse(root);
}

void Tree::traverse(Node* node) {
    if (is_leaf(node))
        std::cout<<"leaf: ";
    else
        std::cout<<"internal: ";
    print_node(*node);
    //std::cout<<node->log_score<<std::endl;
    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
        traverse(temp);
    }
}

Node* Tree::uniform_select() {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    int rand_val = 0;
    if (all_nodes.size() == 0)
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes.size() ==1)
        rand_val = 1;
    else
    {
        std::uniform_int_distribution<> dis(1, all_nodes.size());
        rand_val = dis(gen);
    }

    return all_nodes[rand_val-1];
}

bool Tree::is_leaf(Node* n) const{
    if (n->first_child == nullptr)
        return true;
    else
        return false;
}


void Tree::insert_child(Node* pos, std::unordered_map<u_int, int>&& labels) {

    // create node
    Node* child = new Node();
    // move operator

    child->c_change = labels;
    child->c = child->c_change;


    all_nodes.push_back(child);



    // insert
    if (is_leaf(pos))
    {
        pos->first_child = child;
        update_label(pos->c, child);
    }
    else
    {
        // find the last child
        for (Node* temp = pos->first_child; temp != nullptr; temp=temp->next) {
            if (temp->next != nullptr)
                continue;
            else // add to the last child
                temp->next = child;
            update_label(pos->c, child);
            break;

        }

    }

}

void Tree::destroy() {
    for (auto elem: all_nodes)
    {
        std::cout<<"deleting " << elem->log_score <<std::endl;
        elem->first_child = elem->next = nullptr;
        elem->c.clear();
        delete elem;
        elem = nullptr;
    }
    all_nodes.clear();
    root = nullptr;
}

void Tree::random_insert(std::unordered_map<u_int, int>&& labels)
{
    Node* pos = uniform_select();
    insert_child(pos, std::move(labels));

}

void Tree::insert_at(u_int pos, std::unordered_map<u_int, int> && labels) {
    Node* n = all_nodes[pos];
    insert_child(n, std::move(labels));

}

void Tree::print_node(Node& n) {
    std::cout<<"node.";
    std::cout<<std::endl;
    for (auto i : n.c_change)
        std::cout << " " << i.first << ":" << i.second << std::endl;
}


void Tree::update_label(std::unordered_map<u_int, int>& c_parent, Node *node) {
    // copy assignment
    node->c = c_parent;
    for (auto it=node->c_change.begin(); it!=node->c_change.end(); ++it) {
        node->c[it->first] = node->c[it->first] + it->second;
    }


}


double Tree::average_score() {

    // init max with the smallest value possible
    double max = std::numeric_limits<double>::min();

    for (auto &v2: all_nodes)
        if (v2->log_score > max)
            max = v2->log_score;

    double scores[all_nodes.size()];

    double sum_in_normal_space = 0.0;
    for (int i = 0; i < all_nodes.size(); ++i) {
        scores[i] = all_nodes[i]->log_score - max;
        scores[i] = exp(scores[i]);
        sum_in_normal_space += scores[i];
    }

    return log(sum_in_normal_space / all_nodes.size()) + max;




}


#endif //SC_DNA_TREE_H
