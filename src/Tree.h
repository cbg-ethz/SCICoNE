//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <random>
#include <numeric>
#include <limits>

struct Node{
    std::unordered_map<u_int,int> c;
    std::unordered_map<u_int,int> c_change;
    double log_score = 0.0;
    int z = 0;
    Node* first_child = nullptr;
    Node* next = nullptr;
    Node* parent = nullptr;

    Node() = default;
    ~Node() = default;
    /* TODO: have a move constructor to use on MCMC moves!
     * the move constructor should preserve first_child, next, c_change
     * the log score and c will be re-computed for each different tree
     */
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
    void compute_score(Node* node, int (&D)[N], int& sum_D, int (&r)[N]);

private:

    void traverse(Node*);
    void update_label(std::unordered_map<u_int,int>& c_parent, Node* node);

    template<int N>
    int compute_z(int (&r)[N], std::unordered_map<u_int, int> &c);

    template<int N>
    void compute_stack(Node *node, int (&D)[N], int &sum_D, int (&r)[N]);


};


template<int N>
void Tree::compute_root_score(int (&D)[N], int (&r)[N]) {

    int sum_d;
    sum_d = std::accumulate(D, D + std::size(D), 0);
    int z = 0;

    for (auto const &x : r)
        z += x * this->ploidy;

    root->log_score = sum_d * log(this->ploidy) - sum_d * log(z);
    root->z = z;
}

template<int N>
void Tree::compute_score(Node* node, int (&D)[N], int& sum_D, int (&r)[N]) {


    double val = node->parent->log_score;
    int z = node->parent->z;

    for (auto const &x : node->c_change)
    {

        int cf = node->c[x.first];
        val += D[x.first] * (log(cf+ploidy));



        int cp_f = node->parent->c[x.first];
        val -= D[x.first] * (log(cp_f+ploidy));

        z += r[x.first] * (cf - cp_f);

    }


    val -= sum_D*log(z);
    val += sum_D*log(node->parent->z);

    node->log_score = val;
    node->z = z;
}


template<int N>
void Tree::compute_tree(int (&D)[N], int (&r)[N]) {

    // for the root
    compute_root_score(D,r);
    double root_score = root->log_score;
    int root_z = root->z;
    //reuse the computed sum in each node
    int sum_d = std::accumulate(D, D + std::size(D), 0);
    compute_stack(root->first_child, D, sum_d, r);

}

template<int N>
void Tree::compute_stack(Node *node, int (&D)[N], int &sum_D, int (&r)[N])
{
    // stack based implementation
    std::stack<Node*> stk;
    stk.push(node);

    while (!stk.empty()) {
        Node* top = (Node*) stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        compute_score(top, D, sum_D, r);
    }

    // alternative recursive implementation
//    compute_score(node, D, sum_D, r);
//    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
//        compute_stack(temp, D, sum_D, r);
//    }
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

    // copy
    child->c = child->c_change;

    child->parent = pos;

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
