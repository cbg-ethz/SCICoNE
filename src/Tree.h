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
    std::unordered_map<u_int,int> c = {};
    std::unordered_map<u_int,int> c_change= {};
    double log_score = 0.0;
    int z = 0;
    Node* first_child = nullptr;
    Node* next = nullptr;
    Node* parent = nullptr;

    // copy constructor
    Node(Node& source_node)
    {
        c = source_node.c;
        c_change = source_node.c_change;
        log_score = source_node.log_score;
        z = source_node.z;
        first_child = next = parent = nullptr;
    }
    Node()
    {}
    ~Node()
    {}
};

class Tree {

public:
    Node* root;
    std::vector<Node*> all_nodes; // for uniform selection
    u_int ploidy; // to be added to values of the unordered map for each node

    // constructor
    Tree(u_int ploidy);
    // copy constructor
    Tree(Tree& source, u_int ploidy);
    // destructor
    virtual ~Tree();

    void copy_tree(const Tree& source_tree);
    void recursive_copy(Node *source, Node *destination);

    bool is_leaf(Node*) const;
    Node* uniform_select(bool with_root);
    void random_insert(std::unordered_map<u_int, int>&&);
    void insert_at(u_int pos, std::unordered_map<u_int, int>&&);
    void insert_child(Node *pos, std::unordered_map<u_int, int>&& labels);
    Node* insert_child(Node *pos, Node source);
    double average_score();
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

//    alternative recursive implementation
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

Node* Tree::uniform_select(bool with_root=true) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    int rand_val = 0;
    if (all_nodes.size() == 0)
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes.size() ==1)
        rand_val = 1;
    else
    {
        // 1 is root, 2 is 1st
        std::uniform_int_distribution<> dis(with_root?1:2, all_nodes.size());
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


Node* Tree::insert_child(Node *pos, Node source) {
    /*
     * Creates a child node from the reference node and inserts it to the position pos
     * */

    Node *child = new Node(source);
    child->parent = pos;
    all_nodes.push_back(child);


    // insert
    if (is_leaf(pos))
    {
        pos->first_child = child;
    }
    else
    {
        // find the last child
        for (Node* temp = pos->first_child; temp != nullptr; temp=temp->next) {
            if (temp->next != nullptr)
                continue;
            else // add to the last child
                temp->next = child;
            break;
        }
    }
    return child;

}

void Tree::insert_child(Node* pos, std::unordered_map<u_int, int>&& labels) {

    /*
     * Creates a child node from the labels and inserts it to the position pos
     * */

    // create node
    Node* child = new Node();

    // move operator
    child->c_change = labels;

    // copy (not computed yet)
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
    std::cout<<"destroyed."<<std::endl;
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
    /*
     * Updates the child label (c) based on the parent c
     * */

    // copy assignment
    node->c = c_parent;
    // parent labels are passed to the child c here
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

Tree::Tree(Tree &source, u_int ploidy) {
    if (source.root == nullptr)
        this->root = nullptr;
    else
        this->copy_tree(source);
}

void Tree::copy_tree(const Tree& source_tree) {
    /*
     * TODO
     * have a single copy node method (done)
     * call it recursively (using stack)
     * copy all the other attributes
    */



    this->ploidy = source_tree.ploidy;

    // copy the nodes
    this->root = new Node(*source_tree.root);
    this->all_nodes.push_back(root);
    recursive_copy(source_tree.root, this->root);


//    // stack based implementation
//    std::stack<Node*> stk;
//    stk.push(this->root);
//
//    while (!stk.empty()) {
//        Node* top = (Node*) stk.top();
//        stk.pop();
//        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
//            stk.push(temp);
//        }
//
//        compute_score(top, D, sum_D, r);
//    }


 }

void Tree::recursive_copy(Node* source, Node *destination) {
    /*
     * Copies the tree from source to destination
     * TODO implement it using stack & loop, not to use the function call stack
     * TODO the order of insertion is different (although the tree is the same)
     * */

    for (Node* temp = source->first_child; temp != nullptr; temp=temp->next) {
        auto child = insert_child(destination, *temp);
        recursive_copy(temp, child);
    }


}




#endif //SC_DNA_TREE_H
