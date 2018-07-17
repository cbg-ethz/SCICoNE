//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <numeric>
#include <limits>
#include "Node.h"
#include "MathOp.h"


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

    // TODO overload the assignment operator!


    // moves
    Node * prune_reattach();

    bool is_leaf(Node*) const;
    Node* uniform_select(bool with_root);
    void random_insert(std::unordered_map<u_int, int>&&);
    void insert_at(u_int pos, std::unordered_map<u_int, int>&&);
    void insert_child(Node *pos, std::unordered_map<u_int, int>&& labels);
    Node * insert_child(Node *pos, Node *source);
    Node* insert_child(Node *pos, Node& source);
    std::vector<double> get_scores();
    Tree& operator=(const Tree& other);


    double sum_score();
    void traverse_tree();
    void destroy();
    void print_node(Node&);
    template<int N>
    void compute_tree(const vector<int> &D, int (&r)[N]);
    template<int N>
    void compute_root_score(const vector<int> &D, int& sum_d, int (&r)[N]);

    template<int N>
    void compute_stack(Node *node, const vector<int> &D, int &sum_D, int (&r)[N]);
    u_int get_n_nodes() const;
    int counter = 0;

private:


    u_int n_nodes; //the number of nodes without the root
    void traverse(Node*);
    void update_label(std::unordered_map<u_int,int>& c_parent, Node* node);
    void copy_tree(const Tree& source_tree);
    void recursive_copy(Node *source, Node *destination);
    template<int N>
    void compute_score(Node* node, const vector<int> &D, int& sum_D, int (&r)[N]);
    Node* remove(Node* pos); // does not deallocate, TODO: have a delete method that calls this and deallocates

};


template<int N>
void Tree::compute_root_score(const vector<int> &D, int& sum_d, int (&r)[N]) {

    int z = 0;

    for (auto const &x : r)
        z += x * this->ploidy;

    root->log_score = sum_d * log(this->ploidy) - sum_d * log(z);
    root->z = z;
}

template<int N>
void Tree::compute_score(Node* node, const vector<int> &D, int& sum_D, int (&r)[N]) {


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
void Tree::compute_tree(const vector<int> &D, int (&r)[N]) {


    //reuse the computed sum in each node
    int sum_d = std::accumulate(D.begin(), D.end(), 0);

    // for the root
    compute_root_score(D, sum_d,r);
    double root_score = root->log_score;
    int root_z = root->z;
    compute_stack(root->first_child, D, sum_d, r);

}

template<int N>
void Tree::compute_stack(Node *node, const vector<int> &D, int &sum_D, int (&r)[N])
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
    n_nodes = 0;
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

    int rand_val = 0;
    if (all_nodes.size() == 0)
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes.size() ==1)
        rand_val = 1;
    else
    {
        // 1 is root, 2 is 1st.
        rand_val = MathOp::random_uniform(with_root?1:2,all_nodes.size());
    }

    return all_nodes[rand_val-1];
}

bool Tree::is_leaf(Node* n) const{
    if (n->first_child == nullptr)
        return true;
    else
        return false;
}

Node * Tree::insert_child(Node *pos, Node *source) {

    // set the parent from the child
    source->parent = pos;
    // set the unique id
    if (source->id == 0)
        source->id = ++counter;
    // if the id is set then keep it

    all_nodes.push_back(source);
    n_nodes++;

    // insert
    if (is_leaf(pos))
    {
        pos->first_child = source;
        update_label(pos->c, source);
    }
    else
    {
        // find the last child
        for (Node* temp = pos->first_child; temp != nullptr; temp=temp->next) {
            if (temp->next != nullptr)
                continue;
            else // add to the last child
                temp->next = source;
            update_label(pos->c, source);
            break;
        }
    }
    return source;

}

Node* Tree::insert_child(Node *pos, Node& source) {
    /*
     * Creates a child node from the reference node and inserts it to the position pos
     * Does not preserve the link informations!
     * This is used for copying purposes.
     * */

    // struct copy constructor
    Node *child = new Node(source);
    return insert_child(pos,child);

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
    insert_child(pos,child);

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
    Node* pos = uniform_select(true);
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


double Tree::sum_score() {

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

    return log(sum_in_normal_space) + max;




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

    // copy the nodes using struct copy constructor
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
     * use either 2 stacks or a stack and a tuple
     * */

    for (Node* temp = source->first_child; temp != nullptr; temp=temp->next) {
        auto child = insert_child(destination, *temp);
        recursive_copy(temp, child);
    }


}

Node * Tree::prune_reattach() {
    // returns the pruned node (which is the attached node)

    Node* prune_pos;
    prune_pos = uniform_select(false);

    // copy all nodes
    std::vector<Node*> destination_nodes = this->all_nodes;

    // TODO remove all the descendents of the prune_pos

    std::stack<Node*> stk;
    stk.push(prune_pos);

    while (!stk.empty()) {
        Node *top = (Node *) stk.top();
        stk.pop();
        for (Node *temp = top->first_child; temp != nullptr; temp = temp->next) {
            stk.push(temp);
        }
        destination_nodes.erase(std::remove(destination_nodes.begin(), destination_nodes.end(), top), destination_nodes.end());
    }

    int rand_val = 0;
    rand_val = MathOp::random_uniform(1,destination_nodes.size());
    Node* attach_pos = nullptr;
    attach_pos = destination_nodes[rand_val -1];
    // attach_pos = all_nodes[5]; //TODO remove this after testing is complete
    //do not recompute you attach at the same pos
    if (prune_pos->parent != attach_pos)
    {
        /*
         * 1. Remove prune_pos (done)
         * 2. Insert prune_pos using insert_child function
         * */
        auto pruned_node = remove(prune_pos);
        auto attached_node = insert_child(attach_pos, pruned_node);

        // std sort the all_nodes vector to make sure the indices match between 2 trees
                // e.g. node id 4 will always be at index 4
        std::sort(this->all_nodes.begin(),this->all_nodes.end());
        return attached_node;
    }
    return nullptr;



}

Node *Tree::remove(Node *pos) {

    /*
     * Removes the element at the position but does not deallocate!
     * The caller method has the responsibility to deallocate.
     * */

    assert(pos != root);

    // if pos is the first child then pos' next is the new first child
    if (pos == pos->parent->first_child)
        pos->parent->first_child = pos->next;
    // otherwise find the node among the childs and connect it's previous to its next
    else
    {
        for (Node *temp = pos->parent->first_child; temp != nullptr; temp = temp->next)
        {
            if (temp->next == pos)
            {
                temp->next = temp->next->next;
                break;
            }

        }
    }

    // remove from the all_nodes as well
    all_nodes.erase(std::remove(all_nodes.begin(), all_nodes.end(), pos), all_nodes.end());

    //remove the next link
    pos->next = nullptr;

    n_nodes--;
    return pos;

}

u_int Tree::get_n_nodes() const {
    return n_nodes;
}

Tree &Tree::operator=(const Tree &other) {

    if (this->root != nullptr)
        this->destroy();

    if (other.root == nullptr)
        this->root = nullptr;
    else
        this->copy_tree(other);


    return *this;
}

std::vector<double> Tree::get_scores() {

    vector<double> scores;

    for (auto const &i : all_nodes)
        scores.push_back(i->log_score);
    return scores;
}


#endif //SC_DNA_TREE_H
