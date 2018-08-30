//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <map>
#include <unordered_map>
#include <numeric>
#include <limits>
#include "Node.h"
#include "MathOp.h"
#include <tuple>
#include <random>
#include "SingletonRandomGenerator.h"
#include <set>
#include <iomanip>
#include "Utils.h"

#include <algorithm> // std::remove


class Tree {

public:
    Node* root;
    std::vector<Node*> all_nodes_vec; // for random selection, destructor, insertion by position, iterating without order (e.g. for each node)
    int ploidy; // to be added to values of the unordered map for each node
    u_int n_regions; // number of regions
    double score; // log posterior score of the tree
    int counter = 0;
    u_int n_nodes; //the number of nodes without the root

    // constructor
    Tree(u_int ploidy, u_int n_regions);
    // copy constructor
    Tree(Tree& source);
    // destructor
    virtual ~Tree();

    // moves
    Node* prune_reattach(bool weighted=false, bool validation_test_mode=false);
    std::vector<Node*> swap_labels(bool weighted=false, bool validation_test_mode=false);
    Node* add_remove_events(float lambda_r, float lambda_c, bool weighted=false, bool validation_test_mode=false);
    Node* delete_node();

    bool is_leaf(Node*) const;
    Node* uniform_sample(bool with_root=true);
    Node* weighted_sample();
    void random_insert(std::map<u_int, int>&&);
    void insert_at(u_int pos, std::map<u_int, int>&&);
    void insert_child(Node *pos, std::map<u_int, int>&& labels);

    map<int, double> get_children_id_score(Node *node);
    void destroy();
    void compute_tree(const vector<int> &D, const vector<int>& r);
    void compute_root_score(const vector<int> &D, int& sum_d, const vector<int>& r);
    void compute_stack(Node *node, const vector<int> &D, int &sum_D, const vector<int>& r);

    void compute_weights();
    u_int get_n_nodes() const;
    vector<Node*> get_descendents(Node* n);

    friend std::ostream& operator<<(std::ostream& os, Tree& t);
    Tree& operator=(const Tree& other);

private:


    void update_label(std::map<u_int,int>& c_parent, Node* node);
    void update_desc_labels(Node* node);

    //validation of tree
    bool is_redundant();

    // Validation of subtrees
    bool is_valid_subtree(Node* node);
    bool subtree_contains_negative(Node* n);
    bool zero_ploidy_changes(Node* n);
    bool region_changes(Node *n, u_int region_id);


    void copy_tree(const Tree& source_tree);
    void recursive_copy(Node *source, Node *destination);
    void compute_score(Node *node, const vector<int> &D, int &sum_D, const vector<int> &r, float eta=0.0001f);
    Node* prune(Node *pos); // does not deallocate,
    Node* insert_child(Node *pos, Node *source);
    Node* insert_child(Node *pos, Node& source);
    bool is_ancestor(Node *target, Node *curr);

    bool empty_hashmap(map<u_int, int> &dict);

};

std::ostream& operator<<(std::ostream& os, Tree& t) {

    vector<Node*> nodes = t.get_descendents(t.root);

    os << "Tree score: " << setprecision(8) << t.score << endl;
    for (auto const &x : nodes)
    {
        os << *x << std::endl;
    }
    return os;
}



void Tree::compute_root_score(const vector<int> &D, int& sum_d, const vector<int>& r) {

    int z = 0;
    for (auto const &x : r)
        z += x * this->ploidy;

    root->log_score = sum_d * log(this->ploidy) - sum_d * log(z);
    root->z = z;
}

void Tree::compute_score(Node *node, const vector<int> &D, int &sum_D, const vector<int> &r, float eta) {


    if (node->parent == nullptr)
    {
        compute_root_score(D, sum_D,r);
    }
    else

    {
        double val = node->parent->log_score;
        int z = node->parent->z;

        for (auto const &x : node->c_change)
        {

            // if log zero then use eta value, not to have -infinity
            int cf = (node->c.count(x.first)?node->c[x.first]:0); // use count to check without initializing the not found element
            // the above part can also be done by using map::at and exception handling

            val += D[x.first] * (log((cf+ploidy)==0?(eta):(cf+ploidy)));

            int cp_f = (node->parent->c.count(x.first) ?node->parent->c[x.first] : 0); // use count to check without initializing

            val -= D[x.first] * (log((cp_f+ploidy)==0?(eta):(cp_f+ploidy)));

            z += r[x.first] * (cf - cp_f);

        }

        val -= sum_D*log(z);
        val += sum_D*log(node->parent->z);

        node->log_score = val;
        node->z = z;
    }


}


void Tree::compute_tree(const vector<int> &D, const vector<int>& r) {

    //reuse the computed sum in each node
    int sum_d = std::accumulate(D.begin(), D.end(), 0);
    compute_stack(root, D, sum_d, r);

}

void Tree::compute_stack(Node *node, const vector<int> &D, int &sum_D, const vector<int>& r)
{
    /*
     * Computes the nodes scores in a top-down fashion
     * time complexity: O(n), where n is the number of nodes
     * */

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
        // TODO: use a function pointer to compute_score or traverse or update_desc_labels because rest of the code is the same.
    }

//    alternative recursive implementation
//    compute_score(node, D, sum_D, r);
//    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
//        compute_stack(temp, D, sum_D, r);
//    }
}


Tree::Tree(u_int ploidy, u_int n_regions)
{
    root = new Node();
    this->ploidy = ploidy;
    this->n_regions = n_regions;
    n_nodes = 0;
    score = 0.0;
    // creates a copy of the root ptr and stores it in the vector

    all_nodes_vec.push_back(root);


}

Tree::~Tree() {
    destroy();
}

Node* Tree::uniform_sample(bool with_root) {

    int rand_val = 0;
    if (all_nodes_vec.size() == 0)
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes_vec.size() ==1)
        rand_val = 1;
    else
    {
        // 1 is root, 2 is 1st.
        rand_val = MathOp::random_uniform(with_root?1:2,all_nodes_vec.size());
    }

    return all_nodes_vec[rand_val-1];
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

    all_nodes_vec.push_back(source);
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

void Tree::insert_child(Node* pos, std::map<u_int, int>&& labels) {

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
    for (auto elem: all_nodes_vec)
    {
        //std::cout<<"deleting " << elem->log_score <<std::endl;
        elem->first_child = elem->next = nullptr;
        elem->c.clear();
        delete elem;
        elem = nullptr;
    }
    all_nodes_vec.clear();
    n_nodes = 0;
    root = nullptr;
    //std::cout<<"destroyed."<<std::endl;
}

void Tree::random_insert(std::map<u_int, int>&& labels)
{
    Node* pos = uniform_sample(true);
    insert_child(pos, std::move(labels));

}

void Tree::insert_at(u_int pos, std::map<u_int, int> && labels) {
    Node* n = all_nodes_vec[pos];
    insert_child(n, std::move(labels));

}


void Tree::update_label(std::map<u_int, int>& c_parent, Node *node) {
    /*
     * Updates the child label (c) based on the parent c
     * */

    // copy assignment
    node->c = c_parent;
    // parent labels are passed to the child c here
    for (auto it=node->c_change.begin(); it!=node->c_change.end(); ++it) {
        int new_value = 0;
        new_value = node->c[it->first] + it->second;
        if (new_value == 0)
            node->c.erase(it->first); //erase the zero instead of storing it, if not found by default it is already zero in map.
        else
            node->c[it->first] = new_value;
    }


    // compute and store the hash
    vector<int> keys_values = {};
    for (auto const &it : node->c)  {
        keys_values.push_back(it.first);
        keys_values.push_back(it.second);
    };
    int size_for_hash = size(keys_values) * sizeof(keys_values[0]);
    uint64_t c_hash = Utils::calcul_hash(&keys_values[0], size_for_hash);

    node->c_hash = c_hash;
}



Tree::Tree(Tree &source) {
    if (source.root == nullptr)
        this->root = nullptr;
    else
        this->copy_tree(source);
}

void Tree::copy_tree(const Tree& source_tree) {
    /*
     * TODO
     * call it recursively (using stack)
    */

    this->ploidy = source_tree.ploidy;
    this->counter = source_tree.counter;
    this->score = source_tree.score;
    this->n_regions = source_tree.n_regions;

    // copy the nodes using struct copy constructor
    this->root = new Node(*source_tree.root);
    this->all_nodes_vec.push_back(root); // all nodes cannot be copied since 2 trees cannot have nodes pointing to the same address
    recursive_copy(source_tree.root, this->root); // the nodes of the source tree are inserted into the destination tree
    this->n_nodes = source_tree.n_nodes; // copy this after the insertions are done (insertions change this value).
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

Node * Tree::prune_reattach(bool weighted, bool validation_test_mode) {
    // returns the pruned node (which is the attached node)

    Node* prune_pos = nullptr;

    if (validation_test_mode)
        prune_pos = all_nodes_vec[2];
    else
    {
        if (weighted)
            prune_pos = weighted_sample();
        else
            prune_pos = uniform_sample(false); //without the root
    }


    // copy all nodes
    std::vector<Node*> destination_nodes = this->all_nodes_vec;

    // remove all the descendents of the prune_pos

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

    if (validation_test_mode)
        attach_pos = all_nodes_vec[5];
    else
        attach_pos = destination_nodes[rand_val -1];

    //do not recompute if you attach at the same pos
    if (prune_pos->parent->id != attach_pos->id)
    {
        /*
         * 1. Remove prune_pos (done)
         * 2. Insert prune_pos using insert_child function
         * */
        auto pruned_node = prune(prune_pos);
        auto attached_node = insert_child(attach_pos, pruned_node);

        //update the c vectors of the attached node and its descendents
        update_desc_labels(attached_node);

        if (!is_valid_subtree(attached_node) || is_redundant())
            return nullptr;

        // recompute the weights after the tree structure is changed
        this->compute_weights();

        return attached_node;
    }
    else
    {
        return prune_pos = attach_pos = nullptr;
    }


}

Node *Tree::prune(Node *pos) {

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

    // remove from the all_nodes_vec as well
    all_nodes_vec.erase(std::remove(all_nodes_vec.begin(), all_nodes_vec.end(), pos), all_nodes_vec.end());

    //prune the next link
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



map<int, double> Tree::get_children_id_score(Node *node) {
/*
 * Returns the ids and the log scores of the descendent nodes
 * */

    map<int,double> id_score_pairs;

    // stack based implementation
    std::stack<Node*> stk;
    stk.push(node);

    while (!stk.empty()) {
        Node* top = (Node*) stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        // make sure the id is not in the map before
        assert(id_score_pairs.find(top->id) == id_score_pairs.end());
        id_score_pairs[top->id] = top->log_score;
    }
    return id_score_pairs;
}

void Tree::compute_weights() {

    /*
     * a bottom-up approach
     * time complexity: O(n), where n is the number of nodes
     * */


    // stack to allow the bottom-up computations
    std::stack<Node*> cmp_stack;


    //stack for the tree traversal
    std::stack<Node*> stk;
    stk.push(root); //start with the root

    while (!stk.empty()) {
        Node* top = (Node*) stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        cmp_stack.push(top);
    }

    while (!cmp_stack.empty())
    {
        Node* top = (Node*) cmp_stack.top();
        cmp_stack.pop();
        unsigned n_descendents = 1;

        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next)
            n_descendents += temp->n_descendents;

        top->n_descendents = n_descendents;
    }
}

Node *Tree::weighted_sample() {

    if (all_nodes_vec.size() == 0)
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes_vec.size() ==1)
        throw std::length_error("there is only 1 element which is the root and root cannot be sampled");
    else
    {
        // get the subvector
        vector<Node*>::const_iterator first = all_nodes_vec.begin() + 1;
        vector<Node*>::const_iterator last = all_nodes_vec.end();
        vector<Node*> nodes_to_sample(first, last);

        vector<float> weights;
        for (auto const &x : nodes_to_sample)
        {
            float weight = (1.0f / x->n_descendents); // weights are inversely proportional to n_descendents
            weights.push_back(weight);
        }

        std::mt19937 &generator = SingletonRandomGenerator::get_generator();
        std::discrete_distribution<> d(weights.begin(), weights.end());

        unsigned sample = d(generator);

        return nodes_to_sample[sample];

    }


}

bool Tree::is_ancestor(Node *target, Node *curr) {
    /*
     * Returns true if the target node is an ancestor of the current node
     * Complexity: O(n) where n is the number of ancestors of the curr node
     * */
    Node* p = curr;
    while(p->parent != nullptr)
    {
        if (p == target)
            return true;
        p = p->parent;
    }

    return false;
}

std::vector<Node *> Tree::swap_labels(bool weighted, bool validation_test_mode) {

    /*
     * Swaps the labels between two nodes
     * Requires more than 2 nodes to work
     * */

    if (all_nodes_vec.size() <= 2)
        throw std::logic_error("swapping labels does not make sense when they is only one node besides the root");

    Node *node1, *node2;
    node1 = node2 = nullptr;

    if (weighted)
    {
        node1 = weighted_sample();

        do
            node2 = weighted_sample();
        while (node1 == node2); // make sure you are not swapping the same labels
    }
    else
    {
        if (validation_test_mode)
        {
            node1 = all_nodes_vec[2]; //without the root
            node2 = all_nodes_vec[5];
        }
        else
        {
            node1 = uniform_sample(false); //without the root
            do
                node2 = uniform_sample(false);
            while (node1 == node2);
        }

    }

    // perform std swap on maps
    node1->c_change.swap(node2->c_change);

    vector<Node*> return_nodes;

    // TODO only compute the nodes in between
    if (is_ancestor(node1, node2))
        return_nodes.push_back(node1);
    else if (is_ancestor(node2, node1))
        return_nodes.push_back(node2);
    else
    {
        return_nodes.push_back(node1);
        return_nodes.push_back(node2);
    }

    // updating the labels
    for (auto const &node : return_nodes)
    {
        this->update_desc_labels(node);

        if (!is_valid_subtree(node))
            return {}; // empty vector with list initialization

    }
    if (is_redundant())
        return {};


    return return_nodes;
}


bool Tree::empty_hashmap(map<u_int, int> &dict) {

    for (auto const &it : dict)
    {
        if(it.second != 0)
            return false;
    }
    return true;

}

void Tree::update_desc_labels(Node *node) {

    /*
     * Updates the c vectors of the node and all of its descendents
     * */


    std::stack<Node*> stk;
    stk.push(node);

    while (!stk.empty()) {
        Node* top = (Node*) stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        update_label(top->parent->c, top);
    }
}

bool Tree::is_valid_subtree(Node *node) {
    /*
     * Returns true if the subtree is valid
     * */

    if (subtree_contains_negative(node))
        return false;

    if (zero_ploidy_changes(node))
        return false;

    return true;

}

Node *Tree::add_remove_events(float lambda_r, float lambda_c, bool weighted, bool validation_test_mode) {

    Node* node;

    if (validation_test_mode)
    {
        node = all_nodes_vec[3];
        lambda_r = lambda_c = 0.0f;
    }
    else
    {
        if (weighted)
            node = weighted_sample();
        else
            node = uniform_sample(false); //without the root
    }
    std::mt19937 &generator = SingletonRandomGenerator::get_generator();


    // n_regions from Poisson(lambda_R)+1
    std::poisson_distribution<int> distribution(lambda_r); // the param is to be specified later
    int n_regions_to_sample = distribution(generator) + 1;
    // sample n_regions_to_sample distinct regions uniformly
    int n_regions = this->n_regions;
    int regions_sampled = 0;
    std::set<int> distinct_regions;

    // otherwise we cannot sample distinct uniform regions
    assert(n_regions_to_sample <= n_regions);

    while (regions_sampled < n_regions_to_sample)
    {
        int uniform_val = MathOp::random_uniform(0, n_regions-1);
        if (validation_test_mode)
            uniform_val =3;
        if (distinct_regions.find(uniform_val) == distinct_regions.end())
        {
            distinct_regions.insert(uniform_val);
            regions_sampled++;
        }
    }

    // n_copies from Poisson(lambda_c)+1
    std::poisson_distribution<int> copy_dist(lambda_c); // the param is to be specified later
    // sign
    std::bernoulli_distribution bernoulli(0.5);

    for (auto const& elem : distinct_regions)
    {
        int n_copies = copy_dist(generator) + 1;
        bool sign = bernoulli(generator);
        if (validation_test_mode)
            sign = true;

        node->c_change[elem] += (sign? n_copies : -n_copies);

        if (node->c_change.at(elem) == 0)
            node->c_change.erase(elem); //erase the zero instead of storing it

    }

    if (empty_hashmap(node->c_change))
        return nullptr;
    else
    {
        update_desc_labels(node); // to update the labels of the descendents
        // check if the subtrees are valid after updating the labels
        if (!is_valid_subtree(node) || is_redundant())
            return nullptr;

        return node;
    }

}

bool Tree::subtree_contains_negative(Node* n) {
/*
 * Returns true if any of the descendent nodes contain a value less than -ploidy in the c hashmap
 * */
    int a = -ploidy;
    vector<Node*> descendents = get_descendents(n);
    for (auto const &elem : descendents)
        for (auto const &it : elem->c)
            if (it.second < a)
                return true;
    return false;

}

vector<Node *> Tree::get_descendents(Node *n) {
    /*
     * Returns the descendents of node* n in a list in a BFS fashion.
     * Does preserve the order (e.g. parent is before the children)
     * */
    vector<Node *> descendents;

    std::stack<Node*> stk;
    stk.push(n);
    while (!stk.empty()) {
        Node* top = (Node*) stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next)
            stk.push(temp);
        descendents.push_back(top);
    }
    return descendents;
}

bool Tree::zero_ploidy_changes(Node *n) {

    vector<Node*> descendents = get_descendents(n);
    vector<int> checked_regions;

    for (auto const &node : descendents)
        for (auto const &it : node->c)
            // if the ploidy becomes 0, i.e. hashmap value = -2 for humans, then it cannot change due to biological constraints.
            if(it.second == (-1 * ploidy) && (find(checked_regions.begin(), checked_regions.end(), it.first) == checked_regions.end()))
            {
                bool does_change = region_changes(node, it.first);
                if (does_change)
                    return true;
                else
                    checked_regions.push_back(it.first);
            }
    return false;
}

bool Tree::region_changes(Node *n, u_int region_id) {
    /*
     * Returns true if the region label changes in one of the descendents
     * */

    vector<Node*> descendents = get_descendents(n);

    for (auto const &node : descendents)
        if (node->c_change[region_id] != 0)
            return true;
    return false;

}

bool Tree::is_redundant() {
    /*
     * Returns true if any of the nodes have the same c maps.
     * Worst case complexity O(n), where n: n_nodes
     * */

    unordered_map<uint64_t , unsigned> hash_map;
    for (unsigned i=0; i < all_nodes_vec.size(); i++)
    {
        if (hash_map.count(all_nodes_vec[i]->c_hash))
        {
            // found
            unsigned hash_index = hash_map.at(all_nodes_vec[i]->c_hash);
            bool eq_comparison = all_nodes_vec[i]->c == all_nodes_vec[hash_index]->c;
            return (eq_comparison);
        }
        else
        {
            // not found
            hash_map[all_nodes_vec[i]->c_hash] = i;
        }
    }

    return false;


}

Node* Tree::delete_node() {

    Node* tobe_deleted = uniform_sample(false); // TODO: replace it with the weighted scheme


    return nullptr;
}


#endif //SC_DNA_TREE_H
