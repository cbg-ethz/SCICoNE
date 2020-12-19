//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <numeric>
#include <limits>
#include "Node.h"
#include "MathOp.h"
#include <random>
#include "SingletonRandomGenerator.h"
#include <set>
#include <iomanip>
#include "Utils.h"
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm> // std::remove
#include "globals.cpp"
#include "Lgamma.h"
#include "CustomExceptions.h"

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>


class Tree {
private:
    u_int n_nodes; //the number of nodes without the root
    int ploidy; // to be added to values of the unordered map for each node
    int counter = 0; // counter for the node ids to be given
    vector<int> region_neutral_states;
public:
    Node* root;
    std::vector<Node*> all_nodes_vec; // for random selection, destructor, insertion by position, iterating without order (e.g. for each node)
    u_int n_regions; // number of regions
    double posterior_score; // log posterior score of the tree
    double prior_score; // log prior score of the tree
    double od_score; // overdispersed score of the tree
    double total_attachment_score;

    // overdispersed params
    double nu;
public:
    // constructor
    Tree(u_int ploidy, u_int n_regions, vector<int> &region_neutral_states);
    // copy constructor
    Tree(Tree& source);
    // destructor
    virtual ~Tree();
    // moves
    Node *prune_reattach(bool weighted, bool validation_test_mode=false);
    void genotype_preserving_prune_reattach(double gamma);
    Node* delete_leaf();
    std::vector<Node*> swap_labels(bool weighted=false, bool validation_test_mode=false);
    Node *add_remove_events(bool weighted, bool validation_test_mode=false);
    Node *insert_delete_node(unsigned int size_limit, bool weighted, bool max_scoring);
    Node *condense_split_node(unsigned int size_limit, bool weighted, bool max_scoring);
    std::pair<std::vector<double>, std::vector<std::pair<int, int>>> gibbs_genotype_preserving_scores(double gamma);

    Node* delete_node(Node* node);
    Node* find_node(int id);
    Node* uniform_sample(bool with_root=true) const;
    Node* weighted_sample() const;
    void random_insert(std::map<u_int, int>&&);
    void insert_at(u_int pos, std::map<u_int, int>&&); // uses all_nodes_vec pos
    void insert_child(Node *pos, std::map<u_int, int>&& labels);
    void compute_tree(const vector<double> &D, const vector<int> &r);
    void compute_stack(Node *node, const vector<double> &D, double &sum_D, const vector<int> &r);
    void compute_weights();
    u_int get_n_nodes() const;
    friend std::ostream& operator<<(std::ostream& os, Tree& t);
    Tree& operator=(const Tree& other);
    //validation of tree
    bool is_redundant() const;
    // Validation of subtrees
    bool is_valid_subtree(Node* node) const;// TODO: can be a method of node instead
    bool subtree_out_of_bound(Node *n) const;// TODO: can be a method of node instead
    bool zero_ploidy_changes(Node* n) const;// TODO: can be a method of node instead
    double cost();

    vector<double> omega_condense_split(double lambda_s, bool weighted, bool max_scoring);
    vector<double> chi_condense_split(bool weighted);

    vector<double> omega_insert_delete(double lambda_r, double lambda_c, bool weighted, bool max_scoring);
    vector<double> chi_insert_delete(bool weighted);

    void load_from_file(string file);
    double get_od_root_score(const vector<int> &r, double &sum_D, const vector<double> &D) const;
    double event_prior();
private:
    void update_label(std::map<u_int,int>& c_parent, Node* node);
    void update_desc_labels(Node* node);
    bool region_changes(Node *n, u_int region_id) const;
    void copy_tree(const Tree& source_tree);
    void copy_tree_nodes(Node *destination, Node *source);
    void compute_score(Node *node, const vector<double> &D, double &sum_D, const vector<int> &r, double eta);
    void compute_root_score(const vector<int> &r);
    Node* prune(Node *pos); // does not deallocate,
    Node* insert_child(Node *pos, Node *source);
    Node* insert_child(Node *pos, Node& source);
    bool is_ancestor(Node *target, Node *curr) const;
    void destroy();
};


std::ostream& operator<<(std::ostream& os, Tree& t) {

    /*
     * Overloading the << operator
     * */

    vector<Node*> nodes = t.root->get_descendents(true);

    os << "Tree posterior: " << setprecision(print_precision) << t.posterior_score << std::endl;
    os << "Tree prior: " << setprecision(print_precision) << t.prior_score  << std::endl;
    os << "Event prior: " << setprecision(print_precision )<< t.event_prior() << std::endl;
    os << "Log likelihood: " << setprecision(print_precision) << t.total_attachment_score << std::endl;
    os << "Root score: " << setprecision(print_precision) << t.od_score << std::endl;
    os << "Tree score: " << setprecision(print_precision) << t.posterior_score + t.od_score << std::endl;
    os << "Nu: " << setprecision(print_precision) << t.nu<< std::endl;

    for (auto const &x : nodes)
    {
        os << *x << std::endl;
    }
    return os;
}


void Tree::compute_root_score(const vector<int> &r) {

    /*
     * Computes the score of the root.
     * */

    int z = 0;
    int i = 0;
    for (auto const &x : r) {
        z += x * this->region_neutral_states[i];
        i = i + 1;
    }

    root->attachment_score = 0;
    root->z = z;
}


void Tree::compute_score(Node *node, const vector<double> &D, double &sum_D, const vector<int> &r, double eta) {

    /*
     * Computes the attachment score of a node per cell.
     * Throws std::logic_error
     * */

    if (D.size() != r.size())
        throw std::logic_error("Size of the counts per cell needs to be equal to the number of regions!");

    if (node->parent == nullptr)
        compute_root_score(r);
    else
    {
        double log_likelihood = node->parent->attachment_score;
        double z = node->parent->z;
        double z_parent = node->parent->z;

        int p = ploidy;

        // for debug purposes
        double z_orig = z;

        // update z first
        for (auto const &x : node->c_change)
        {
            int cf = (node->c.count(x.first)?node->c[x.first]:0); // use count to check without initializing the not found element
            int cp_f = (node->parent->c.count(x.first) ?node->parent->c[x.first] : 0); // use count to check without initializing

            // to prevent log(0)
            p = this->region_neutral_states[x.first];
            double node_cn = (cf+p)==0?(eta):(cf+p);
            double parent_cn = (cp_f+p)==0?(eta):(cp_f+p);

            z += r[x.first] * (node_cn - parent_cn);
        }

        for (auto const &x : node->c_change)
        {
            // if log zero then use eta value, not to have -infinity
            int cf = (node->c.count(x.first)?node->c[x.first]:0);
            int cp_f = (node->parent->c.count(x.first) ?node->parent->c[x.first] : 0);

            // to prevent log(0)
            p = this->region_neutral_states[x.first];
            double node_cn = (cf+p)==0?(eta):(cf+p);
            double parent_cn = (cp_f+p)==0?(eta):(cp_f+p);

            // option for the overdispersed version
            if (is_overdispersed)
            {
                log_likelihood += lgamma(D[x.first] + nu*(node_cn)*r[x.first]);
                log_likelihood -= lgamma(D[x.first] + nu*(parent_cn)*r[x.first]);

                log_likelihood -= lgamma(nu*(node_cn)*r[x.first]);
                log_likelihood += lgamma(nu*(parent_cn)*r[x.first]);
            }
            else
                log_likelihood += D[x.first] * (log(node_cn) - log(parent_cn));
        }

        if (is_overdispersed)
        {
            // updates c independent parts
            log_likelihood += lgamma(nu*z);
            log_likelihood -= lgamma(nu*z_parent);

            log_likelihood -= lgamma(sum_D+nu*z);
            log_likelihood += lgamma(sum_D+nu*z_parent);
        }
        else
        {
            // this quantity is smt. to subtract, that's why signs are reversed
            log_likelihood -= sum_D*log(z);
            log_likelihood += sum_D*log(node->parent->z);
        }

        node->attachment_score = log_likelihood;
        node->z = z;
    }

}


void Tree::compute_tree(const vector<double> &D, const vector<int> &r) {

    /*
     * Computes the tree.
     * */

    double sum_d = std::accumulate(D.begin(), D.end(), 0.0);
    //reuse the computed sum in each node
    compute_stack(root, D, sum_d, r);

}


void Tree::compute_stack(Node *node, const vector<double> &D, double &sum_D, const vector<int> &r)
{
    /*
     * Computes the nodes scores in a top-down fashion
     * time complexity: O(n), where n is the number of nodes
     * */

    // stack based implementation
    std::stack<Node*> stk;
    stk.push(node);

    while (!stk.empty()) {
        Node* top = static_cast<Node*> (stk.top());
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        compute_score(top, D, sum_D, r, eta);
        // TODO: reuse this part of the code (the recursive iteration of nodes)
    }

//    alternative recursive implementation
//    compute_score(node, D, sum_D, r);
//    for (Node* temp = node->first_child; temp != nullptr; temp=temp->next) {
//        compute_stack(temp, D, sum_D, r);
//    }
}


Tree::Tree(u_int ploidy, u_int n_regions, vector<int> &region_neutral_states)
{
    // Tree constructor
    root = new Node();
    // compute and store the hash
    vector<int> keys_values = {};

    uint64_t size_for_hash = keys_values.size() * sizeof(keys_values[0]);
    uint64_t c_hash = Utils::calculate_hash(&keys_values[0], size_for_hash);
    root->c_hash = c_hash;

    this->region_neutral_states = region_neutral_states;

    this->ploidy = ploidy;
    this->n_regions = n_regions;
    n_nodes = 0;
    prior_score = 0.0;
    total_attachment_score = 0.0;
    posterior_score = 0.0;
    od_score = 0.0;
    // creates a copy of the root ptr and stores it in the vector

    nu = 1.0 / static_cast<double>(n_regions);
    all_nodes_vec.push_back(root);

}


Tree::~Tree() {
    destroy();
}


Node* Tree::uniform_sample(bool with_root) const{

    /*
     * Returns a uniformly sampled node pointer.
     *
     * */

    int rand_val = 0;
    if (all_nodes_vec.empty())
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes_vec.size() ==1)
    {
        if (not with_root)
            throw std::logic_error("root is the only node and cannot be sampled");
    }
    else
    {
        int min_val = with_root ? 0 : 1;
        rand_val = MathOp::random_uniform(min_val, all_nodes_vec.size() - 1);
    }

    return all_nodes_vec[rand_val];
}


Node * Tree::insert_child(Node *pos, Node *source) {
/*
 * Inserts the source node to the pos node.
 * Updates the labels by calling update_label with the parent node.
 * */

    // set the parent from the child
    source->parent = pos;

    // set the unique id
    if (source->id == 0)
        source->id = ++counter;
    // if the id is set then keep it

    all_nodes_vec.push_back(source);
    n_nodes++;

    // insert
    if (pos->is_leaf())
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
    /*
     * Destroys the tree
     * */

    for (auto elem: all_nodes_vec)
    {
        //std::cout<<"deleting " << elem->attachment_score <<std::endl;
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
    /*
     * Inserts a node to a uniform random position.
     * */

    Node* pos = uniform_sample(true);
    insert_child(pos, std::move(labels));

}


void Tree::insert_at(u_int pos, std::map<u_int, int> && labels) {
    /*
     * Inserts the node by its all_nodes_vector position
     * */
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
    unsigned long size_for_hash = keys_values.size() * sizeof(keys_values[0]);
    uint64_t c_hash = Utils::calculate_hash(&keys_values[0], size_for_hash);

    node->c_hash = c_hash;
}


Tree::Tree(Tree &source) {
    /* Copy constructor */

    if (source.root == nullptr)
        this->root = nullptr;
    else
        this->copy_tree(source);
}


void Tree::copy_tree(const Tree& source_tree) {
    /*
     * Copies the source tree
     */

    this->region_neutral_states = source_tree.region_neutral_states;
    this->ploidy = source_tree.ploidy;
    this->total_attachment_score = source_tree.total_attachment_score;
    this->counter = source_tree.counter;
    this->prior_score = source_tree.prior_score;
    this->posterior_score = source_tree.posterior_score;
    this->od_score = source_tree.od_score;
    this->n_regions = source_tree.n_regions;
    this->nu = source_tree.nu;

    // copy the nodes using struct copy constructor
    this->root = new Node(*source_tree.root);
    this->all_nodes_vec.push_back(root); // all nodes cannot be copied since 2 trees cannot have nodes pointing to the same address
    copy_tree_nodes(this->root, source_tree.root); // the nodes of the source tree are inserted into the destination tree
    this->n_nodes = source_tree.n_nodes; // copy this after the insertions are done (insertions change this value).

 }

void Tree::copy_tree_nodes(Node *destination, Node *source) {
    /*
     * Copies the tree from source to destination in a BFS fashion
     * Uses a queue of Node* pairs to keep the copying order.
     * */

    deque<pair<Node*,Node*>> to_copy;

    for (Node* temp = source->first_child; temp != nullptr; temp=temp->next) // initialize with the first children
        to_copy.push_back({destination,temp});

    while(!to_copy.empty())
    {
        pair<Node*,Node*> top = to_copy.back();
        to_copy.pop_back();
        Node* new_destination = insert_child(top.first, *top.second);
        for (Node* temp = top.second->first_child; temp != nullptr; temp=temp->next)
            to_copy.push_back({new_destination,temp});
    }
}

Node* Tree::prune_reattach(bool weighted, bool validation_test_mode) {
    /*
     * Prunes a node and reattaches it to another node which is not among the descendents of the pruned node.
     * Returns the pruned node (which also happens to be the attached node)
     * Requires more than two nodes to perform.
     * Throws InvalidMove exception.
     *
     * */

    if (all_nodes_vec.size() <= 2)
        throw InvalidMove("prune and reattach move does not make sense when there is only one node besides the root");

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
        Node* top = static_cast<Node*> (stk.top());
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

        std::map<u_int,int> pruned_c;

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
     * Throws std::logic_error
     * */

    if (pos == root)
        throw std::logic_error("root cannot be pruned!");

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


void Tree::load_from_file(string file) {
    /*
     * Loads the tree from file
     * Throws std::logic_error, std::runtime_error, InvalidTree exception
     * */

    // first destroy the tree if it is not empty (or 1 node only)
    if (this->root->first_child != nullptr)
    {
        throw std::logic_error("Tree to read into must be empty!");
    }


    //string file = "10nodes_0regions_100reads_size_limit_test_tree_inferred_segmented.txt";
    std::ifstream infile(file);
    std::string line;

    std::getline(infile, line); // tree posterior
    std::getline(infile, line); // tree prior
    std::getline(infile, line); // event prior
    std::getline(infile, line); // log likelihood
    std::getline(infile, line); // root score
    std::getline(infile, line); // tree score
    std::getline(infile, line); // nu value
    std::getline(infile, line); // pass the first 5 lines, including the root
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        string a, b, c ;
        if (!(iss >> a >> b >> c)) { std::cout << "Error reading tree from file" << std::endl; break; } // error
        b.pop_back();
        int node_id = stoi(b);
        if (node_id > this->counter)
            this->counter = node_id + 1;
        string del3 = ",[";
        string delim_cp = "]";
        string token = c.substr(0, c.find(del3));
        string token_r = c.substr(c.find(del3)+1, c.length());
        // TODO: if token_r does not start and end with [ and ] then throw

        if (token_r.front() != '[' || token_r.back() != ']')
            throw std::runtime_error(" Incorrect tree format to parse! \nThe rows should be space separated.");

        int parent_id = stoi(token.substr(5,token.size()));
        string s = token_r.substr(1, token_r.find(delim_cp)-1);
        // cout <<" node id: " << node_id << " parent id: " << parent_id <<"s=" <<s<<  endl;
        // process pair (a,b)
        size_t pos = 0;
        vector<string> pairs;
        string delimiter= ",";
        while ((s.find(delimiter)) != std::string::npos)
        {
            pos = s.find(delimiter);
            pairs.push_back(s.substr(0, pos));
            s.erase(0, pos+ delimiter.length());
        }
        pairs.push_back(s);

        // create node
        Node* child = new Node();
        child->id = node_id;

        // fill the c_change map
        for (auto const& s : pairs)
        {
            u_int key = stoi(s.substr(0,s.find(':')));
            int value = stoi(s.substr(s.find(':')+1, s.length()));
            child->c_change[key] = value;
        }


        // find the node with a given id and insert there
        //t.in
        Node* parent = nullptr;
        for (Node* ptr : all_nodes_vec)
        {
            if (ptr->id == parent_id)
            {
                parent = ptr;
                break;
            }
        }
        insert_child(parent, child);
    }

    if (!is_valid_subtree(this->root) || is_redundant())
        throw InvalidTree("The loaded tree is invalid!");

    compute_weights();

}

void Tree::compute_weights() {

    /*
     * Computes the weights of all of the nodes.
     * a bottom-up approach
     * time complexity: O(n), where n is the number of nodes
     * */


    // stack to allow the bottom-up computations
    std::stack<Node*> cmp_stack;


    //stack for the tree traversal
    std::stack<Node*> stk;
    stk.push(root); //start with the root

    while (!stk.empty()) {
        Node* top = static_cast<Node*> (stk.top());
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        cmp_stack.push(top);
    }

    while (!cmp_stack.empty())
    {
        Node* top = static_cast<Node*> (cmp_stack.top());
        cmp_stack.pop();
        unsigned n_descendents = 1;

        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next)
            n_descendents += temp->n_descendents;

        top->n_descendents = n_descendents;
    }
}


Node *Tree::weighted_sample() const{

    /* Samples and returns the node proportional to the weight
     * Throws std::length_error
     * */

    if (all_nodes_vec.empty())
        throw std::length_error("length of nodes must be bigger than zero, in order to sample from the tree");
    else if (all_nodes_vec.size() ==1)
        throw std::length_error("there is only 1 element which is the root and root cannot be sampled");
    else
    {
        // get the subvector
        auto first = all_nodes_vec.begin() + 1;
        auto last = all_nodes_vec.end();
        vector<Node*> nodes_to_sample(first, last);

        vector<float> weights;
        for (auto const &x : nodes_to_sample)
        {
            float weight = (1.0f / x->n_descendents); // weights are inversely proportional to n_descendents
            weights.push_back(weight);
        }

        std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;
        boost::random::discrete_distribution<> d(weights.begin(), weights.end());

        unsigned sample = d(generator);

        return nodes_to_sample[sample];

    }

}


bool Tree::is_ancestor(Node *target, Node *curr) const{
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
     * Throws InvalidMove, InvalidTree and std::logic_error.
     * */

    if (all_nodes_vec.size() <= 2)
        throw InvalidMove("swapping labels does not make sense when they is only one node besides the root");

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
            node1 = all_nodes_vec[3]; //without the root
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

    // if the labels to be swapped are the same, then reject the move
    if(node1->c_change == node2->c_change)
        throw std::logic_error("swapping 2 nodes with the same labels does not make sense, the move will be rejected");

    // if two sibling leaves are swapped, then reject
    if(node1->parent == node2->parent && node1->first_child == nullptr && node2->first_child == nullptr)
        throw std::logic_error("swapping the sibling leaves will have no impact, move will be rejected.");

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
            throw InvalidTree("Swap labels move created an invalid tree.");

    }
    if (is_redundant())
        throw InvalidTree("Swap labels move created a redundant tree");


    return return_nodes;
}


void Tree::update_desc_labels(Node *node) {

    /*
     * Updates the c vectors of the node and all of its descendents
     * */


    std::stack<Node*> stk;
    stk.push(node);

    while (!stk.empty()) {
        Node* top = static_cast<Node*>(stk.top());
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next) {
            stk.push(temp);
        }
        if (top == this->root) // delete move can cause this to happen (if a F.O. child of root is deleted) and top->parent will be NULL
            continue;
        update_label(top->parent->c, top);
    }
}


bool Tree::is_valid_subtree(Node *node) const{
    /*
     * Returns true if the subtree is valid
     * */

    if (subtree_out_of_bound(node)) // does it for the subtree
        return false;

    if (zero_ploidy_changes(node)) // does it for the subtree
        return false;

    return true;

}


Node* Tree::add_remove_events(bool weighted, bool validation_test_mode) {

    /*
     * Adds and removes events to a node.
     * Throws logical error if the "n_regions_to_sample" value, (drawn from a poisson) is bigger than total number of regions available.
     * Returns the pointer to the node being affected.
     * */

    if (all_nodes_vec.size() <= 1)
        throw InvalidMove("Adding or removing events does not make sense when there is 1 node or less. Root has to be neutral.");

    Node* node;

    if (validation_test_mode)
    {
        node = all_nodes_vec[5];
        node->c_change = {{3,-2}};
    }
    else
    {
        if (weighted)
            node = weighted_sample();
        else
            node = uniform_sample(false); //without the root

        std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;

        // n_regions from Poisson(lambda_R)+1
        boost::random::poisson_distribution<int> poisson_dist(lambda_r); // the param is to be specified later
        int n_regions_to_sample = poisson_dist(generator) + 1;
        // sample n_regions_to_sample distinct regions uniformly
        int n_regions = this->n_regions;
        int regions_sampled = 0;
        std::set<u_int> distinct_regions;

        // otherwise we cannot sample distinct uniform regions
        if (n_regions_to_sample > n_regions)
            throw std::logic_error("the number of distinct regions to sample cannot be bigger than the total number of distinct regions available");

        while (regions_sampled < n_regions_to_sample)
        {
            u_int uniform_val = static_cast<u_int> (MathOp::random_uniform(0, n_regions-1));
            if (distinct_regions.find(uniform_val) == distinct_regions.end())
            {
                distinct_regions.insert(uniform_val);
                regions_sampled++;
            }
        }

        // n_copies from Poisson(lambda_c)+1
        boost::random::poisson_distribution<int> copy_dist(lambda_c); // the param is to be specified later
        // sign
        boost::random::bernoulli_distribution<double> bernoulli(0.5);
        for (auto const& elem : distinct_regions)
        {
            int n_copies = copy_dist(generator) + 1;
            bool sign = bernoulli(generator);

            node->c_change[elem] += (sign? n_copies : -n_copies);

            if (node->c_change.at(elem) == 0)
                node->c_change.erase(elem); //erase the zero instead of storing it
        }
    }

    if (Utils::is_empty_map(node->c_change))
        return nullptr; //TODO: maybe throw a certain exception here
    else
    {
        update_desc_labels(node); // to update the labels of the descendents
        // check if the subtrees are valid after updating the labels
        if (!is_valid_subtree(node) || is_redundant())
            return nullptr; //TODO: maybe throw a certain exception here

        return node;
    }
}


bool Tree::subtree_out_of_bound(Node *n) const{
/*
 * Returns true if any of the descendent nodes contain a value less than -ploidy in the c hashmap
 * */

    int lb = -ploidy;
    int ub = copy_number_limit; // global variable
    vector<Node*> descendents = n->get_descendents(true);
    for (auto const &elem : descendents)
        for (auto const &it : elem->c) {
            lb = -this->region_neutral_states[it.first];
            if (it.second < lb || it.second > ub)
                return true;
        }
    return false;

}


bool Tree::zero_ploidy_changes(Node *n) const{
/*
 * Returns true if the ploidy value changes after hitting zero.
 * This is a domain knowledge constraint.
 * */

    vector<Node*> descendents = n->get_descendents(true);

    for (auto const &node : descendents)
    {
        if (node->id == 0)
            continue; // root cannot have events

        for (auto const &it : node->parent->c)
            if(it.second <= (-1 * this->region_neutral_states[it.first])) {
                bool does_change = region_changes(node, it.first);
                if (does_change)
                    return true;
            }
    }

    return false;

}


bool Tree::region_changes(Node *n, u_int region_id) const{
    /*
     * Returns true if the region label changes in one of the descendents
     * */

    vector<Node*> descendents = n->get_descendents(true);

    for (auto const &node : descendents)
        if (node->c_change.count(region_id) != 0)
            return true;
    return false;

}


bool Tree::is_redundant() const {
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


Node *Tree::insert_delete_node(unsigned int size_limit, bool weighted, bool max_scoring) {
    /*
     * Adds or deletes nodes move, that takes the mcmc transition probabilities into account.
     * Returns the node to perform partial score computation on.
     * Throws std::logic_error
     * */

    Node* return_node = nullptr;

    vector<double> chi = chi_insert_delete(weighted); // add weights
    vector<double> omega = omega_insert_delete(lambda_r, lambda_c, weighted, max_scoring); // delete weights

    // sample the number of regions to be affected with Poisson(lambda_r)+1
    std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;

    // 0.5 prob bernoulli
    boost::random::bernoulli_distribution<double> bernoulli_05(0.5);

    int K = this->n_regions;

    boost::random::uniform_real_distribution<double> prob_dist(0.0,1.0);
    double rand_val = prob_dist(generator); // to be btw. 0 and 1

    if (rand_val < 0.5)
    {
        if (verbosity > 2)
            std::cout<<"insert is chosen"<<std::endl;
        // add is chosen
        if (all_nodes_vec.size() >= size_limit)
            throw std::logic_error("Tree size limit is reached, insert node move will be rejected!");

        boost::random::discrete_distribution<>* dd;
        dd = new boost::random::discrete_distribution<>(chi.begin(),chi.end());

        u_int pos_to_insert = (*dd)(generator); // this is the index of the all_nodes_vector.
        delete dd;

        // create the map
        // create a map, fill it properly with r amount of labels
        map<u_int, int> distinct_regions;

        Utils::random_initialize_labels_map(distinct_regions, K); // modifies the distinct_regions

        Node* parent = all_nodes_vec[pos_to_insert];
        // retrieve the first order children of pos_to_insert node
        vector<Node*> first_order_children;
        for (Node* temp = parent->first_child; temp != nullptr; temp=temp->next)
        {
            first_order_children.push_back(temp);
        }

        this->insert_at(pos_to_insert, static_cast<map<u_int, int> &&>(distinct_regions)); // cast the expression to r_value
        Node* new_node = all_nodes_vec.back(); // the last inserted elem, e.g. new node

        // foreach first_order_child, with 0.5 prob. change their parent to become the last element.
        // Change the parent with 0.5 prob
        for (Node* elem : first_order_children)
        {
            bool rand = bernoulli_05(generator); // either true or false with 0.5 prob
            if (rand)
            {
                Node* pruned_child = prune(elem); // prune from the old parent
                insert_child(new_node, pruned_child); // insert into next parent
            }
        }
        return_node = new_node;
    }

    else // delete is chosen
    {
        if (verbosity > 2)
            std::cout<<"delete is chosen"<<std::endl;
        if (all_nodes_vec.size() <= 1)
            throw std::logic_error("Root cannot be deleted, delete move will be rejected");

        boost::random::discrete_distribution<>* dd;
        dd = new boost::random::discrete_distribution<>(omega.begin(),omega.end());

        u_int64_t idx_tobe_deleted = (*dd)(generator);
        vector<Node*> descendents_of_root = this->root->get_descendents(false); // without the root
        Node* tobe_deleted = descendents_of_root[idx_tobe_deleted];

        if (verbosity > 2)
            std::cout<<"node to be deleted:" << tobe_deleted->id <<std::endl;
        delete dd;
        return_node = delete_node(tobe_deleted); // returns the parent of the deleted node
    }

    //update the c vectors of the parent and its new descendents
    update_desc_labels(return_node);
    // check if the subtrees are valid after updating the labels
    if (!is_valid_subtree(return_node) || is_redundant())
        return nullptr;
    // recompute the weights after the tree structure is changed
    this->compute_weights();
    return return_node;
}


Node *Tree::condense_split_node(unsigned int size_limit, bool weighted, bool max_scoring) {
    /*
     * Condenses two nodes into one or splits a node into two.
     * Throws std::logic_error, InvalidMove
     * */

    if (all_nodes_vec.size() <= 1)
        throw InvalidMove("condense or split does not make sense when there is 1 node or less. ");


    Node* return_node = nullptr;
    vector<double> chi = this->chi_condense_split(weighted); // split weights
    vector<double> omega = this->omega_condense_split(lambda_s, weighted, max_scoring); // condense weights

    std::mt19937 &generator = SingletonRandomGenerator::get_instance().generator;
    // n_regions from Poisson(lambda_S)+1
    boost::random::poisson_distribution<int> poisson_s(lambda_s); // the param is to be specified later
    boost::random::bernoulli_distribution<double> bernoulli_05(0.5);

    boost::random::uniform_real_distribution<double> prob_dist(0.0,1.0);
    double rand_val = prob_dist(generator); // to be btw. 0 and 1

    vector<Node*> descendents_of_root = this->root->get_descendents(false); // without the root
    if (rand_val < 0.5)
    {
        // split is chosen

        if (all_nodes_vec.size() >= size_limit)
            throw std::logic_error("Tree size limit is reached, split move will be rejected!");


        boost::random::discrete_distribution<>* dd;

        dd = new boost::random::discrete_distribution<>(chi.begin(),chi.end());

        u_int pos_to_insert = static_cast<u_int>((*dd)(generator)); // this is the index of the descendents_of_root.
        delete dd;

        Node* parent = descendents_of_root[pos_to_insert]; // this node will be split

        if (parent->c_change.size() == 1)
            throw logic_error("Nodes with single events in a single region cannot be split");

        // compute v_p and v_c
        map<u_int,int> v_p;
        map<u_int,int> v_c;

        for (auto const& x : parent->c_change) // for a node // TODO: make a method that takes 2 maps by reference and updates them, reuse it if split move is chosen
        {
            bool sign = bernoulli_05(generator);
            int beta = 2 * sign -1;
            int lambda = poisson_s(generator);
            int v_k = x.second;
            if (v_k % 2 == 0) // even value
                v_p[x.first] = v_k/2 + beta*(0.5 + lambda); // results in integer
            else // odd value
                v_p[x.first] = v_k/2 + beta*lambda;
            v_c[x.first] = v_k - v_p[x.first]; // v_c = v - v_p

            // make sure no zero values are created
            if (v_p.at(x.first) == 0)
                v_p.erase(x.first);
            if (v_c.at(x.first) == 0)
                v_c.erase(x.first);
        }
        parent->c_change = v_p; // the parent c_change becomes v_p

        // retrieve the first order children of parent, because some will be passed to the new children
        vector<Node*> first_order_children;
        for (Node* temp = parent->first_child; temp != nullptr; temp=temp->next)
        {
            first_order_children.push_back(temp);
        }

        this->insert_child(parent, static_cast<map<u_int,int>&&>(v_c)); // insert the child with v_c
        Node* new_node = all_nodes_vec.back(); // the last inserted elem, e.g. new node
        // foreach first_order_child, with 0.5 prob. change their parent to become the last element.
        // Change the parent with 0.5 prob
        for (Node* elem : first_order_children)
        {
            bool rand = bernoulli_05(generator); // either true or false with 0.5 prob
            if (rand)
            {
                Node* pruned_child = prune(elem); // prune from the old parent
                insert_child(new_node, pruned_child); // insert into next parent
            }
        }
        return_node = parent;
    }
    else
    {
        // condense is chosen
        // condense is delete with its c_change values passed to the parent

        if (all_nodes_vec.size() <= 2)
            throw std::logic_error("condensing nodes does not make sense when there are 2 or less nodes. Root has to be neutral.");


        boost::random::discrete_distribution<>* dd;
        dd = new boost::random::discrete_distribution<>(omega.begin(),omega.end());
        u_int64_t idx_tobe_deleted = (*dd)(generator); // this is the index of the descendents_of_root,
        delete dd;

        Node* tobe_deleted = descendents_of_root[idx_tobe_deleted];

        // check if suitable for condensing
        if (Utils::is_empty_map(tobe_deleted->c_change) || Utils::is_empty_map(tobe_deleted->parent->c_change)) // if one of them is an empty map then reject the move
            throw std::logic_error("cannot condense a node with an empty node (root)");

        map<u_int,int> condensed_c_change = tobe_deleted->parent->c_change;

        for (auto const& x : tobe_deleted->c_change)
        {
            int new_val = x.second;
            if (condensed_c_change.count(x.first))
                new_val += condensed_c_change[x.first];
            if (new_val == 0)
                condensed_c_change.erase(x.first);
//                throw std::logic_error("Any of the events cannot cancel completely upon condense");
            else
                condensed_c_change[x.first] = new_val;
        }

        if (condensed_c_change.size() == 1)
        {
            for (auto const &i : condensed_c_change)
            {
                if(i.second == 1)
                    throw std::logic_error("Nodes cannot be combined if they result in a single event");
            }
        }

        tobe_deleted->parent->c_change = condensed_c_change;
        return_node = delete_node(tobe_deleted); // returns the parent of the deleted node

    }
    //update the c vectors of the parent and its new descendents
    update_desc_labels(return_node);
    // check if the subtrees are valid after updating the labels
    if (!is_valid_subtree(return_node) || is_redundant())
        return nullptr;
    // recompute the weights after the tree structure is changed
    this->compute_weights();
    return return_node;
}

Node* Tree::delete_node(Node *node) {
    /*
     * Deletes the node.
     * Assigns the first order children of the deleted node to the parent of the deleted node.
     * Returns the pointer to the parent of the deleted node.
     * */

    Node* tobe_deleted = node;

    tobe_deleted = prune(tobe_deleted); // the node is pruned
    Node* parent_of_deleted = tobe_deleted->parent;

    // retrieve the first order children and prune & reattach them to parent_of_deleted
    vector<Node*> first_order_children;
    for (Node* temp = tobe_deleted->first_child; temp != nullptr; temp=temp->next)
    {
        first_order_children.push_back(temp);
    }
    for (Node* elem : first_order_children)
    {
        Node* pruned_child = prune(elem); // prune removes the next link, you cannot prune while iterating over children
        pruned_child = insert_child(parent_of_deleted, pruned_child);
    }

    delete tobe_deleted; // deallocate it.

    return parent_of_deleted; // return the node that is going to be used for the partial trees scoring
}


vector<double> Tree::omega_condense_split(double lambda_s, bool weighted, bool max_scoring) {
    /*
     * Returns the omega probabilities computed for a tree for the condense/split move.
     * Omega vector is representing the condense weights
     * */

    vector<double> omega; // condense weights
    vector<Node*> descendents_of_root = this->root->get_descendents(false); // without the root

    for (auto const &node : descendents_of_root)
    {
        double omega_val = 1.0;
        if (not max_scoring)
          omega_val = MathOp::compute_omega_condense_split(node, lambda_s, this->n_regions);
        if (weighted)
            omega_val /= (node->parent->n_descendents-1);
        omega.push_back(omega_val);
    }

    return omega;
}


vector<double> Tree::chi_condense_split(bool weighted) {
    /*
     * Returns the chi probabilities computed for a tree for the condense/split move.
     * Chi vector is representing the split weights
     * */

    vector<double> chi; // split weights

    vector<Node*> descendents_of_root = this->root->get_descendents(false); // without the root
    // compute chi and xi
    for (auto const &node : descendents_of_root) // all nodes without the root
    {
        if (node->c_change.size() <= 1)
            continue;
        else
        {
            double chi_val = pow(2, node->get_n_children());
            if (weighted)
                chi_val /= (node->n_descendents + 1);

            chi.push_back(chi_val);
        }

    }

    return chi;
}


vector<double> Tree::chi_insert_delete(bool weighted) {
    /*
     * Returns the chi probabilities computed for the insert/delete move.
     * Chi vector is representing the insert weights
     * */

    vector<double> chi; // add weights
    // iterate over all_nodes, use node->n_descendents-1 as n_children
    // using all_nodes_vec is safe here since neither its order nor its size will change until the action and the node are chosen
    for (auto const &node : all_nodes_vec)
    {
        double chi_val;
        if (weighted)
            chi_val = pow(2, node->get_n_children()+1) / (node->n_descendents+1);
        else
            chi_val = pow(2, node->get_n_children()); // chi is to be computed for the n_first order children

        chi.push_back(chi_val);
    }

    return chi;

}


vector<double> Tree::omega_insert_delete(double lambda_r, double lambda_c, bool weighted, bool max_scoring) {
    /*
     * Returns the omega probabilities computed for the insert/delete move.
     * Omega vector is representing the delete weights
     * */

    vector<double> omega; // delete weights
    u_int K = this->n_regions;

    vector<Node*> all_nodes = root->get_descendents(false); // without root
    for (auto const &node : all_nodes) { // computes the omega vector
        double omega_val = 1.0;
        if (not max_scoring)
          omega_val = MathOp::compute_omega_insert_delete(node, lambda_r, lambda_c, K);
        if (weighted)
            omega_val = omega_val/node->n_descendents;
        omega.push_back(omega_val);
    }

    return omega;
}


double Tree::cost() {
    /*
     * Returns the cost of recomputing the children nodes in case of a move
     * */

    vector<Node*> nodes = root->get_descendents(false); // get all of the children of root

    // vector<double> weights;
    double zeta = 0.0;
    for (auto const &x : nodes)
    {
        float weight = (1.0 / x->n_descendents); // weights are inversely proportional to n_descendents
        // weights.push_back(weight);
        zeta += weight;
    }

    return zeta;
}


double Tree::get_od_root_score(const vector<int> &r, double &sum_D, const vector<double> &D) const{
    /*
     * Returns the overdispersed root score
     * */

    double od_root_score = 0.0;

    if (is_overdispersed)
    {
        int z = 0;
        int k = 0;
        for (auto const &x : r) {
            z += x * this->region_neutral_states[k];
            k = k + 1;
        }
        od_root_score += lgamma(nu*z);
        od_root_score -= lgamma(sum_D+nu*z);

        for (u_int i = 0; i < r.size(); ++i)
        {
            od_root_score += lgamma(D[i] + nu*this->region_neutral_states[i]*r[i]);
            od_root_score -= lgamma(nu*this->region_neutral_states[i]*r[i]);
        }
    }
    else
    {
        int sum_r = std::accumulate(r.begin(),r.end(),0);
        double term1 = -sum_D*std::log(sum_r);
        double term2 = 0.0;
        for (u_int i = 0; i < r.size(); ++i)
            term2 += D[i]*std::log(r[i]);

        od_root_score += term1+term2;
    }
    return od_root_score;

}


double Tree::event_prior() {
    /*
     * Computes and returns the tree event prior
     * Throws std::logic_error
     * */

    int n = this->get_n_nodes(); //n_nodes

    vector<double> p_v;
    for (auto it = this->all_nodes_vec.begin()+1; it != this->all_nodes_vec.end(); ++it) // without the root
    {
        Node* node = *it;
        double pv_i = node->compute_event_prior(this->n_regions);
        p_v.push_back(pv_i);

    }

    if (n != static_cast<int>(p_v.size()))
        throw std::logic_error("number of nodes needs to be equal to the size of the event prior vector");

    double PV = 0.0;
    PV += std::accumulate(p_v.begin(), p_v.end(), 0.0);
    PV -= Lgamma::get_val(n+1);


    return PV;
}


void Tree::genotype_preserving_prune_reattach(double gamma) {
    /*
     * Prunes a node and reattaches it to another node which is not among the descendents of the pruned node.
     * Preserves the genotypes of all nodes.
     * Performs gibbs sampling and returns the tree.
     * Requires more than two nodes to perform.
     * Throws InvalidMove, std::logic_error
     * */

    if (this->all_nodes_vec.size() <= 2)
        throw InvalidMove("prune and reattach move does not make sense when there is only one node besides the root");


    std::vector<double> all_possible_scores; // event priors of all valid attachments
    std::vector<std::pair<int,int>> prune_attach_indices;

    std::tie(all_possible_scores, prune_attach_indices) = gibbs_genotype_preserving_scores(gamma);


    double max_tree_score = *std::max_element(all_possible_scores.begin(), all_possible_scores.end());
    for (u_int j = 0; j < all_possible_scores.size(); ++j)
        all_possible_scores[j] = std::exp(all_possible_scores[j] - max_tree_score);

    // sample from the tree scores
    std::mt19937 &gen = SingletonRandomGenerator::get_instance().generator;
    boost::random::discrete_distribution<> d(all_possible_scores.begin(), all_possible_scores.end());
    unsigned sampled_tree_index = d(gen);

    // perform prune & reattach if needed
    if (sampled_tree_index != 0)
    {
        std::pair<int,int> prune_attach_pos = prune_attach_indices[sampled_tree_index];
        Node* to_prune = this->find_node(prune_attach_pos.first);
        if (to_prune == nullptr)
            throw std::logic_error("Node to prune could not be found in the tree. Move will be rejected.");

        Node* pruned_node = this->prune(to_prune);

        Node* attach_pos = this->find_node(prune_attach_pos.second);
        if (attach_pos == nullptr)
            throw std::logic_error("Node to attach could not be found in the tree. Move will be rejected.");

        // update c_change
        pruned_node->c_change = Utils::map_diff(pruned_node->c, attach_pos->c);

        if (pruned_node->c_change.empty()) //reject the attachment if c_change becomes empty
            throw std::logic_error("Empty c_label node created. Move will be rejected.");

        this->insert_child(attach_pos, pruned_node);

        // recompute the weights after the tree structure is changed
        this->compute_weights();

        to_prune = pruned_node = attach_pos = nullptr;

    }
    else
    {
        throw std::logic_error("The node is reattached to the same position.");
    }
}


Node *Tree::find_node(int id) {
    /*
     * Finds and returns the pointer to the node with the given id.
     * Returns nullptr if the node is not found.
     * id: The node id to look for
     * */

    vector<Node*> nodes = this->root->get_descendents(true);

    for (Node* const x : nodes)
        if (x->id == id)
            return x;

    return nullptr; // not found
}


std::pair<std::vector<double>, std::vector<std::pair<int, int>>> Tree::gibbs_genotype_preserving_scores(double gamma)
{
    /*
     * Computes the event priors for all possible attachments for every node.
     * Returns the pair of all scores and their corresponding prune and attach indices
     * */

    std::vector<double> all_possible_scores; // event priors of all valid attachments
    std::vector<std::pair<int,int>> prune_attach_indices;

    // the score for the original tree
    all_possible_scores.push_back(0.0); //exp(0) is 1.
    prune_attach_indices.emplace_back(std::make_pair(std::nan(""),std::nan("")));

    for (u_int i = 1; i < this->all_nodes_vec.size(); ++i) // i=1, excluding the root
    {
        // i: prune position index
        int prune_pos_idx = this->all_nodes_vec[i]->id;
        Node* prune_pos = this->all_nodes_vec[i];
        double original_node_event_prior = prune_pos->compute_event_prior(this->n_regions);

        // copy all nodes
        std::vector<Node*> destination_nodes = this->all_nodes_vec;

        // remove all the descendents of the prune_pos
        std::stack<Node*> stk;
        stk.push(prune_pos);

        while (!stk.empty()) {
            Node* top = static_cast<Node*> (stk.top());
            stk.pop();
            for (Node *temp = top->first_child; temp != nullptr; temp = temp->next) {
                stk.push(temp);
            }
            destination_nodes.erase(std::remove(destination_nodes.begin(), destination_nodes.end(), top), destination_nodes.end());
        }

        for (Node* attach_pos : destination_nodes)
        {

            int attach_pos_idx = attach_pos->id;
            //do not recompute if you attach at the same pos
            if (prune_pos->parent->id != attach_pos_idx)
            {

                // Think as if a node is pruned and reattached somewhere without having to create the entire tree
                // copy the nodes
                Node copy_pruned(*prune_pos);
                Node copy_attach_pos(*attach_pos);

                // create parent-child relationship
                copy_attach_pos.first_child = &copy_pruned;
                copy_pruned.parent = &copy_attach_pos;

                // update c_change
                copy_pruned.c_change = Utils::map_diff(copy_pruned.c, copy_attach_pos.c);

                if (copy_pruned.c_change.empty()) //reject the attachment if c_change becomes empty
                    continue;


                // compute event prior
                double event_prior = copy_pruned.compute_event_prior(this->n_regions);

                // check if invalid tree occurs

                // 1. insert the modified node for a while
                this->insert_child(attach_pos, &copy_pruned);
                // 2. check for zero ploidy changes
                bool zero_ploidy_changes = this->zero_ploidy_changes(attach_pos);
                // 3. prune the copy_pruned back
                this->prune(&copy_pruned);
                // 4. decide
                if (zero_ploidy_changes)
                    continue;
                else
                {
                    prune_attach_indices.emplace_back(std::make_pair(prune_pos_idx,attach_pos_idx));
                    double score_diff = (event_prior - original_node_event_prior) * gamma;
                    all_possible_scores.push_back(score_diff);
                }
            }
        }
    }

    return std::make_pair(all_possible_scores, prune_attach_indices);
}


Node *Tree::delete_leaf() {
    /*
     * Uniformly samples a leaf node and deletes it.
     * Throws InvalidMove exception if there is only one node besides the root.
     * This move is intended to delete leaves that no cells attach.
     * Returns the parent of the deleted node
     * */

    if (all_nodes_vec.size() <= 2)
        throw InvalidMove("delete leaf move does not make sense when there is only one node besides the root");

    vector<Node*> all_leaves;
    // get all the leaves
    for (Node* node : all_nodes_vec) {
        if (node->is_leaf())
            all_leaves.push_back(node);
    }

    // sample uniformly
    int rand_val = MathOp::random_uniform(0, all_leaves.size() - 1);

    Node* to_delete = all_leaves[rand_val];

    // perform delete
    Node* parent_of_deleted = this->delete_node(to_delete);

    return parent_of_deleted;
}


#endif //SC_DNA_TREE_H
