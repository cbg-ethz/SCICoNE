//
// Created by Tuncel  Mustafa Anil on 7/16/18.
//

#ifndef SC_DNA_NODE_H
#define SC_DNA_NODE_H

#include <map>
#include <iostream>
#include <stack>
#include <set>
#include <cmath>
#include <cassert>
#include "globals.cpp"

using namespace std;

struct Node{
    int id = 0;
    std::map<u_int,int> c = {};
    uint64_t  c_hash = 0;
    std::map<u_int,int> c_change= {};
    double attachment_score = 0.0;
    int z = 0;
    unsigned n_descendents = 1; // including itself
    Node* first_child = nullptr;
    Node* next = nullptr;
    Node* parent = nullptr;

    inline bool operator<(const Node &rhs) const {
        return id < rhs.id;
    }

    inline bool operator>(const Node &rhs) const {
        return rhs < *this;
    }

    inline bool operator<=(const Node &rhs) const {
        return !(rhs < *this);
    }

    inline bool operator>=(const Node &rhs) const {
        return !(*this < rhs);
    }

    inline bool operator==(const Node &rhs) const {
        return id == rhs.id;
    }

    inline bool operator!=(const Node &rhs) const {
        return !(rhs == *this);
    }

    inline int get_n_children() const;
    inline bool is_leaf() const;
    inline vector<Node*> get_descendents(bool with_n=true) const;
    inline bool children_repeat_genotype() const;
    inline double compute_event_prior(u_int n_regions) const;

    // copy constructor
    Node(Node& source_node): c(source_node.c), c_hash(source_node.c_hash), c_change(source_node.c_change)
    {
        id = source_node.id;
        // log scores are not copied since they rely on cells
        attachment_score = 0.0;
        z = source_node.z;
        n_descendents = source_node.n_descendents;
        first_child = next = parent = nullptr;
    }
    Node()
    {}
    ~Node()
    {}


    friend std::ostream& operator<<(std::ostream& os, Node& n);

};

inline std::ostream& operator<<(std::ostream& os, Node& n) {
    os << "node " << n.id << ": ";

    if (n.parent == nullptr)
        os << "p_id:NULL";
    else
        os << "p_id:" << n.parent->id;

    os << ',';
    os << '[';

    if (! n.c_change.empty())
    {
        auto last_elem_id = n.c_change.rbegin()->first;
        for (auto i : n.c_change)
        {
            os << i.first << ":" << i.second;
            if (i.first != last_elem_id)
                os << ',';
        }
    }
    os << ']';
    return os;
}

inline int Node::get_n_children() const{
/*
 * Returns the number of first order children.
 * Number of children does not have to be equal to number of descendents!
 * */

    int n_children = 0;
    for (Node* temp = this->first_child; temp != nullptr; temp=temp->next)
    {
        n_children++;
    }
    return n_children;
}

inline bool Node::is_leaf() const{
    /*
     * Returns true if the node is a leaf node, e.g. has no children
     * */
    return (this->first_child == nullptr);
}

inline vector<Node *> Node::get_descendents(bool with_n) const {
    /*
     * Returns the descendents of node* n in a list in a BFS fashion.
     * If with_n, then the descendents contain the node itself, otherwise not.
     * Does preserve the order (e.g. parent is before the children)
     *
     * */
    vector<Node *> descendents;

    std::stack<Node*> stk;
    stk.push(const_cast<Node*> (this)); // because this pointer is constant
    while (!stk.empty()) {
        Node* top = stk.top();
        stk.pop();
        for (Node* temp = top->first_child; temp != nullptr; temp=temp->next)
            stk.push(temp);
        descendents.push_back(top);
    }

    if (!with_n)
        descendents.erase(descendents.begin()); // erase the first node, which is n

    return descendents;
}

bool Node::children_repeat_genotype() const {
    /*
     * Returns true if any first order children of a node repeat the same change sign on the same region.
     * */

    std::set<std::pair<u_int, bool>> region_signs;
    // iterate over the first order children
    for (Node* temp = this->first_child; temp != nullptr; temp=temp->next)
    {

        for (auto const& x : temp->c_change)
        {
            std::pair<u_int, bool> region_sign = std::make_pair(x.first,std::signbit(x.second));
            if (!region_signs.count(region_sign)) // use count to check without initialising
            {
                // cannot be found
                region_signs.insert(region_sign);
            }
            else
            {
                // found
                return true;
            }
        }
    }

    return false;
}

double Node::compute_event_prior(u_int n_regions) const {
    /*
     * Computes and returns the event prior of the node.
     * */

    int repetition_count = 0; // the repetition count to be used in the penalisation
    double c_penalisation = c_penalise; // the penalisation coefficient, global var

    const map<u_int,int>& c_change = this->c_change;
    int v = 0;
    int v_prev = 0; // the first region is zero
    int i_prev = -1; // the initial index is -1, it'll be updated later

    auto last_elem_id = c_change.rbegin()->first;

    for (auto const &event_it : c_change)
    {
        // penalisation for repetition
        int parent_state = 0;
        try
        {
            parent_state = this->parent->c.at(event_it.first);
            int c_change_val = event_it.second;

            if (signbit(c_change_val) != signbit(parent_state))
                repetition_count++;
        }
        catch (const std::out_of_range& e)
        {
            // pass
        }
        int diff;
        if (static_cast<int>(event_it.first) - 1 != i_prev) // if the region is adjacent to its previous
        {
            int diff_right = 0 - v_prev; // the right hand side change at the end of the last consecutive region
            if (diff_right > 0)
                v += diff_right;
            v_prev = 0;
        }
        diff = event_it.second - v_prev;
        if (diff > 0)
            v += diff;
        v_prev = event_it.second;
        i_prev = event_it.first;

        if (event_it.first == last_elem_id)
        {
            int diff_last = 0 - v_prev;

            if (diff_last > 0)
            {
                v += diff_last;
                assert(v>0);
            }
        }
    }

    double pv_i = 0.0;

    /* K: max region index  */
    int K = n_regions;
    pv_i -= v*log(2*K); // the event prior
    pv_i -= c_penalisation*repetition_count; // penalise the repetitions

    return pv_i;
}

#endif //SC_DNA_NODE_H
