//
// Created by Tuncel  Mustafa Anil on 7/16/18.
//

#ifndef SC_DNA_NODE_H
#define SC_DNA_NODE_H

#include <map>
#include <iostream>
#include <stack>

using namespace std;

struct Node{
    int id = 0;
    std::map<u_int,int> c = {};
    uint64_t  c_hash = 0;
    std::map<u_int,int> c_change= {};
    double log_score = 0.0;
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

    // copy constructor
    Node(Node& source_node): c(source_node.c), c_hash(source_node.c_hash), c_change(source_node.c_change)
    {
        id = source_node.id;
        // log scores are not copied since they rely on cells
        log_score = 0.0;
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





//    os << endl << "\t c values:";
//    for (auto i : n.c)
//        os << "\t  " << i.first << ":" << i.second << ',';
//
//    os << endl << "\t z value: " << n.z;
//    os << endl << "\t n_descendents: " << n.n_descendents;


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

#endif //SC_DNA_NODE_H
