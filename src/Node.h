//
// Created by Tuncel  Mustafa Anil on 7/16/18.
//

#ifndef SC_DNA_NODE_H
#define SC_DNA_NODE_H

#include <unordered_map>
#include <iostream>

struct Node{
    int id = 0;
    std::unordered_map<u_int,int> c = {};
    std::unordered_map<u_int,int> c_change= {};
    double log_score = 0.0;
    int z = 0;
    unsigned n_descendents = 1; // including itself
    Node* first_child = nullptr;
    Node* next = nullptr;
    Node* parent = nullptr;

    bool operator<(const Node &rhs) const {
        return id < rhs.id;
    }

    bool operator>(const Node &rhs) const {
        return rhs < *this;
    }

    bool operator<=(const Node &rhs) const {
        return !(rhs < *this);
    }

    bool operator>=(const Node &rhs) const {
        return !(*this < rhs);
    }

    bool operator==(const Node &rhs) const {
        return id == rhs.id;
    }

    bool operator!=(const Node &rhs) const {
        return !(rhs == *this);
    }

    // copy constructor
    Node(Node& source_node)
    {
        id = source_node.id;
        c = source_node.c;
        c_change = source_node.c_change;
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

std::ostream& operator<<(std::ostream& os, Node& n) {
    os << "node " << n.id << ": ";
    for (auto i : n.c_change)
        os << " " << i.first << ":" << i.second << ',';

    if (n.parent == nullptr)
        os << "parent: NULL";
    else
        os << "parent: " << n.parent->id;

    os << endl << "\t c values:";
    for (auto i : n.c)
        os << "\t  " << i.first << ":" << i.second << ',';

    os << endl << "\t z value: " << n.z;
    os << endl << "\t n_descendents: " << n.n_descendents;


    return os;
}

#endif //SC_DNA_NODE_H
