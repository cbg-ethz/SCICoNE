//
// Created by Tuncel  Mustafa Anil on 7/16/18.
//

#ifndef SC_DNA_NODE_H
#define SC_DNA_NODE_H

#include <unordered_map>

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

#endif //SC_DNA_NODE_H
