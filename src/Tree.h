//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>

struct Node{
    std::unordered_map<std::string,int> c;
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

    bool is_leaf(Node*) const;
    Node* uniform_select();
    void random_insert(std::unordered_map<std::string, int>&&);
    void insert_at(u_int pos, std::unordered_map<std::string, int>&&);
    void insert_child(Node*, std::unordered_map<std::string, int>&&);
    Tree();



    virtual ~Tree();
    void traverse_tree();
    void destroy();
    void print_node(Node&);
private:
    void traverse(Node*);




};




#endif //SC_DNA_TREE_H
