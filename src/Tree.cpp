//
// Created by Tuncel  Mustafa Anil on 7/4/18.
//

#include "Tree.h"

Tree::Tree(u_int ploidy)
{
    root = new Node();
    this->ploidy = ploidy;
    // creates a copy of the root ptr and stores it in the vector
    // TODO std::move can be used here
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
    child->c = child->c_change;


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







