//
// Created by Tuncel  Mustafa Anil on 6/28/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>

//TODO have a random insert method & uniform random delete method
// a uniform random select method and a delete Node method

struct Node{
    unordered_map<std::string,int> c;
    double log_score;
    Node* first_child;
    Node* next;
};

class TreeNode {
public:
    Node* root;
    std::vector<Node*> all_nodes;
    bool isLeaf();

    /*
    void insertChild(TreeNode*);
    void insertNextSibling(TreeNode*);
    void removeFirstChild();
    void removeNextSibling();
    */

    TreeNode();
    ~TreeNode();
    void traverse(Node*);


};

void TreeNode::traverse(Node* root) {
    if (root->first_child == nullptr)
        std::cout << "leaf: ";
    else
        std::cout << "internal: ";
    std::cout<<root->log_score << std::endl;
    for (Node* temp = root->first_child; temp != nullptr; temp=temp->next)
    {
        // tail recursive
        traverse(temp);
    }
}


bool TreeNode::isLeaf() {
    if (root->first_child == nullptr)
        return true;
    else
        return false;
}
/*
void TreeNode::insertChild(unordered_map<std::string,int>& c) {
    if (root->firstChild == nullptr)
    {
        root->firstChild = node;
    }
    else
    {
        for (TreeNode* temp = this->firstChild; temp != nullptr; temp = temp->rightSibling) {
            if (temp->rightSibling == nullptr)
            {
                temp->rightSibling = node;
                break;
            }

        }
    }
}

template<class T>
void TreeNode<T>::removeFirstChild() {

    if (this->firstChild == nullptr)
        return;

    else
    {
        TreeNode<T>* temp = this->firstChild;
        this->firstChild = temp->rightSibling;
        delete temp;
        temp = nullptr;

    }

}

template<class T>
void TreeNode<T>::insertNextSibling(TreeNode* node) {

    if (this->rightSibling == nullptr)
    {
        this->rightSibling = node;
    }
    else
    {
        for (TreeNode<T>* temp = this->rightSibling; temp != nullptr; temp = temp->rightSibling) {
            if (temp->rightSibling == nullptr)
            {
                temp->rightSibling = node;
                break;
            }

        }
    }

}


void TreeNode::removeNextSibling() {

}

TreeNode::TreeNode()
{
    root = new Node;
    root->log_score = 0.0;
    root->first_child = root->next = nullptr;
}

template<class T>
TreeNode<T>::~TreeNode() {
    if (this->is_leaf()) {
        std::cout << "leaf: ";
        delete (this);
        //this = nullptr;
    } else {
        std::cout << this->value << std::endl;
        for (TreeNode *temp = this->firstChild; temp != nullptr; temp = temp->rightSibling) {
            delete temp;
        }
        delete (this);
        //this = nullptr;
    }
}

*/


#endif //SC_DNA_TREE_H
