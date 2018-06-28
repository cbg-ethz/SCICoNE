//
// Created by Tuncel  Mustafa Anil on 6/28/18.
//

#ifndef SC_DNA_TREE_H
#define SC_DNA_TREE_H

#include <iostream>

template <class T>
class TreeNode {
public:
    T value;
    bool isLeaf();
    TreeNode* parent;
    TreeNode* firstChild;
    TreeNode* rightSibling;
    void insertChild(TreeNode*);
    void insertNextSibling(TreeNode*);
    void removeFirstChild();
    void removeNextSibling();

    TreeNode(T val);
    ~TreeNode();

    static void traverse(TreeNode* root) {
        if (root->isLeaf())
            std::cout << "leaf: ";
        else
            std::cout << "internal: ";
        std::cout<<root->value << std::endl;
        for (TreeNode* temp = root->firstChild; temp != nullptr; temp=temp->rightSibling)
        {
            traverse(temp);
        }
    }

};



#endif //SC_DNA_TREE_H
