//
// Created by Tuncel  Mustafa Anil on 6/28/18.
//

#include "TreeNode.h"


template<class T>
bool TreeNode<T>::isLeaf() {
    if (this->firstChild == nullptr)
        return true;
    else
        return false;
}

template<class T>
void TreeNode<T>::insertChild(TreeNode *node) {
    if (this->firstChild == nullptr)
    {
        this->firstChild = node;
    }
    else
    {
        for (TreeNode<T>* temp = this->firstChild; temp != nullptr; temp = temp->rightSibling) {
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

template<class T>
void TreeNode<T>::removeNextSibling() {

}

template<class T>
TreeNode<T>::TreeNode(T val)
{
    this->value = val;
    this->firstChild = this->parent = this->rightSibling= nullptr;
}

template<class T>
TreeNode<T>::~TreeNode() {

}



//
template TreeNode<int>::TreeNode(int);
template TreeNode<int>::~TreeNode();
template void TreeNode<int>::insertNextSibling(TreeNode*);
template void TreeNode<int>::insertChild(TreeNode*);

template bool TreeNode<int>::isLeaf();