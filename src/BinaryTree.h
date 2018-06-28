//
// Created by Tuncel  Mustafa Anil on 6/23/18.
//

#ifndef SC_DNA_BINARYTREE_H
#define SC_DNA_BINARYTREE_H

// TreeNode node definition
template <class T>
struct BinaryTreeNode
{
    T data;
    BinaryTreeNode<T> *llink;
    BinaryTreeNode<T> *rlink;
};

template <class T>
class BinaryTree {

public:
    BinaryTree() {};
    virtual ~BinaryTree() {};

private:

};


#endif //SC_DNA_BINARYTREE_H
