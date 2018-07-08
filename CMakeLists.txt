cmake_minimum_required(VERSION 3.10)
project(sc_dna)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

add_executable(sc_dna
        src/BinaryTree.cpp
        src/BinaryTree.h
        src/main.cpp
        src/MathOp.cpp
        src/MathOp.h
        src/Tree.cpp
        src/Tree.h
        src/TreeNode.h)
