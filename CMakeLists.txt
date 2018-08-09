cmake_minimum_required(VERSION 3.10)
project(sc_dna)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

add_executable(sc_dna
        src/main.cpp
        src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h tests/validation.h)
#src/Inference.h