cmake_minimum_required(VERSION 3.8)
project(sc_dna)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "") # TODO: put -Wall -Wextra later to see

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS_DEBUG "-g") # icc needs -O3 after -g

include_directories(src)
include_directories(include)
include_directories(simulations)

add_executable(sc_dna
        src/MathOp.cpp
        src/main.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h simulations/Simulation.cpp simulations/Simulation.h)

add_executable(tests tests/run_tests.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h tests/validation.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h)

add_executable(simulations simulations/run_simulations.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h simulations/run_simulations.cpp simulations/Simulation.cpp simulations/Simulation.h)


#src/Inference.h
