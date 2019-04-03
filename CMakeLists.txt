cmake_minimum_required(VERSION 3.6)
project(sc_dna)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG  "-g") # icc needs -O3 after -g
# set(CMAKE_CXX_FLAGS_DEBUG "-Wpedantic")
# set(CMAKE_CXX_FLAGS_DEBUG "-Wall")
# set(CMAKE_CXX_FLAGS_DEBUG "-Wextra")



# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for boost
find_package(Boost 1.47 COMPONENTS random REQUIRED)



include_directories(src)
include_directories(include)
include_directories(simulations)

add_executable(inference
        src/MathOp.cpp
        src/infer_trees.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h)

add_executable(tests tests/run_tests.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h tests/validation.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h)

add_executable(simulation simulations/run_simulations.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h simulations/run_simulations.cpp simulations/Simulation.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h)

add_executable(score src/score_tree.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h)

add_executable(breakpoint_detection
        src/breakpoint_detection.cpp src/MathOp.cpp src/MathOp.h  src/Utils.cpp src/Utils.h  src/SignalProcessing.cpp src/SignalProcessing.h src/Lgamma.cpp src/Lgamma.h include/xxhash.h include/xxhash.c )

#src/Inference.h
