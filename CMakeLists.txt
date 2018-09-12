cmake_minimum_required(VERSION 3.10)
project(sc_dna)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(include)

add_executable(sc_dna
        src/main.cpp
        src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h /usr/local/Cellar/gperftools/2.7/include/gperftools/profiler.h)

add_executable(tests tests/run_tests.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h tests/validation.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h)


#src/Inference.h