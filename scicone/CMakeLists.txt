cmake_minimum_required(VERSION 3.9)
project(sc_dna)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-fsanitize=undefined")
set(CMAKE_CXX_FLAGS "-fsanitize=address")

# set(CMAKE_CXX_FLAGS "--coverage")
# set(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS_DEBUG  "-g2") # icc needs -O3 after -g
set(CMAKE_CXX_FLAGS_DEBUG "-Wuninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "-Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-Wextra")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")

# Include test data
configure_file(Config.h.in ${CMAKE_BINARY_DIR}/Config.h)
include_directories(${CMAKE_BINARY_DIR})

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for OpenMP for parallelization
find_package(OpenMP)

if (OPENMP_FOUND)
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

# Search for boost and nlopt
find_package(Boost 1.47 COMPONENTS random REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(NLopt REQUIRED)
include_directories(${NLOPT_INCLUDE_DIRS})

include_directories(src)
include_directories(include)
include_directories(simulations)

add_executable(inference
        src/MathOp.cpp
        src/infer_trees.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h src/CustomExceptions.h)

add_executable(simulation simulations/run_simulations.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h simulations/run_simulations.cpp simulations/Simulation.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h src/CustomExceptions.h)

add_executable(score src/score_tree.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h src/CustomExceptions.h)

add_executable(breakpoint_detection
        src/breakpoint_detection.cpp src/MathOp.cpp src/MathOp.h  src/Utils.cpp src/Utils.h  src/SignalProcessing.cpp src/SignalProcessing.h src/Lgamma.cpp src/Lgamma.h include/xxhash.h include/xxhash.c src/CustomExceptions.h)

# tests
enable_testing()

add_executable(tests tests/run_tests.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h  src/SignalProcessing.cpp src/SignalProcessing.h  src/SingletonRandomGenerator.h tests/validation.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h src/CustomExceptions.h)

add_executable(save_root_per_nu tests/save_root_score_per_nu.cpp src/MathOp.cpp
        src/MathOp.h
        src/Tree.h src/Node.h  src/Inference.h src/SingletonRandomGenerator.h tests/validation.h include/xxhash.h include/xxhash.c src/Utils.cpp src/Utils.h src/globals.cpp src/Lgamma.cpp src/Lgamma.h src/CustomExceptions.h)

# unit tests
add_test(unit_tests tests)

if (OPENMP_FOUND)
  target_link_libraries(inference NLopt::nlopt OpenMP::OpenMP_CXX)
  target_link_libraries(simulation NLopt::nlopt OpenMP::OpenMP_CXX)
  target_link_libraries(score NLopt::nlopt OpenMP::OpenMP_CXX)
  target_link_libraries(breakpoint_detection NLopt::nlopt OpenMP::OpenMP_CXX)
  target_link_libraries(tests NLopt::nlopt OpenMP::OpenMP_CXX)
  target_link_libraries(save_root_per_nu NLopt::nlopt OpenMP::OpenMP_CXX)
else()
  target_link_libraries(inference NLopt::nlopt)
  target_link_libraries(simulation NLopt::nlopt)
  target_link_libraries(score NLopt::nlopt)
  target_link_libraries(breakpoint_detection NLopt::nlopt)
  target_link_libraries(tests NLopt::nlopt)
  target_link_libraries(save_root_per_nu NLopt::nlopt)
endif(OPENMP_FOUND)
