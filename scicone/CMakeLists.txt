cmake_minimum_required(VERSION 3.16)
project(SCICoNE VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CMake 4 removed the bundled FindBoost module; Boost ships its own config.
if (POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

include(GNUInstallDirs)

# Default to an optimised build when the caller does not say otherwise.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g2 -gdwarf-3 -Wall -Wextra -Wpedantic -Wuninitialized")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")

# Sanitizers are a development aid, not something to ship: they cost an order of
# magnitude in runtime and make the binaries depend on the sanitizer runtime.
option(SCICONE_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
if (SCICONE_ENABLE_SANITIZERS)
    add_compile_options(-fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined)
endif()

# Include test data
configure_file(Config.h.in ${CMAKE_BINARY_DIR}/Config.h)
include_directories(${CMAKE_BINARY_DIR})

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for OpenMP for parallelization
find_package(OpenMP)

# Search for boost and nlopt.
# Only Boost.Random's header-only distributions are used, so the headers suffice.
# Header-only Boost installations do not all ship a CMake package configuration
# (conda's libboost-headers among them), so fall back to locating the includes.
find_package(Boost 1.47 QUIET)
if (NOT TARGET Boost::headers)
    find_path(Boost_INCLUDE_DIR NAMES boost/random/uniform_int_distribution.hpp)
    if (NOT Boost_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find the Boost headers. Install Boost >= 1.47 "
                            "or point CMAKE_PREFIX_PATH at it.")
    endif()
    add_library(Boost::headers INTERFACE IMPORTED)
    set_target_properties(Boost::headers PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}")
    message(STATUS "Found Boost headers: ${Boost_INCLUDE_DIR}")
endif()

find_package(NLopt REQUIRED)

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

# The executables end up in a shared bin/ directory once installed, so give them
# names that cannot collide with unrelated tools.
set(SCICONE_EXECUTABLES inference simulation score breakpoint_detection tests save_root_per_nu)
foreach(target IN LISTS SCICONE_EXECUTABLES)
    set_target_properties(${target} PROPERTIES OUTPUT_NAME "scicone-${target}")
    target_link_libraries(${target} PRIVATE Boost::headers NLopt::nlopt)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endforeach()

# `tests` and `save_root_per_nu` resolve their input data relative to the source
# tree (see Config.h.in), so they are only meaningful inside a build directory
# and are deliberately left out of the installation.
install(TARGETS inference simulation score breakpoint_detection
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
