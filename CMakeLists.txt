cmake_minimum_required(VERSION 3.10)
project(sc_dna)

# using cpp14
set(CMAKE_CXX_STANDARD 14)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

add_executable(sc_dna src/main.cpp)
target_link_libraries(sc_dna ${CONAN_LIBS})
