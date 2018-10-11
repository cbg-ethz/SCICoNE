//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#ifndef SC_DNA_UTILS_H

#include "xxhash.h"
#include <iostream>
#include <map>
#include <random>
#include "SingletonRandomGenerator.h"
#include "MathOp.h"
#include <unistd.h>
#include <fstream>
#include <sstream>

using namespace std;

#define SC_DNA_UTILS_H


class Utils {

public:

    static uint64_t calcul_hash(const void* buffer, size_t length);
    static bool is_empty_map(std::map<u_int, int>& dict);
    static void random_initialize_labels_map(std::map<u_int, int> &distinct_regions, int n_regions, double lambda_r,
                                             double lambda_c);

    static void read_counts(vector<vector<double>> &mat, const string path);
};


#endif //SC_DNA_UTILS_H
