//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#ifndef SC_DNA_UTILS_H

#include "xxhash.h"
#include <iostream>
#include <map>
#include <set>
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
    static map<u_int, int> map_diff(map<u_int, int> a, map<u_int, int> b);
    static void random_initialize_labels_map(std::map<u_int, int> &distinct_regions, int n_regions, double lambda_r,
                                             double lambda_c);
    static void read_counts(vector<vector<double>> &mat, const string &path);
    static void read_vector(vector<int> &vec, const string &path);
    static vector<vector<double>> condense_matrix(vector<vector<double>>& D, vector<int>& region_sizes);
    static vector<vector<int>> regions_to_bins_cnvs(vector<vector<int>>& cnvs, vector<int>& region_sizes);

};


#endif //SC_DNA_UTILS_H
