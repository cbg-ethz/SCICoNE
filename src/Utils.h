//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#ifndef SC_DNA_UTILS_H

#include "xxhash.h"
#include <iostream>
#include <map>

#define SC_DNA_UTILS_H


class Utils {

public:

    static uint64_t calcul_hash(const void* buffer, size_t length);
    static bool is_empty_map(std::map<u_int, int>& dict);
};


#endif //SC_DNA_UTILS_H
