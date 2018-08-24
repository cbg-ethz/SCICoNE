//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#ifndef SC_DNA_UTILS_H

#include "xxhash.h"
#include <iostream>

#define SC_DNA_UTILS_H


class Utils {

public:

    static uint64_t calcul_hash(const void* buffer, size_t length);
};


#endif //SC_DNA_UTILS_H
