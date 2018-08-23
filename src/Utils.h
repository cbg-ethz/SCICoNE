//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#ifndef SC_DNA_UTILS_H

#include "xxhash.h"

#define SC_DNA_UTILS_H


class Utils {

public:

    static unsigned long long calcul_hash(const void* buffer, size_t length);
};


#endif //SC_DNA_UTILS_H
