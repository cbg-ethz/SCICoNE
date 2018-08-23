//
// Created by Tuncel  Mustafa Anil on 8/23/18.
//

#include "Utils.h"

unsigned long long Utils::calcul_hash(const void *buffer, size_t length) {

    /*
     * Calculates the hash using xxhash.
     * */

    unsigned long long const seed = 0;   /* or any other value */
    unsigned long long const hash = XXH64(buffer, length, seed);
    return hash;

}
