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

bool Utils::is_empty_map(std::map<u_int, int> &dict){
    /*
     * Returns true if the map is empty, e.g. all zero
     * */
    for (auto const &it : dict)
    {
        if(it.second != 0)
            return false;
    }
    return true;
}
