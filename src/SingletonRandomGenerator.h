//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_SINGLETONRANDOMGENERATOR_H
#define SC_DNA_SINGLETONRANDOMGENERATOR_H

#include <random>

class SingletonRandomGenerator
{
    /*
     *  Singleton random generator used throughout the program except for the xxhash part.
     *  xxhash requires another seed than the seed used here.
     * */


public:


    std::mt19937 generator; //Standard mersenne_twister_engine

    static SingletonRandomGenerator & get_instance(int seed = -1)
    {
        static SingletonRandomGenerator instance(seed); // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }
private:
    explicit SingletonRandomGenerator(int seed) {

        if (seed == -1)
        {
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            generator = gen;
        }
        else
        {
            std::mt19937 gen(seed);
            generator = gen;
        }

    }

public:
    SingletonRandomGenerator(SingletonRandomGenerator const&) = delete;
    void operator=(SingletonRandomGenerator const&)  = delete;
};


#endif //SC_DNA_SINGLETONRANDOMGENERATOR_H
