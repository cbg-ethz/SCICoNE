//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_SINGLETONRANDOMGENERATOR_H
#define SC_DNA_SINGLETONRANDOMGENERATOR_H

#include <random>

class SingletonRandomGenerator
{
public:


    std::mt19937 generator; //Standard mersenne_twister_engine

    static std::mt19937& get_generator(int seed=-1)
    {
        static SingletonRandomGenerator instance(seed); // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance.generator;
    }
private:
    SingletonRandomGenerator(int seed) {

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
