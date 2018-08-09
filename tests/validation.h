//
// Created by Tuncel  Mustafa Anil on 8/9/18.
//

#ifndef SC_DNA_VALIDATION_H
#define SC_DNA_VALIDATION_H


#include "Inference.h"
#include <vector>

using namespace std;




void test_swap_label()
{
    // counts per region per cell
    vector<vector<int>> D = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

    // region sizes
    vector<int> r = {4,2,3,5,2};

    // move probabilities
    vector<float> move_probs = {0.0f,0.0f,1.0f,0.0f};

    Inference mcmc;
    mcmc.initialize_worked_example();
    mcmc.compute_t_table(D,r);

    mcmc.apply_swap(D,r,false,true);

    assert(int(mcmc.t_prime_sums[0]) == -552);
    assert(int(mcmc.t_prime_sums[1]) == -413);
    assert(int(mcmc.t_prime_sums[2]) == -670);
    assert(int(mcmc.t_prime_sums[3]) == -547);
    assert(int(mcmc.t_prime_sums[4]) == -406);

    cout<<"Swap label validation test passed!"<<endl;


}



#endif //SC_DNA_VALIDATION_H
