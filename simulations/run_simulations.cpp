//
// Created by Tuncel  Mustafa Anil on 9/27/18.
//

#include <iostream>
#include "Simulation.h"

using namespace std;



int main(int argc, char* argv[])
{

    int n_regions = 50;
    int n_nodes = 50;
    double lambda_r = 0.1;
    double lambda_c = 0.2;
    int n_cells = 500;
    int n_reads = 10000;
    int max_region_size = 10;
    int ploidy = 2;

    // create the D matrix
    vector<vector<double>> D(n_cells, vector<double>(n_regions)); // initialize with the default value

    vector<int> region_sizes(n_regions); // sampling the region sizes

    Simulation sim(n_regions, n_nodes, lambda_r, lambda_c, n_cells, n_reads, max_region_size, ploidy);


    cout <<'d';
    //simulate(n_regions, n_nodes, lambda_r, lambda_c, n_cells, n_reads, D, region_sizes); // initializes D and region_sizes

    // initialize the tree and infer the CNV profiles




}

