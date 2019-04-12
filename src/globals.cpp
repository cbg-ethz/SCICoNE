//
// Created by Tuncel  Mustafa Anil on 10/19/18.
//

#ifndef SC_DNA_GLOBALS_H
#define SC_DNA_GLOBALS_H

#include <map>
#include <string>
using namespace std;
extern int print_precision; // the precision for the printed double/float values
extern int copy_number_limit; // the limit on max copy number state a region/bin can have

extern double v; // v param from the manuscript, for the size changing moves
extern double lambda_s;
extern double lambda_r; // lambda param for the poisson that generates the number of regions
extern double lambda_c; // lambda param for the poisson that generates the copy number state of a region
extern double c_penalise; // the penalisation term for cancelling events
extern int verbosity; // the verbosity of the output
extern string f_name_postfix; // postfix for the output files created
#endif //SC_DNA_GLOBALS_H
