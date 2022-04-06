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

extern double lambda_s;
extern double lambda_r; // lambda param for the poisson that generates the number of regions
extern double lambda_c; // lambda param for the poisson that generates the copy number state of a region
extern double cf; // the cluster fraction variable to be used in tree prior
extern double c_penalise; // the penalisation term for cancelling events
extern unsigned is_overdispersed; // param to specify if tree scoring will use the overdispersed model
extern int verbosity; // the verbosity of the output
extern string f_name_posfix; // posfix for the output files created
extern double eta; // a small value to replace 0 when needed, for log(0) case for example
extern double beta; // penalty for number of genotypes
#endif //SC_DNA_GLOBALS_H
