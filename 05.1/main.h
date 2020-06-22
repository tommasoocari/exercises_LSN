#ifndef __MAIN__
#define __MAIN__

#include <vector>
#include <string> 

//Random Number Generator
#include "../Parallel_Random_Number_Generator/random.h"
Random rnd;

//statistics
int N_throws;
int N_blocks;
int N_print;

//simulation
std::vector<double> pos_input(3);
double r_cutoff;
int type_T;

std::vector<double> pos_new(3);
std::vector<double> pos_old(3);
std::vector<double> origin(3,0);

//output
std::string namefile_points;
std::string namefile_data;
std::string namefile_data_acceptance;
std::string namefile_analysis;
std::string namefile_analysis_acceptance;

#endif
