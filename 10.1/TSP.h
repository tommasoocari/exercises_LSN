#ifndef __TSP__
#define __TSP__

#include <vector>

//Random numbers
#include "../Parallel_Random_Number_Generator/random.h"
int seed[4];
Random rnd;

//properties
const int n_cities = 32;
std::vector<double> path;
std::vector<std::vector<double>> annealing_schedule;
std::vector<std::vector<double>> coordinates;

//simulation
int n_steps, n_blocks;
int n_print = n_steps / 100;
double accepted, attempted;
int istep;
double beta, n_steps_block;

//functions
void Input(void);
void Coordinates(void);
void Normalize(void);
void Step(double);
void Reset(int);
void Measure(int);
void End_Block(int,double,double);
void Write_Best_Path(void);

std::vector<double> RandomPath(int);
double Distance(std::vector<double>, std::vector<double>);
double Fitness(std::vector<double>);

#endif
