#ifndef __TSP__
#define __TSP__

#include <vector>

//Random numbers
#include "../Parallel_Random_Number_Generator/random.h"
int seed[4];
Random rnd;

//population
const int n_indiv = 1000, n_genes = 32, n_coord = 2;
std::vector<std::vector<double>> population, coordinates;

//Simulation
int n_steps;
double probability_crossover, probability_mutation;

//parallel computing
int size_par, rank_par;

//function
void Input(void);
void Coordinates(void);
void New_Population(void);
void Normalize(void);
void Measure(int);
void Next_Generation(double,double);
void Write_Best_Path(void);
void Migration(void);

//reproduction
std::vector<std::vector<double>> Cross_Over(int,int);
int Selection(double,double);

//others
std::vector<double> Random_Indiv(int);
double Distance(std::vector<double>, std::vector<double>);
double Fitness(std::vector<double>);

//comparisons
bool Compare_Indiv(std::vector<double> v1, std::vector<double> v2);

#endif
