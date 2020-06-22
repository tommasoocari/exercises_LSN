#ifndef __SA__
#define __SA__

#include <vector>

//Random numbers
#include "../Parallel_Random_Number_Generator/random.h"
int seed[4];
Random rnd;

//properties
double mu, sigma;
std::vector<std::vector<double>> annealing_schedule;


//simulation
int n_steps, n_blocks;
int n_print = n_steps / 100;
double accepted, attempted;
int istep;
double beta, n_steps_block;

//functions
void Input(void);
void Step(double);
void Measure(int);
void End_Block(int);


//Variational_MC.h
/*
//Random numbers
#include "../Parallel_Random_Number_Generator/random.h"
extern int seed1[4];
extern Random rnd1;

//parameters, observables
extern const int m_props=1000;
extern int n_props, iv, ie, ik;
extern double mass;
extern double walker[m_props];
extern const int m_bins = 1000;
extern int nbins, histo[m_bins];
extern double x_min, x_max, bin_size;

// averages
extern double blk_av[m_props],blk_norm,accepted1,attempted1;
extern double glob_av[m_props],glob_av2[m_props];
extern double stima_pot,stima_tot,stima_kin,err_pot,err_tot,err_kin;

//configuration
extern double x;

// simulation
extern int nstep = 10000, nblk = 1;
extern double delta;

//pigreco
extern const double pi=3.1415927;
*/
#endif

