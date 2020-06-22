/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "../Parallel_Random_Number_Generator/random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ie, ik;
double mass;
double walker[m_props];
const int m_bins = 1000;
int nbins, histo[m_bins];
double x_min, x_max, bin_size;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_tot,stima_kin,err_pot,err_tot,err_kin;

//configuration
double x;

// simulation
int nstep, nblk;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Print(void);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void Measure(void);
double Prob(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
