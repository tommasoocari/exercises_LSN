/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <string>

//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,igofr;
double nbins;
double walker[m_props];

// averages
double acc,att;
double histo[m_props];
double bin_size;

double blk_av[m_props], blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot, stima_kin, stima_etot, stima_temp;
double err_pot,err_etot,err_kin,err_temp,err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint, seed;
double delta;
bool start_with_old;
bool rescale_temp_bool;
double desired_temp;

// error
int n_blocks;

//print
std::string namefile;

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
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
