/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __VAR__
#define __VAR__



//functions
double Get_Energy(double mu, double sigma);
void Input1(double, double);
void Reset(int);
void Print(void);
void Accumulate(void);
double Averages(int);
void Move(double, double);
void ConfFinal(void);
void Measure(double, double);
double Prob(double, double, double);
double Error(double,double,int);

double Gauss(double x, double mu, double sigma);

double Wavefunction(double x, double mu, double sigma);
double Wavefunction_prime(double x, double mu, double sigma);
double Wavefunction_second(double x, double mu, double sigma);

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
