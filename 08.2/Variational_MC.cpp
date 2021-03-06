/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Variational_MC.h"
#include "Wavefunction.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Print();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;

  cout << "Variational Quantum Monte Carlo simulation (one particle)" << endl << endl;
  cout << "External 1D potential v(x) = x^4 - 2.5 * x^2" << endl << endl;  
  cout << "The program uses natural units " << endl;
  
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> mass;

  cout << "Particle mass = " << mass << endl;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
 
  n_props = 3; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("config.0");
  ReadConf >> x;
  cout << "Initial position = " << x << endl;
  ReadConf.close();

//Preparing histogram
  nbins = 100;
  x_min = -5;
  x_max = 5;
  bin_size = (x_max - x_min) / ((double) nbins);
  
  for(int i=0; i<nbins; i++) histo[i] = 0;
  
//Evaluate energy of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial kinetic energy = " << walker[ik] << endl;
  cout << "Initial potential energy = " << walker[iv] << endl;
  cout << "Initial total energy = " << walker[ie] << endl << endl;

}

void Move(void)
{
  double p, prob_old, prob_new;
  double xold, xnew;

  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)

  //Old
  xold = x;

  prob_old = Prob(xold);

  //New
  xnew = x + delta*(rnd.Rannyu() - 0.5);
  
  prob_new = Prob(xnew);

  //Metropolis test
  if(prob_old == 0) p = 1;
  else p = prob_new/prob_old;
  
  if(p >= rnd.Rannyu())  
  {
  //Update
    x = xnew;
    accepted = accepted + 1.0;
  }
  attempted = attempted + 1.0;

}

void Measure()
{
  int bin;
  
  walker[iv] = pow(x,4) - 2.5*pow(x,2);   
  walker[ik] = - 1.0/(2 * mass) * Wavefunction_second(x) / Wavefunction (x);
  walker[ie] = walker[ik] + walker[iv];

  if(x > x_min and x < x_max)
  {
    bin = (int)((x - x_min) / bin_size);
    histo[bin] = histo[bin] + 1;
  }
  
  
}

double Prob(double x)
{
  //not normalized!

  return pow(Wavefunction(x),2);
  
}


void Print(void)
{
  ofstream Epot, Etot, Ekin, Points;
    
   Epot.open("all.epot.0",ios::app);
   Ekin.open("all.ekin.0",ios::app);
   Etot.open("all.etot.0",ios::app);
   Points.open("all.points.0",ios::app);

   Epot << walker[iv] << endl;
   Ekin << walker[ik] << endl;
   Etot << walker[ie] << endl;
   Points << x << endl;
   
   Epot.close();
   Ekin.close();
   Etot.close();
   Points.close();
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
   
  ofstream Epot, Etot, Ekin;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output.epot.0",ios::app);
    Ekin.open("output.ekin.0",ios::app);
    Etot.open("output.etot.0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_kin = blk_av[ik]/blk_norm; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_tot = blk_av[ie]/blk_norm; //Total energy
    glob_av[ie] += stima_tot; 
    glob_av2[ie] += stima_tot*stima_tot;
    err_tot=Error(glob_av[ie],glob_av2[ie],iblk);

    //if(iblk == nblk){

    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
    Etot << setw(wd) << iblk <<  setw(wd) << stima_tot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_tot << endl;

    //}

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  WriteConf << x << endl;
  
  WriteConf.close();

  rnd.SaveSeed();

}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
    
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
