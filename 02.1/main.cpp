#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include "../lib/error.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>

using namespace std;

double f1(double x) {return M_PI/2* cos( M_PI * x / 2 );} //integrand
double f2(double x) {return M_PI/4 * cos( M_PI * x / 2) / (1-x); } //integrand2

int main(){

  //Random Number Generator
  Random rnd;

  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
   
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  //exercise
  int N = 100; // number of blocks
  int M = 10000; // montecarlo steps in each block


  //Read input informations
  ifstream ReadInput;
    
  ReadInput.open("input.dat");

  ReadInput >> N;
  ReadInput >> M;

  
  cout << "Monte Carlo integration" << endl;

  cout << endl << "Number of blocks = " << N << endl;
  cout << "Number of steps in one block = " << M << endl << endl;

  cout << endl
       << "The calculated integral is int_0^1 (pi / 2) * cos(pi * x / 2)"
       << endl;
  
  ReadInput.close();

  
  ofstream file_out;
  int L = N * M;  

  
  vector<double> r;
  vector<double> f;
  for(int i=0; i<L; i++){
    r.push_back(rnd.Rannyu());
    f.push_back(f1(r[i]));
  }

  
  file_out.open("02.1_montecarlo_ud.csv");

  cout << endl << "Using uniform distribution" << endl; 
  
  for(int i=0; i<L; i=i+M){ //using uniform distribution and integrand (f1)
    file_out << accumulate(f.begin()+i,f.begin()+i+M,0.0) / M << endl;
   }
  
  file_out.close();
  
  do_blocking_method("02.1_montecarlo_ud.csv","02.1_results_ud.csv");


  // inverse cumulative method
  
  for(int i=0; i<r.size(); i++){
    r[i] = 1+ sqrt(1-r[i]);
    f[i] = f2(r[i]);
  }

  file_out.open("02.1_montecarlo_is.csv");

  cout << endl << "Using importance sampling: p(x) = 1 - x" << endl; 
  
  for(int i=0; i<L; i=i+M){ //using p(x) = 1 - x and integrand2 (f2) 
    file_out << accumulate(f.begin()+i,f.begin()+i+M,0.0) / M << endl;
   }
  
  file_out.close();
  
  do_blocking_method("02.1_montecarlo_is.csv","02.1_results_is.csv");

  rnd.SaveSeed();
  
  return 0;
  
} 
