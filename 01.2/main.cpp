#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include "lib.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

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

  //esercizio 1.2b
  int L = 10000;
  vector<int> vect_N = {1,2,10,100};

  int max_N = * max_element(vect_N.begin(), vect_N.end());

  vector<double> r_u;
  vector<double> r_exp;
  vector<double> r_cauchy;

  //generation
  for(int i=0; i<max_N*L; i++){
    r_u.push_back(rnd.Rannyu());
    r_exp.push_back(rnd.Rannyexp(1));
    r_cauchy.push_back(rnd.Rannycauchy(0,1));
  };

  //accumuation
  do_sum(vect_N, r_u, "01.2_results_u.csv"); 
  do_sum(vect_N, r_exp, "01.2_results_exp.csv");
  do_sum(vect_N, r_cauchy, "01.2_results_cauchy.csv");
  
  rnd.SaveSeed();
  
  return 0;
  
}


