#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include "../lib/error.h"
#include <vector>
#include <cmath>
#include <fstream>

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
  
  //exercise
  int Nblocks;
  int M;

  double d; //side of the simulation space
  double L; //length of the stick


  //Read input informations
  ifstream ReadInput;
    
  ReadInput.open("input.dat");

  ReadInput >> Nblocks;
  ReadInput >> M;
    
  ReadInput >> d;
  ReadInput >> L;

  
  cout << "Simulation of Buffon's experiment" << endl;

  cout << endl << "Side of the simualtion space = " << d << endl;
  cout << "Length of the stick = " << L << endl;

  cout << "Number of blocks = " << Nblocks << endl;
  cout << "Number of steps in one block = " << M << endl << endl;
  ReadInput.close();
  
  
  ofstream file_out;
  file_out.open("01.3_montecarlo.csv");

  int N_tot = M*Nblocks;
  vector<double> extr1(2,0);
  vector<double> extr2(2,0);
  bool hit=false;

  vector<double*> r_point(N_tot,0);
  vector<double*> r_angle;
 
  for(int i=0; i<N_tot; i++){ //fill vector coordinates of the center
    r_point[i] = new double[2];
    r_point[i][0] = rnd.Rannyu(0,d);
    r_point[i][1] = rnd.Rannyu(0,d);
  }
  
  for(int i=0; i<N_tot; i++) //fill vector of angles 
    r_angle.push_back(rnd.Uniform_Angle());

  //Buffon experiment
  int N_hit=0;
  int k=0;

  for(int block = 0; block < Nblocks; block++){ //for each block

    N_hit = 0;
    
    for(int i=0; i<M; i++){

      k = i + block*M;
      
      for(int q=0; q<2; q++){
	extr1[q] = r_point[k][q] + (L/2) * r_angle[k][q];
	extr2[q] = r_point[k][q] - (L/2) * r_angle[k][q];
      }

      hit = !(extr1[1] > 0 and extr1[1] < d and extr2[1] > 0 and extr2[1] < d);
      if(hit) N_hit += 1;
    }

    file_out << 2*L*M*1.0/(d*N_hit) << endl;

  }

  //blocking method for uncertainties
  do_blocking_method("01.3_montecarlo.csv","01.3_results.csv");

  file_out.close();
  
  rnd.SaveSeed();
  
  return 0;
  
} 
