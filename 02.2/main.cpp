#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include "../lib/error.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

vector<int> discrete_step(int n);

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


  int N_blocks = 100;
  int N_simul = 100;
  int L = 100;

  //Read input informations
  ifstream ReadInput;
    
  ReadInput.open("input.dat");

  ReadInput >> N_blocks;
  ReadInput >> N_simul;
  ReadInput >> L; //number of RW steps
  
  cout << "RW Simulations" << endl;

  cout << endl << "Number of RW steps =" << L << endl;
  cout << "Number of blocks = " << N_blocks << endl;
  cout << "Number of Monte Carlo steps in one block = "
       << N_simul << endl << endl;
  
  ReadInput.close();

  int N_perstep = N_blocks*N_simul;
  int M = N_simul * N_blocks * L;


  //Discrete RW
  vector<int> r; //vector of random numbers
  for(int i=0; i<M; i++) r.push_back(rnd.Rannyu(0,6));

  vector<vector<int>> positions; //vector positions of the walker
  for (int i=0; i<N_perstep; i++) positions.push_back({0,0,0});

  vector<int> appo;
  vector<double> blocks_results(N_blocks,0);
  
  double sum=0;
  int pos=0;
  int k=0;

  ofstream file_out;
  file_out.open("02.2_results_discrete.csv");

  cout << endl << "Discete RW " << endl;
  cout << "Uniform probability of step of length 1 along an axis" << endl;
  
  for(int step=0; step<L; step++){ // for RW steps

    pos = 0;

    for(int block = 0; block<N_blocks; block++){ // for blocks

      sum = 0;

      for(int i=0; i<N_simul; i++){ //simulation
	k = i + block * N_simul + step * N_blocks * N_simul;
	pos = i + block * N_simul;
	appo = discrete_step(r[k]);
	for(int d=0; d<3; d++) positions[pos][d] += appo[d]; //update positions
	sum += pow(positions[pos][0],2)+pow(positions[pos][1],2)
	  +pow(positions[pos][2],2); //calculate distance
      }

      blocks_results[block] = sqrt(sum/N_simul); //average distance of a block

    }

    file_out << step+1 << ","
	     << average(blocks_results) << "," //block average
	     << error(blocks_results) << endl; //block MSD for error
    
  }
  
  file_out.close();

  
  //continuum RW
  vector<double*> r2; //vector of uniform solid angles
  for(int i=0; i<M; i++) r2.push_back(rnd.Uniform_Solid_Angle());

  vector<vector<double>> positions2; //vector positions of N_perstep walkers
  for (int i=0; i<N_perstep; i++) positions2.push_back({0.,0.,0.});

  sum=0;
  pos =0;
  k=0;

  double * appo2;

  file_out.open("02.2_results_continuum.csv");

  cout << endl << "Continuum RW " << endl;
  cout << "Steps of length 1 choosing an uniform solid angle" << endl;
  
  for(int step=0; step<L; step++){// for RW step

    pos = 0; 

    for(int block = 0; block<N_blocks; block++){ // for block

      sum = 0;

      for(int i=0; i<N_simul; i++){ // simulation
	k = i + block * N_simul + step * N_blocks * N_simul;
	pos = i + block * N_simul;
	appo2 = r2[k];
	for(int d=0; d<3; d++) positions2[pos][d] += appo2[d]; //update position
	sum += pow(positions2[pos][0],2)+pow(positions2[pos][1],2)
	  +pow(positions2[pos][2],2); //calculate distance
      }
      

      blocks_results[block] = sqrt(sum/N_simul); //average distance of a block

    }

    file_out << step+1 << ","
	     << average(blocks_results) << "," //block average
	     << error(blocks_results) << endl; //block MSD for error
  }
  
  file_out.close();
  
  rnd.SaveSeed();
  
  return 0;

}

vector<int> discrete_step(int n){
  switch (n){
  case 0:
    return {1,0,0};
  case 1:
    return {-1,0,0};
  case 2:
    return {0,1,0};
  case 3:
    return {0,-1,0};
  case 4:
    return {0,0,1};
  case 5:
    return {0,0,-1};
  default:
    return {0,0,0};
  }
}
    
      
      
      

      

  
  
