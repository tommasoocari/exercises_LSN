#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include "../lib/error.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

double GBM_step(double S_prec, double delta_t, double mu,
		double sigma, double Z);

//double pp(double n) {return n+1;}

int main(){

  //Read seed for random numbers
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

  //reading parameters
  ifstream file_param;
  file_param.open("parameters.dat");

  string name_param;
  vector<double> params;
  double appo;
  
  while(file_param >> name_param >> appo) params.push_back(appo);

  file_param.close();


  //define parameters
  double S_0 = params[0], T = params[1], K = params[2];
  double r = params[3], sigma = params[4];

  int N_blocks = params[5], N_simulpb = params[6];

  int N_intervals = params[7];

  int N_simultot = N_blocks*N_simulpb;

  cout << "Program to sample the price of European options using Geometric Brownian Motion" << endl;
  cout << endl << "Asset price at t = 0 = " << S_0 << endl;
  cout << "Delivery time = " << T << endl;
  cout << "Strike price = " << K << endl;
  cout << "Risk-free interests rate = " << r << endl;
  cout << "Volatility = " << endl;

  cout << endl << "Number of blocks = " << N_blocks << endl;
  cout << "Number of simulation per block = " << N_simulpb << endl;
  cout << "Number of time intervals = " << N_intervals << endl;
  
  
  //directly sampling call options
  
  double * r1 = rnd.Gauss(0,1,N_simultot);
  vector<double> simul_results(N_simulpb,0);
  double S_t = S_0;
  double C = 0;
  int k=0;

  ofstream file_out;
  file_out.open("03.1_montecarlo_directcall.csv");
  
  for(int block=0; block<N_blocks; block++){

    for(int i=0; i<N_simulpb; i++){
      S_t = S_0;
      k = i + N_simulpb*block;
      S_t = GBM_step(S_t, T, r, sigma, r1[k]);
      C = exp(-r*T)*max(0.0,S_t-K);
      simul_results[i] = C;
    }

    file_out << average(simul_results) << endl;
    
  }

  file_out.close();

  do_blocking_method("03.1_montecarlo_directcall.csv",
		     "03.1_results_directcall.csv");

  
  //directly sampling put options

  double * r2 = rnd.Gauss(0,1,N_simultot);
  double P = 0;
  k=0;

  file_out.open("03.1_montecarlo_directput.csv");
  
  for(int block=0; block<N_blocks; block++){

    for(int i=0; i<N_simulpb; i++){
      S_t = S_0;
      k = i + N_simulpb*block;
      S_t = GBM_step(S_t, T, r, sigma, r2[k]);
      P = exp(-r*T)*max(0.0,-S_t+K);
      simul_results[i] = P;
    }

    file_out << average(simul_results) << endl;
    
  }

  file_out.close();

  do_blocking_method("03.1_montecarlo_directput.csv",
		     "03.1_results_directput.csv");


  // indirect sampling call options

  double * r3 = rnd.Gauss(0,1,N_simultot*N_intervals);
  C = 0;
  k=0;

  file_out.open("03.1_montecarlo_indirectcall.csv");
  
  for(int block=0; block<N_blocks; block++){

    for(int i=0; i<N_simulpb; i++){
      S_t = S_0;
      for(int j=0; j<N_intervals; j++){
        k = j + i*N_intervals + N_simulpb*block*N_intervals;
	S_t = GBM_step(S_t, T/N_intervals, r, sigma, r3[k]);
      }
      C = exp(-r*T)*max(0.0,S_t-K);
      simul_results[i] = C;
    }

    cout << block+1 << "%" << endl;     

    file_out << average(simul_results) << endl;

  }

  file_out.close();

  do_blocking_method("03.1_montecarlo_indirectcall.csv",
		     "03.1_results_indirectcall.csv");


  // indirectly sampling put options
  
  double * r4 = rnd.Gauss(0,1,N_simultot*N_intervals);
  P = 0;
  k=0;

  file_out.open("03.1_montecarlo_indirectput.csv");
  
  for(int block=0; block<N_blocks; block++){

    for(int i=0; i<N_simulpb; i++){
      S_t = S_0;
      for(int j=0; j<N_intervals; j++){
        k = j + i*N_intervals + N_simulpb*block*N_intervals;
	S_t = GBM_step(S_t, T/N_intervals, r, sigma, r4[k]);
      }
      P = exp(-r*T)*max(0.0,-S_t+K);
      simul_results[i] = P;      
    }

    cout << block+1 << "%" << endl; 

    file_out << average(simul_results) << endl;

  }

  file_out.close();

  do_blocking_method("03.1_montecarlo_indirectput.csv",
		     "03.1_results_indirectput.csv");

  rnd.SaveSeed();
  
  return 0;
  
}


double GBM_step(double S_prec, double delta_t, double mu,
		double sigma, double Z){
  return S_prec * exp((mu- 0.5 * sigma*sigma) * delta_t + sigma*Z* sqrt(delta_t));
}





    
      
      
      

      

  
  
