#include <iostream>
#include "../Parallel_Random_Number_Generator/random.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;

double error(vector<double> AV, vector<double> AV2, int n);

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


  // esercizio 01.1
  int M = 100000;
  int N = 100;
  int L = (M/N);
  vector<double> r;
  vector<double> ave(N);
  vector<double> av2(N);
  vector<double> sum_prog(N);
  vector<double> su2_prog(N);
  vector<double> err_prog(N);
  double sum=0;
  int k=0;

  ofstream file_out;

  for(int i=0; i<M; i++) r.push_back(rnd.Rannyu());
  
  //esercizio 01.1_1
  file_out.open("01.1_1_results.dat");

  for(int i=0; i<N; i++){
    sum = 0;
    for(int j=0; j<L; j++){
      k = j+ i*L;
      sum += r[k];
      }
    ave[i] = (sum/L);
    av2[i] = (pow(ave[i],2));
  }

  for(int i=0; i<N; i++){
    for(int j =0; j<i+1; j++){
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }
    sum_prog [i] = sum_prog [i]/ (i+1);
    su2_prog [i] = su2_prog [i]/ (i+1);
    err_prog[i] = error(sum_prog, su2_prog, i);

    file_out << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
  }

  file_out.close();

  
  //esercizio 01.1_2
  fill(ave.begin(), ave.end(), 0);
  fill(av2.begin(), av2.end(), 0);
  fill(sum_prog.begin(), sum_prog.end(), 0);
  fill(su2_prog.begin(), su2_prog.end(), 0);
  fill(err_prog.begin(), err_prog.end(), 0);
  
  file_out.open("01.1_2_results.dat");

  for(int i=0; i<N; i++){
    sum = 0;
    for(int j=0; j<L; j++){
      k = j+ i*L;
      sum = sum +pow((r[k] - 0.5),2);
      }
    ave[i] = (sum/L);
    av2[i] = (pow(ave[i],2));
  }

  for(int i=0; i<N; i++){
    for(int j =0; j<i+1; j++){
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }
    sum_prog [i] = sum_prog [i]/ (i+1);
    su2_prog [i] = su2_prog [i]/ (i+1);
    err_prog[i] = error(sum_prog, su2_prog, i);

    file_out << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
  }

  file_out.close();

  //esercizio 01.1_3  
  int M1=100;
  int n=10000;
  vector<int> occupation(M1,0);
  int n_cell;
  double chi_square=0;
  double expec = ((double) n) /M1;
  
  for(int i=M-1; i<M1*n; i++) r.push_back(rnd.Rannyu()); //fill random vector
  
  file_out.open("01.1_3_results.dat");

  for(int j=0; j<100; j++){

    fill(occupation.begin(), occupation.end(), 0); //reset vector of occupations
    
    for(int i=0; i<n; i++){ //fill vector of occupations
      k = i + j*n;
      n_cell = (int) (r[k] * M1);
      occupation[n_cell] = occupation[n_cell] +1;
      }

    chi_square = 0; 
    
    for(int i=0; i<M1; i++) //Chi square test
      chi_square = chi_square + pow(occupation[i] - expec,2)/expec;

    file_out << j+1 << " " << chi_square <<endl;
    
  }

  file_out.close();

  rnd.SaveSeed();
  
  return 0;
  
}


double error(vector<double> AV, vector<double> AV2, int n){
  if(n == 0) return 0;
  else return sqrt((AV2[n] - pow(AV[n],2))/n);
}


