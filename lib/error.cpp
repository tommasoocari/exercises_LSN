#include "error.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iterator>

using namespace std;


double average(std::vector<double>::iterator begin, std::vector<double>::iterator end){
  return accumulate(begin, end, 0.0)/distance(begin,end);
}

double error(std::vector<double>::iterator begin, std::vector<double>::iterator end){
  int N = distance(begin,end);
  double sum = 0;
  for(auto it = begin; it != end; it = it+1) sum += pow(*it,2);
  if(N == 1) return 0;
  else return sqrt( (sum/N - pow(average(begin,end),2))/(N-1) );
}

double average(vector<double> v){ return average(v.begin(), v.end());}

double error(vector<double> v){ return error(v.begin(), v.end());}




void do_blocking_method(string filename_in, string filename_out){

  cout << endl << "Blocking method to calculate uncertainties" << endl;
  
  ifstream file_in;
  ofstream file_out;
  double n;
  vector<double> vec;
  vector<double> ve2;
  double prog_sum;
  double prog_su2;
  double error;
  
  file_in.open(filename_in);

  while(!file_in.eof()){
    file_in >> n;
    vec.push_back(n);
    ve2.push_back(n*n);
  };

  vec.erase(vec.end()-1);
  ve2.erase(ve2.end()-1);
  
  file_in.close();
  
  file_out.open(filename_out);

  for(int i=0; i<vec.size(); i++){
    prog_sum=0;
    prog_su2=0;
    for(int j=0; j<i+1; j++){
      prog_sum = prog_sum + vec[j];
      prog_su2 = prog_su2 + ve2[j];
    }
    prog_sum = prog_sum /(i+1);
    prog_su2 = prog_su2 /(i+1);
    if(i==0) error = 0;
    else error = sqrt((prog_su2 - prog_sum*prog_sum)/i);
    file_out << i+1 << "," << prog_sum << "," << error << "," << endl; 
  }

  cout << "Results written in " << filename_out << endl;  
  
  file_out.close();
}


void do_raw_blocking_method(string filename_in, string filename_out, int n_blocks){

  cout << endl << "Blocking method to calculate uncertainties" << endl;
  
  vector<double> raw_data;
  ifstream file_in;
  ofstream file_out;
  int steps_per_block;
  double n,appo;
  vector<double> vec;
  vector<double> ve2;
  double prog_sum;
  double prog_su2;
  double error;
  
  file_in.open(filename_in);

  while(!file_in.eof()){
    file_in >> n;
    raw_data.push_back(n);
  };

  raw_data.erase(raw_data.end()-1);

  steps_per_block = raw_data.size()/n_blocks;

  for(int i=0; i<n_blocks; i++){
    appo = average(raw_data.begin()+i*steps_per_block, raw_data.begin()+(i+1)*steps_per_block);
    vec.push_back(appo);
    ve2.push_back(appo*appo);
  }

  //vec.erase(vec.end()-1);
  //ve2.erase(ve2.end()-1);
  
  file_in.close();
  
  file_out.open(filename_out);

  for(int i=0; i<vec.size(); i++){
    prog_sum=0;
    prog_su2=0;
    for(int j=0; j<i+1; j++){
      prog_sum = prog_sum + vec[j];
      prog_su2 = prog_su2 + ve2[j];
    }
    prog_sum = prog_sum /(i+1);
    
    prog_su2 = prog_su2 /(i+1);
    if(i==0) error = 0;
    else error = sqrt((prog_su2 - prog_sum*prog_sum)/i);
    file_out << i+1 << "," << prog_sum << "," << error << "," << endl;  
  }

  cout << "Results written in " << filename_out << endl;
  
  file_out.close();
  
}
