#include "lib.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>


using namespace std;

void do_sum( vector<int> vect_N, vector<double> random, string filename){

  double sum =0;
  int k=0;
  int L=0;
  ofstream file_out;
  int max_N = * max_element(vect_N.begin(), vect_N.end());

  file_out.open(filename);
  L = random.size() / max_N; 
  
  for(vector<int>::iterator it = vect_N.begin(); it != vect_N.end(); ++it){
    //for # elements in the sum
    
    for(int i=0; i<L; i++){ //for each data-block
      sum = 0;
      
      for (int j=0; j<*it; j++){ //performing the sum
	  k = j +i*max_N;
	  sum = sum+random[k];
	}

      sum = sum / (*it);
      
      file_out << sum << ",";
    }

    file_out << endl;
  }

  file_out.close();
  
}
