#include <cmath>
#include "lib.h"

using namespace std;

double length(vector<double> v1, vector<double> v2){
  if(v1.size() != v2.size()) return 0;
  double su2=0;
  for(int i=0 ; i<v1.size(); i++) su2 = su2 + pow(v1[i] - v2[i],2);
  return sqrt(su2);
}
