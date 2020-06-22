#include "distribution.h"

//2p hydrogen orbital

double prob(std::vector<double> x){
  std::vector<double> origin(3,0);
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2] * x[2]);
  return x[2]*x[2] * 1/(32*M_PI) * exp(-r);
};


