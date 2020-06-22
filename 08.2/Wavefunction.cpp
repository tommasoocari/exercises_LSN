#include "Wavefunction.h"

double mu = 0.81;
double sigma = 0.36;

double Gauss(double x, double m_mu, double m_sigma){
  return exp(-pow(x - m_mu, 2) / (2 * m_sigma));
}

double Wavefunction(double x){
  return Gauss(x, mu, sigma) + Gauss(x, -mu, sigma);
}

double Wavefunction_prime(double x){
  return (1.0 / sigma) * ((mu - x) * Gauss(x, mu, sigma) - (mu + x) * Gauss (x, -mu, sigma));
}

double Wavefunction_second(double x){
  return (1.0 / pow(sigma, 2) ) * ( Gauss(x, mu, sigma) * (mu * mu - 2 * mu * x + x * x - sigma) + Gauss(x, -mu, sigma) * ( mu * mu + 2 * mu * x + x * x - sigma) ); 
}
