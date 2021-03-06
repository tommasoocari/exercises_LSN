#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "TSP.h"
#include "../Parallel_Random_Number_Generator/random.h"
#include <algorithm>
#include <iomanip>
#include <numeric>

using namespace std;

int main(){
  Input();
  Coordinates();

  path = RandomPath(n_cities);

  Normalize();
  Measure(0);

  istep = 0;
  
  for(int iblk = 1; iblk < n_blocks + 1; iblk++)
    {
      Reset(iblk);
      
      for(int i = 0; i < n_steps_block; i++)
	{
	  Step(beta);
	  istep++;
	  Normalize();
	  if(istep % 100 == 0)
	    Measure(istep);
	}

      End_Block(iblk,beta,n_steps_block);
      
    }

  Write_Best_Path();

}

void Input(void)
{
  cout << endl << "Travel Salesman Problem" << endl;
  cout << endl << "Simulated Annealing" << endl;
  
  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  //Read input for information
  ifstream ReadInput;
  ReadInput.open("input.dat");

  cout << endl << "Reading input from file input.dat" << endl;

  cout << endl << "Number of cities: " << n_cities << endl;

  ReadInput.close();

  ReadInput.open("schedule.dat");

  cout << endl << "Reading annealing schedule from file schedule.dat" << endl;

  string appo;
  double appo1, appo2;
  n_steps = 0;

  ReadInput >> appo >> appo;
  
  while(!ReadInput.eof())
    {
      ReadInput >> appo1 >> appo2;
      annealing_schedule.push_back({appo1,appo2});
      n_steps = n_steps + appo1; // update #steps
    }
  annealing_schedule.pop_back();
  n_steps = n_steps - 1; 
  
 n_blocks = annealing_schedule.size();
  
 cout << endl << "Number of algorithm steps: " << n_steps << endl;
 cout << "Number of blocks: " << n_blocks << endl;

 cout << endl << "Annealing schedule" << endl;
 cout << setw(24) << "Number of steps per block" << setw(24) << "beta" << endl;
 for (auto i : annealing_schedule)
   {
     for(auto j : i) cout << setw(24) << j;
     cout << endl;
   }
 cout << endl;
}

void Coordinates(void)
{
  //Read the coordinates of the cities
  ifstream ReadCoordinates;
  ReadCoordinates.open("input.coordinates");

  cout << endl << "Reading coordinates of the cities from the file input.coordinates" << endl;
  
  for(int i = 0; i < n_cities; i++)
    {
      double x, y;
      ReadCoordinates >> x;
      ReadCoordinates >> y;
      coordinates.push_back({x,y});
    }
  
  ReadCoordinates.close();
}

vector<double> RandomPath(int n_cities)
{
  vector<double> indiv;
  indiv.push_back(0.0);

  indiv.push_back(1 + rand() % (n_cities-1)); //random number beetween 1 and n_genes - 1
  
  while(indiv.size() != n_cities)
    {
      int random = 1 + rand() % (n_cities-1);
      if(find(indiv.begin(), indiv.end(), random) == indiv.end())
	{
	  indiv.push_back(random);
	}
    }
  indiv.push_back(0.0);
  return indiv;
}

double Fitness(vector<double> indiv)
{
  double sum = 0.0;
  
  for(int i = 0; i < n_cities; i++)
    {
      sum = sum + Distance(coordinates[(int)indiv[i]], coordinates[(int)indiv[(i+1)%n_cities]]);
    }

  return sum;
}

double Distance(vector<double> v1, vector<double> v2)
{
  if(v1.size() != v2.size()) return 0;
  double su2=0;
  for(int i=0 ; i<v1.size(); i++) su2 = su2 + pow(v1[i] - v2[i],2);
  return sqrt(su2);
}

void Normalize(void)
{
  if(path[1] > path[n_cities - 1])
      reverse(path.begin()+1, path.end()-1);
  path[n_cities] = Fitness(path);
}

void Reset(int iblk)
{
  accepted = 0.0;
  attempted = 0.0;
  
  beta =annealing_schedule[iblk-1][1];
  n_steps_block = annealing_schedule[iblk-1][0];
}

void Measure(int istep)
{
  double fitness;
  ofstream Fit;

  Fit.open("output.fitness.0",ios::app);
  
  fitness = path[n_cities];

  Fit << istep << " " << fitness << endl;

  Fit.close();
}

void Step(double beta)
{
  vector<double> old_path = path;
  vector<double> new_path = path;
  double p;
  int r,r1,r2;
  int gene1, gene2;
  int length_max, length;

  r = rnd.Rannyu(0,4); //mutations

  if(r == 0) //swap
    {
      r1 = rnd.Rannyu(1, n_cities- 1);
      r2 = rnd.Rannyu(1, n_cities- 1);
      gene1 = min(r1,r2);
      gene2 = max(r1,r2);

      swap(new_path[gene1], new_path[gene2]);
    }
  if(r == 1) //reverse
    {
      r1 = rnd.Rannyu(1, n_cities);
      r2 = rnd.Rannyu(1, n_cities);
      gene1 = min(r1,r2);
      gene2 = max(r1,r2);

      reverse(new_path.begin() + gene1, new_path.begin() + gene2);
    }
  if(r == 2) //blockswap
    {
      r1 = rnd.Rannyu(1, n_cities- 1);
      r2 = rnd.Rannyu(1, n_cities- 1);
      gene1 = min(r1,r2);
      gene2 = max(r1,r2);
      length_max = n_cities-gene2 -1;
      length = rnd.Rannyu(1,length_max + 1);

      for(int i = 0; i < length + 1; i++)
	swap(new_path[gene1+i], new_path[gene2+i]);
    }
  else //blockshift
    {
      r1 = rnd.Rannyu(1, n_cities- 1);
      r2 = rnd.Rannyu(1, n_cities- 1);
      gene1 = min(r1,r2);
      gene2 = max(r1,r2);
      length_max = n_cities-gene2;
      length = rnd.Rannyu(1,length_max + 1);

      rotate(new_path.begin()+gene1, new_path.begin()+gene2, new_path.begin()+gene2+length);
      
    }
  
  p = exp (- beta * (Fitness(new_path) - Fitness(old_path)));
  if(p > 1) p = 1;
  attempted = attempted + 1;
  
  if(rnd.Rannyu() < p) //metropolis acceptation
    {
      path = new_path;
      accepted = accepted + 1;
    }
  
}

void End_Block(int iblk, double beta, double n_steps_block)
{
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
  cout << "beta = " << beta << endl;
  cout << "Number of steps = " << n_steps_block << endl;
  cout << "----------------------------" << endl << endl;
}

void Write_Best_Path(void)
{
  vector<double> best_path;
  ofstream Path;

  Path.open("output.path.0");
  best_path = path;

  for(int i = 0; i < n_cities; i++)
    Path << coordinates[(int)best_path[i]][0] << " "
	 << coordinates[(int)best_path[i]][1] << endl;

  Path << coordinates[(int)best_path[0]][0] << " "
	 << coordinates[(int)best_path[0]][1] << endl;

  Path.close();
  
  rnd.SaveSeed();
}

