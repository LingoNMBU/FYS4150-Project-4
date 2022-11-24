/*
g++ -O3 problem8.cpp src/lattice.cpp -I include -fopenmp -o problem8.exe -larmadillo 
*/

#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include "lattice.hpp"

using namespace std;

int main()
{
  // Read command line arguments
  const double T_min = 2.20; //min temperature
  const double T_max = 2.35; //max temperature
  const int n_T = 28; // number of steps
  const int burnin = 20000;
  const int N = 80;
  const double J = 1;
  const bool ordered = false;
  const int n_iter = 2*10e5;
  const string filename = "TME_N80";
  const double delta_T = (T_max - T_min) / (n_T - 1);  // nT points
  const int N_threads = 7;

  std::cout << "N" << std::endl;
  std::cout << N << std::endl;
  std::cout << std::endl;

  std::cout << "N_iter" << std::endl;
  std::cout << n_iter << std::endl;
  std::cout << std::endl;
  
  arma::mat TME = arma::mat(n_T,5).fill(0.); //Initialize matrix outside of parallel block

  #pragma omp parallel num_threads(N_threads)// Start parallel region
  {//Start prallel block
    // Here we start the parallelized loop over T
    #pragma omp for
    for (int i = 0; i < n_T; ++i)
    {
        arma::vec E,M,E2,M2;
        double T;
        
        T = T_min + i * delta_T;

        Lattice lattice = Lattice(N,J,ordered,T);
        lattice.markov_mc(n_iter);

        E = lattice.ME.col(1);
        M = arma::abs(lattice.ME.col(0));
        M2 = arma::pow(M,2);
        E2 = arma::pow(E,2);

        // if burnin
        //E.shed_rows(0,burnin);
        //M.shed_rows(0,burnin);

        #pragma omp critical  //For safety
        {
        TME(i,0) = T;
        TME(i,1) = arma::mean(M);
        TME(i,2) = arma::mean(E);
        TME(i,3) = arma::mean(M2);
        TME(i,4) = arma::mean(E2);
        }

    } // End parallelized loop 
  } // End parallel\ block

  //save file
  TME.save(filename, arma::csv_ascii);
  return 0;
}



