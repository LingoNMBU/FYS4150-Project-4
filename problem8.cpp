/*
See description in README.txt
Build: 
g++ -O3 main_omp_outer_loop.cpp -fopenmp -o main_omp_outer_loop.exe
Run: 
./main_omp_outer_loop.exe <A_min> <A_max> <n_A> <output_file_name>


 g++ -O3 main_omp_outer_loop.cpp src/lattice.cpp -I include -fopenmp -o main_omp_outer_loop.exe -larmadillo

 
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
  const double T_min = 2.18; //min temperature
  const double T_max = 2.28; //max temperature
  const int n_T = 20; // number of steps
  const int burnin = 20000;
  const int N = 100;
  const double J = 1;
  const bool ordered = false;
  const int n_iter = 10e5;

  std::cout << "N" << std::endl;
  std::cout << N << std::endl;
  std::cout << std::endl;

  std::cout << "N_iter" << std::endl;
  std::cout << n_iter << std::endl;
  std::cout << std::endl;
  
  // 
  // Outer loop over A values
  // 

  const double delta_T = (T_max - T_min) / (n_T - 1);  // n_A points correspond to (n_A - 1) intervals

  arma::mat TME = arma::mat(n_T,5);

  #pragma omp parallel // Start parallel region
  {

    // Here we start the parallelized loop over A
    #pragma omp for
    for (int i = 0; i < n_T; ++i)
    {
        double T = T_min + i * delta_T;

        Lattice lattice = Lattice(N,J,ordered,T);
        lattice.markov_mc(n_iter);

        arma::vec E = lattice.ME.col(1);
        arma::vec M = lattice.ME.col(0);

        E.shed_rows(0,burnin);
        M.shed_rows(0,burnin);

        TME(i,0) = T;
        TME(i,1) = arma::mean(M);
        TME(i,2) = arma::mean(E);
        TME(i,3) = arma::mean(M%M);
        TME(i,4) = arma::mean(E%E);

    } // End parallelized loop over 

  } // End entire parallel region

  TME.save("TME_N100_20T", arma::csv_ascii);

  // And all was well
  return 0;
}



