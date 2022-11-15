#ifndef __lattice_hpp__  
#define __lattice_hpp__


#include <armadillo>
#include <random>
#include <chrono>
#include <algorithm>

class Lattice
{
  private:

  public:
    int N, N_spins;
    double J, T, e, m, kb, beta;
    bool ordered;
    arma::imat spins;
    std::map<int, double> deltaE;
    arma::mat MECX;


    //constructor
    Lattice(int N, double J, bool ordered=false, double T=1.0);

    void change_spin(int i, int j);

    int total_energy(arma::imat spins_m);

    int energy_single(int i, int j);

    int total_magnetization(arma::imat spins_m);

    void markov_mc(int n_iter);

    void markov_mc2(int n_iter);

    double sample_boltzmann(double T, arma::imat spin_m, double Z);

    double normalization_constant(double T);

};
#endif  