#ifndef __lattice_hpp__  
#define __lattice_hpp__


#include <armadillo>
#include <random>

class Lattice
{
  private:

  public:
    int N, N_spins;
    double J, T, e, m, kb;
    bool ordered;
    arma::imat spins;
    std::map<int, double> deltaE;
    arma::mat MECX;


    //constructor
    Lattice(int N, double J, bool ordered=false, double T=1.0);

    void change_spin(int i, int j);

    arma::vec total_energy(arma::imat spins_m);

    arma::vec total_magnetization(arma::imat spins_m);

    void markov_mc(int n_iter);

    double sample_boltzmann(double T, arma::imat spin_m, double Z);

    double normalization_constant(double T);

};
#endif 