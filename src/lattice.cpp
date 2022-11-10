#ifndef __lattice_cpp__  
#define __lattice_cpp__

#include "lattice.hpp"

// Constructor
Lattice::Lattice(int N_in, double J_in, bool ordered_in, double T_in)
{
  N = N_in;               // Lattice size
  J = J_in;             //COupling constant
  ordered = ordered_in  ;   // if true makes lattice of only 1s 
  T = T_in;             //Temperature
  N_spins = N*N;
  kb = 1.380649e-23;

  if (!ordered)
  {
    spins = arma::randi(N, N, arma::distr_param(0, 1));
    spins.elem(find(spins == 0)).fill(-1);
  }
  else
  {
    spins = arma::imat(N,N).fill(1);
  }

  //Storing values of exp(-kb * deltaE) for different deltaE
  deltaE[-8] = exp(8./T);          // -8:1
  deltaE[-4] = exp(4./T);          // -4:1
  deltaE[0] = 1.;           // 0:1
  deltaE[4] = exp(-4./T);    // 4:exp(-4/T)
  deltaE[8] = exp(-8./T);   // 8: exp(-8/T)
}

// Change the spin of a single particle in the lattice
void Lattice::change_spin(int i, int j)
{
  spins(i,j) = spins(i,j) * (-1);
}


// Calculate the total energy of the system
arma::vec Lattice::total_energy(arma::imat spins_m)
{
    //spins
    // {-1, -1, +1}
    // {-1, +1, -1}
    // {+1, -1, +1}

    //indices
    // {00, 01, 02}
    // {10, 11, 12}
    // {20, 21, 22}

    //bonds
    //   |    |    |  
    // - 00 - 01 - 02
    //   |    |    |
    // - 10 - 11 - 12
    //   |    |    |
    // - 20 - 21 - 22

    // can iterate two backwards neighbours per particle to cover whole lattice without double counting


    //Applying periodic boundary condition by adding the last rows and columns to the start of the matrix
    arma::imat spins_o;
    spins_o = spins_m;

    spins_o.insert_rows(0,spins_o.row(N-1));
    spins_o.insert_cols(0,spins_o.col(N-1));

    spins_o(0,0) = 0; //Zero for clarity, wont be used

    arma::vec Es = arma::vec(2);

  double E_tot =  0;
  double E2_tot = 0;
  for(int i = 1; i < N+1; i++)
  {
    for (int j = 1; j < N+1; j++)
    {
        //neighbour one row up
        E_tot += spins_o(i-1, j)*spins_o(i, j);
        E2_tot += E_tot*E_tot;
        //neighbour one column back
        E_tot += spins_o(i, j-1)*spins_o(i, j);
        E2_tot += E_tot*E_tot;
    }

  }
  Es(0) = -J*E_tot;
  Es(1) = -J*E2_tot;

  return Es;
}
arma::vec Lattice::total_magnetization(arma::imat spins_m)
{
    double M_tot = 0;
    double M_2_tot = 0;
    arma::vec Ms = arma::vec(2);
    for (int i = 0; i < N; i++)
     {
        for (int j = 0; j < N; j++)
        {
            M_tot += spins(i,j);
            M_2_tot += M_tot*M_tot;
        }        
     }
     Ms(0) = abs(M_tot);
     Ms(1) = M_2_tot;
    return Ms;
}

//double Lattice::sample_boltzmann(double T, arma::mat spin_m, double Z)
//{
//    double E_s, beta, kb;
//
//    kb = 1.380649 * pow(10.,-23.);
//
//    beta = 1./(kb*T);
//
//    E_s = total_energy(spin_m);
//
//    p_s = 1./Z * exp(-beta*E_s);
//
//    return p_s;
//}


//double Lattice::normalization_constant(double T)
//{
    //All possible states

    // 1  1  1  1

    // 1  1  1 -1
    // 1  1 -1  1
    // 1 -1  1  1
    //-1  1  1  1

    // 1  1 -1 -1
    // 1 -1 -1  1
    //-1 -1  1  1    
    //-1  1  1 -1
    //-1  1 -1  1
    // 1 -1  1 -1

    // 1 -1 -1 -1
    //-1  1 -1 -1
    //-1 -1  1 -1
    //-1 -1 -1  1

    //-1 -1 -1 -1


    // 1  1  1  1 , flip 0
    // 1  1  1 -1 , flip 1
    // 1  1 -1 -1 , flip 0
    // 1  1 -1  1 , flip 1

    // 1 -1 -1  1 , flip 2
    // 1 -1 -1 -1 , flip 1
    // 1 -1  1 -1 , flip 0
    // 1 -1  1  1 , flip 3

    //-1 -1  1  1, flip 0
    //-1 -1  1 -1, flip 1
    //-1 -1 -1  1, flip 0
    //-1 -1 -1 -1, flip 2

    //-1  1 -1 -1, flip 0
    //-1  1 -1  1, flip 1
    //-1  1  1  1, flip 0
    //-1  1  1 -1, 

    

    //-1 -1 -1 -1
    //-1 -1 -1  1
    //-1 -1  1  1
    //-1 -1  1 -1


    
    //-1  1 -1 -1
    //-1  1 -1  1
    //-1  1  1  1
    //-1  1  1 -1


//}



void Lattice::markov_mc(int n_iter)
{
    //We need a probability distribution from which we sample our states

    //We need a method to propose a new state for the system
    //we need a method for accepting a new state

    //Metropolis Hasting acceptance rule
    //P(x -> x') = T(xi -> x')A(xi -> x')
    //1. Generate x' according to T(xi -> x')
    //2. Compute acceptance probability as A = min(1, p(x')/p(xi) * T(x' ->  xi)/ T(xi -> x'))

    //   If T(x' -> xi) = T(xi -> x'), simplifies to Metropolis rule

    //Hint: When computing P(x')/p(xi), computing the normalization constant can be quite expensive, and so calculating the ratio is not necessary. Things cancel


    //We need to calculate different things for the state and store them

    //We need to update our state

  //initial values
  MECX = arma::mat(n_iter+1, 4).fill(0.);
  MECX(0,0) = total_magnetization(spins)(0);
  MECX(0,1) = total_energy(spins)(0);

  int i = 0;
  while(i < n_iter)
  {
      i++;
    //Suggested algorithm

    arma::imat state0 = spins;

    //1. Generate candidate state by flipping one random spin

    int rand_i = arma::randi(arma::distr_param(0, N-1));
    int rand_j = arma::randi(arma::distr_param(0, N-1));

    change_spin(rand_i,rand_j);

    arma::imat state1 = spins;

    //2. Calculate the ratio p(s')/p(si), do this efficiently!(stuff cancels)

    int E0 = total_energy(state0)(0);
    int E1 = total_energy(state1)(0);

    double p1_p0 = deltaE[(E1 - E0)];
    double A = std::min(1.0, p1_p0);

    //p(s')/p(si) = (exp(-kb*E(s'))/Z) / (exp(-kb*E(si))/Z) = exp(-kb*E(s') - -kb*E(si) ) = exp(-kb*(E(s') - E(si)) ) 

    //3. Generate r form uniform distribution, accept if A > r, reject if A < r
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double r = distribution(generator);

    if (A < r)
    {
      //New state rejected
      spins = state0;

      double M0 = total_magnetization(spins)(0)/N_spins;
      //C = (1./(kb*T*T)) * (E^2)

      MECX(i,0) = M0;
      MECX(i,1) = E0/N_spins;


    }
    else
    {
      //New state accepted
      double M1 = total_magnetization(spins)(0)/N_spins;
      //C = (1./(kb*T*T)) * (E^2)

      MECX(i,0) = M1;
      MECX(i,1) = E1/N_spins;

      std::cout <<std::endl;
      std::cout << total_energy(spins)(0) << std::endl;
      std::cout << std::endl;   
    }

     
  }
    //4. Repeat 1-> N times, is one MC cycle.
    //5. Calculate quantities
}






#endif 