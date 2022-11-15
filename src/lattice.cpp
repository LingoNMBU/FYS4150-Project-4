#ifndef __lattice_cpp__  
#define __lattice_cpp__

#include "lattice.hpp"

// Constructor
Lattice::Lattice(int N_in, double J_in, bool ordered_in, double T_in)
{
  
  kb = 1.380649e-23;
  N = N_in;                 // Lattice size
  J = J_in;                 //COupling constant
  ordered = ordered_in;   // if true makes lattice of only 1s 
  T = (J*T_in)/kb;               //Temperature
  N_spins = N_in*N_in;
  beta = 1/(kb*T);


  if (!ordered)
  {
    arma::arma_rng::set_seed_random();
    spins = arma::randi(N, N, arma::distr_param(0, 1));
    spins.elem(find(spins == 0)).fill(-1);
  }
  else
  {
    spins = arma::imat(N,N).fill(1);
  }

  //Storing values of exp(-beta * deltaE) for different deltaE
  deltaE[-8] = 1.0;          // -8:1 , in reality ~2980, but alwways higher than 1
  deltaE[-4] = 1.0;          // -4:1 , in reality ~54, but always higher than 1
  deltaE[0] = 1.0;                         // 0:1
  deltaE[4] = exp(-(4.)/(beta));            // 4:exp(-4/T)
  deltaE[8] = exp(-(8.)/(beta));            // 8: exp(-8/T)
}

// Change the spin of a single particle in the lattice
void Lattice::change_spin(int i, int j)
{
  spins(i,j) = -1*spins(i,j);
}


// Calculate the total energy of the system
int Lattice::total_energy(arma::imat spins_m)
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

    spins_o(0,0) = 0; //Zero for clarity, wont be used, corner.

    int E = 0;

  int E_tot =  0;
  //double E2_tot = 0;
  for(int i = 1; i < N+1; i++)
  {
    for (int j = 1; j < N+1; j++)
    {
        //neighbour one row up and neighbour one column back
        E_tot +=  spins_o(i-1, j)*spins_o(i, j) + spins_o(i, j-1)*spins_o(i, j);;
    }

  }
  E = -J*E_tot;
  //Es(1) = abs(-J*E2_tot);

  return E;
}

int Lattice::energy_single(int i, int j)
{
  int E_single = 0;
  int s_mid   = spins(i, j);
  int s_left  = spins(i, (j - 1 + N) % N);  // if j = 0 -> L-1%L = L-1
  int s_right = spins(i, (j + 1) % N);      // if j = L -1 -> L%L = 0
  int s_up    = spins((i - 1 + N) % N, j);  // if i = 0 -> -1 + L%L = L-1
  int s_down  = spins((i + 1) % N, j);      // if i = L-1 -> L%L = 0

  //std::cout << "spin left" << std::endl;
  //std::cout << s_left << std::endl;
  //std::cout << std::endl;
  //std::cout << "spin right" << std::endl;
  //std::cout << s_right << std::endl;
  //std::cout << std::endl;
  //std::cout << "spin up" << std::endl;
  //std::cout << s_up << std::endl;
  //std::cout << std::endl;
  //std::cout << "spin down" << std::endl;
  //std::cout << s_down << std::endl;
  //std::cout << std::endl;
  //std::cout << "spin mid" << std::endl;
  //std::cout << s_mid << std::endl;
  //std::cout << std::endl;

  E_single = -J*s_mid*(s_left + s_right + s_up + s_down); // E = -J*sum over bonds(E_s*E_neighbour)

  //std::cout << "E_single" << std::endl;
  //std::cout << E_single << std::endl;
  //std::cout << std::endl;

  return E_single;
}

int Lattice::total_magnetization(arma::imat spins_m)
{
    int M_tot = 0;
    for (int i = 0; i < N; i++)
     {
        for (int j = 0; j < N; j++)
        {
            M_tot += spins(i,j);
            //M_2_tot += M_tot*M_tot;
        }        
     }

    return M_tot;
}

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
  //Random seed

  //initial values
  MECX = arma::mat(n_iter+1, 2).fill(0.);
  int E0, E1, M0, M1;

  E0 = total_energy(spins);      
  M0 = total_magnetization(spins);

  MECX(0,0) = M0/N_spins;
  MECX(0,1) = E0/N_spins;

  //Generate N_spins random indexes
  arma::arma_rng::set_seed_random();
  arma::imat rand_is = arma::randi(n_iter, N_spins, arma::distr_param(0, N-1));
  arma::imat rand_js = arma::randi(n_iter, N_spins, arma::distr_param(0, N-1));

  //Generate N_pins random floats between 0 and 1(
  arma::mat rs(n_iter, N_spins, arma::fill::randu);
  int iter = 0;
  for (int i = 0; i < n_iter-1; i++)
  {
    //Algorithm
    for (int j = 0; j < N_spins; j++)
    {
      //1. Generate candidate state by flipping one random spin
      change_spin(rand_is(i,j),rand_js(i,j));

      //calculate energy
      E1 =  total_energy(spins);

      //2. Calculate the ratio p(s')/p(si)
      M1 = total_magnetization(spins);

      int dE = (E1 - E0);  

      double p1_p0 = deltaE[dE];
      double A = std::min(1.0, p1_p0);

      //3. Generate r form uniform distribution
      double r = rs(i,j);


      //4. Accept if A > r, reject if A < r
      if (A < r)
      {
        //New state rejected
        change_spin(rand_is(i,j),rand_js(i,j));
      }
      else
      {
        //New state accepted
        E0 += dE;
        M0 += (M1 - M0);
      }
    }
    MECX(i+1,0) = M0/N_spins;
    MECX(i+1,1) = E0/N_spins;

    iter++;
    if (iter%20000 == 0)
    {
      std::cout << std::endl;
      std::cout << iter << std::endl;
      std::cout << std::endl;
    }
  }
}


void Lattice::markov_mc2(int n_iter)
{

  //initial values
  std::cout << std::endl;
  std::cout << T << std::endl;
  std::cout << std::endl;

  MECX = arma::mat(n_iter, 2).fill(0.);
  arma::ivec Ms0 = arma::ivec(2).fill(0);
  arma::ivec Ms1 = arma::ivec(2).fill(0);

  int E = total_energy(spins);      
  int M = total_magnetization(spins);

  MECX(0,0) = M/N_spins;
  MECX(0,1) = M/N_spins;

  //Generate n_iter x N_spins random indexes
  arma::arma_rng::set_seed_random();
  arma::imat rand_is = arma::randi(n_iter, N_spins, arma::distr_param(0, N-1));
  arma::imat rand_js = arma::randi(n_iter, N_spins, arma::distr_param(0, N-1));

  //Generate n_iter x N_spins random floats between 0 and 1
  arma::mat rs(n_iter, N_spins, arma::fill::randu);

  //Test
  int rejected = 0;
  int accepted = 0;

  //Start MC cycle
  int iter = 0;
  for (int i = 0; i < n_iter-1; i++)
  {  
    //Inner loop
    for (int j = 0; j < N_spins; j++)
    {
      //1. Generate candidate state by flipping one random spin
      int i_ind = rand_is(i,j);
      int j_ind = rand_js(i,j);

      int M_s0 = spins(i_ind, j_ind);
      int E_s0 = energy_single(i_ind, j_ind );  

      //std::cout << "spin0" << std::endl;
      //std::cout << M_s0 << std::endl;
      //std::cout << std::endl;  

      //std::cout << "Es0" << std::endl;
      //std::cout << E_s0 << std::endl;
      //std::cout << std::endl;  

      change_spin(i_ind, j_ind);
      int M_s1 = spins(i_ind, j_ind);
      int E_s1 = energy_single(i_ind, j_ind );

      //std::cout << "spin1" << std::endl;
      //std::cout << M_s1 << std::endl;
      //std::cout << std::endl;  

      //std::cout << "Es1" << std::endl;
      //std::cout << E_s1 << std::endl;
      //std::cout << std::endl;
      int dE = E_s1 - E_s0; 
      int dM = M_s1 - M_s0;

      //2. Calculate the ratio p(s')/p(si)
      double p1_p0 = deltaE[dE];
      double A = std::min(1.0, p1_p0);

      //std::cout << "dE" << std::endl;
      //std::cout << dE << std::endl;
      //std::cout << std::endl;  

      //std::cout << "p1/p0" << std::endl;
      //std::cout << p1_p0 << std::endl;
      //std::cout << std::endl;  
      //    
      //std::cout << "Acceptance Rate" << std::endl;
      //std::cout << A << std::endl;
      //std::cout << std::endl;  

      //3. Generate r form uniform distribution
      //4. Accept if A > r, reject if A < r

      if (A < rs(i,j))
      {
        //New state rejected
        change_spin(i_ind, j_ind);
        rejected += 1;
      }
      else
      {
        M += dM;
        E += dE;
        //std::cout << std::endl;
        //std::cout << "change" << std::endl;
        //std::cout << std::endl;
        accepted += 1;

      }
    }
    MECX(i,0) = M/N_spins;
    MECX(i,1) = E/N_spins;
    iter++;

    //std::cout << "accepted" << std::endl;
    //std::cout << accepted << std::endl;
    //std::cout << std::endl;  

    //std::cout << "rejected" << std::endl;
    //std::cout << rejected << std::endl;
    //std::cout << std::endl;   

    if (iter%20000 == 0)
    {
      std::cout << std::endl;
      std::cout << iter << std::endl;
      std::cout << std::endl;
    }
  }
}




#endif 