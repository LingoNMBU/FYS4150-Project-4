
#include "lattice.hpp"

int main()
{
    Lattice lattice4 = Lattice(2,1, false);
    lattice4.markov_mc(10e3);
    lattice4.MECX.save("MECX_4.csv", arma::csv_ascii);




}