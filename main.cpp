
#include "lattice.hpp"

int main()
{

    bool prob4 = true;
    if (prob4)
    {
    Lattice lattice4 = Lattice(2,1, false, 1.0);
    lattice4.markov_mc(10e5);
    lattice4.MECX.save("MECX_4.csv", arma::csv_ascii);
    }

    bool prob5 = true;
    if (prob5)
    {
    Lattice lattice51 = Lattice(20, 1, false, 1.0);
    lattice51.markov_mc(10e4);
    lattice51.MECX.save("MECX_5_10_rand20.csv", arma::csv_ascii);

    Lattice lattice52 = Lattice(20, 1, false, 2.4);
    lattice52.markov_mc(10e4);
    lattice52.MECX.save("MECX_5_24_rand20.csv", arma::csv_ascii);

    Lattice lattice53 = Lattice(20, 1, true, 1.0);
    lattice53.markov_mc(10e4);
    lattice53.MECX.save("MECX_5_10_ord20.csv", arma::csv_ascii);

    Lattice lattice54 = Lattice(20, 1, true, 2.4);
    lattice54.markov_mc(10e4);
    lattice54.MECX.save("MECX_5_24_ord20.csv", arma::csv_ascii);
    }





}