
#include "lattice.hpp"

using namespace std;

int main()
{

    bool prob4 = true;
    if (prob4)
    {
        Lattice lattice4 = Lattice(40, 1, false, 1.0);
        lattice4.markov_mc(10e5);
        lattice4.ME.save("ME_4.csv", arma::csv_ascii);
    }


    bool prob5 = true;
    if (prob5)
    {
        Lattice lattice51 = Lattice(20, 1, false, 1.0);
        lattice51.markov_mc(10e5);
        lattice51.ME.save("ME_5_10_rand20.csv", arma::csv_ascii);

        Lattice lattice52 = Lattice(20, 1, false, 2.4);
        lattice52.markov_mc(10e5);
        lattice52.ME.save("ME_5_24_rand20.csv", arma::csv_ascii);

        Lattice lattice53 = Lattice(20, 1, true, 1.0);
        lattice53.markov_mc(10e5);
        lattice53.ME.save("ME_5_10_ord20.csv", arma::csv_ascii);

        Lattice lattice54 = Lattice(20, 1, true, 2.4);
        lattice54.markov_mc(10e5);
        lattice54.ME.save("ME_5_24_ord20.csv", arma::csv_ascii);
    }






}