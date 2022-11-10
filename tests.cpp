
#include "lattice.hpp"
#include <cassert>

int main()
{
    std::cout.precision(17);

    //Energy test
    Lattice lattice0 = Lattice(2,1, true);
    double E0 = lattice0.total_energy(lattice0.spins)(0);
    assert(E0 == -8);

    Lattice lattice1 = Lattice(3,1, true);
    double E1 = lattice1.total_energy(lattice1.spins)(0);
    // -1 per neighbour pair, N*N*2 neihgbours giving -18 total
    assert(E1 == -18);

    lattice1.change_spin(2,2);
    lattice1.change_spin(1,1);
    double E12 = lattice1.total_energy(lattice1.spins)(0);
    //Changing two non-neighbour spins, change four neighbours each, which changes the sum of spins +8 each, giving +16 in total
    assert(E12 == -2); 

    Lattice lattice2 = Lattice(10,1, true);  
    lattice2.markov_mc(1000);

    //std::cout << std::endl;
    //std::cout << lattice2.spins << std::endl;
    //std::cout << std::endl;   

    //lattice2.markov_mc(10);

    //std::cout << std::endl;
    //std::cout << lattice2.spins << std::endl;
    //std::cout << std::endl;   

    //lattice2.markov_mc(10);

    std::cout << std::endl;
    std::cout << lattice2.spins << std::endl;
    std::cout << std::endl;   




    return 0;
}