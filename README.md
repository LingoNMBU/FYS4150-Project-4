# FYS4150-Project-4
Project 4 computational physics: Jupyter is ok i guess

\n

problem 4,5 and 6: data is generated with:

g++ -O3 problem4_5_6.cpp src/lattice.cpp -I include -o problem4_5_6.cpp.exe -larmadillo 

&

./problem4_5_6.exe

\n

problem 8: data is generated with the following:

g++ -O3 problem8.cpp src/lattice.cpp -I include -fopenmp -o problem8.exe -larmadillo 

&

./problem8.exe 

different lattice sizes are modified in the cpp file

\n

problem 9 is completed in the Jupyter Notebook problem9.ipynb

Plotting is done in the relevant jupyter notebooks, filepaths for imports may need modification
