# Variational Monte Carlo (NOT DONE WITH THE PROJECT YET...)
(The structure of the program is based on the example structure given by the teaching staff in FYS4411. See: https://github.com/Schoyen/variational-monte-carlo-fys4411.)

In this project we study a Bose gas trapped in a potential trap. We start by looking at a non-interacting system, where the there is no forces between the particles in the system. After that we will introduce a new wavefunction where the interactions between the particles are considered.   

The main part of the repository is system.cpp file. This is the class that runs all the Metropolis algorithms and the gradient descent algorithm. It is also responsible for connecting the wavefunctions with the hamiltonians. All the simulation-methods write to file, and the data-analysis and plotting is done in the jupyter-notebooks. The name of the jupyter notebooks describes which system it analyses. 
## Compilling and running the simulation
To run the desired simulation, you can set up the system in main.cpp. When the desired system is initialized you can compile the code using
`./compile_project` and then run the simluations using `./vmc`.  

## Testing
To test the code, compile the code with (while in the testing folder) 
`g++  TEST.cpp  ../sampler.cpp ../system.cpp ../Hamiltonians/*.cpp ../InitialStates/*.cpp ../Math/*.cpp ../Wavefunctions/*.cpp ../particle.cpp  -o test` and then run the tests with `./test "[TEST]"`. If no arguments are given, all tests are run. To run specific tests, change "[TEST]" with one of the following:
* Test numerical simple Gaussian energy calculation with Metropolis Importance Sampling: "[Numerical Simple Gaussian Energy Calculation Metropolis Importance Sampling]".
* Test numerical simple Gaussian energy calculation with Metropolis brute-force: "[Numerical Simple Gaussian Energy Calculation Metropolis Brute Force]".
* Test simple gaussian energy calculation with Metropolis Importance Sampling: "[Simple Gaussian Energy Calculation Metropolis Importance Sampling]".
* Test simple gaussian energy calculation with Metropolis brute-force: "[Simple Gaussian Energy Calculation Metropolis Brute-Force]"
* Test Gradient Descent on Simple Gaussian: "[Gradient Descent]". 

## Benchmarks
In the folders /Data and /Plots various results are presented. In the /Data folder you will find local energies from various runs of the simluations. The files are named according to the simulated systems. In the /Plots folder you will find plots from some of the simluations, also named according to the system simulated. 


 
