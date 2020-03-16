# Variational Monte Carlo
(The structure of the program is based on the example structure given by the teaching staff in FYS4411. See: https://github.com/Schoyen/variational-monte-carlo-fys4411.)

In this project we study particles trapped in different potentials. We start by looking at a non-interacting system, where the there is no forces between the particles in the system. After that we will introduce a new wavefunction where the interactions between the particles are considered. 

## Compilling and running the project
To run the desired simulation, you can set up the system in main.cpp. When the desired system is initialized you can compile the code using
`./compile_project` and then run the simluations using `./vmc`.  

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

#### Windows
Compilation of the project using Windows is still an open question to me, but please include a pull-request if you've got an example. CMake should be OS-independent, but `make` does not work on Windows.

