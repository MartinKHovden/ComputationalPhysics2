# Using Machine Learning (ML) and Markov Chain Monte Carlo (MCMC) for calculating the ground state energy of electrons confined to move in a harmonic oscillator trap.
GitHub repo for second project in FYS4411: Computational Physics II at the University of Oslo. In this project, ML and MCMC are used for calculating the ground state energy of electrons confined to move in a harmonic oscillator trap. Both the interacting and the non-interacting system is studied. The methods are inspired by the work of G. Carleo and M. Troyer in their article __"Solving the quantum many-body problem with artificial neural networks"__, __Science 355, Issue 6325, pp. 602-606 (2017)__, found at: https://science.sciencemag.org/content/355/6325/602


# Structure of the repo
The repo contains the main library.jl file where all the code for running the simulations are found. Plots and benchmarks are found in the folders "/plots" and "/output". The plots are named after what part of the exercise it solves, and all plots are titled with some more info. In addition, there are jupyter notebooks containing the analysis of the data.  

To run the main library file, a Julia installation is needed. See more at: https://julialang.org/. When Julia works properly, the files can be executed:
```
julia filename.jl
```

# Testing of the code
Unit-tests are found in tests.jl. To run the tests, use 
```
julia tests.jl
```
The test-library contains test of the functions for calculating values. The tests of the sampling algorithms are visual tests in the report and jupyter notebooks that aim to see if they produce good chain of samples and if they converge to the correct values. 
