include("library.jl")
using .library
using Test
using Random

function main()
    ############## INITIALIZE THE VARIOUS PARAMETERS #################################################################
    num_particles = 1                          # Number of particles
    num_dims = 2                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    mc_burn_in = 0.2                           # Fraction of steps before sampling
    num_mc_cycles = 500000                   # Number of steps in the MC algorithm
    num_optimization_steps = 100                # Number of optimization steps

    brute_force_step_length = 0.5              # Step-length in the Brute-force Metropolis
    importance_sampling_step_length = 0.005     # Time-step in the importance sampling algorithm
    learning_rate = 1.0                        # Learning rate in the optimization algorithm 
    D = 0.5                                    # Diffusion constant for importance sampling
    ##################################################################################################################

    ############## SET UP THE SYSTEM #################################################################################
    # system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    system_gibbs = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared_gibbs, interacting)
    ##################################################################################################################

    ############## CALCULATE THE ENERGIES FOR EACH OPTIMIZATION STEP #################################################
    start = time()
    print("Test")
    # @time runOptimizationBruteForce(system, num_optimization_steps, num_mc_cycles, mc_burn_in, brute_force_step_length, learning_rate)
    # @time runOptimizationImportanceSampling(system, num_optimization_steps, 100000, mc_burn_in, importance_sampling_step_length, D, learning_rate)
    @time runOptimizationGibbsSampling(system_gibbs, num_optimization_steps, num_mc_cycles, mc_burn_in, learning_rate)
    runtime = time() - start
    ##################################################################################################################


    ############## RUN A FINAL SIMULATION WITH THE OPTIMAL WEIGHTS FROM ABOVE ########################################
    # runMetorpolisBruteForce(system, num_mc_cycles, mc_burn_in, brute_force_step_length, true)
    runMetropolisImportanceSampling(system, num_mc_cycles, mc_burn_in, importance_sampling_step_length, D, true)
    ##################################################################################################################

    ############## RUN A GRID SEARCH FOR GIVEN METHOD ################################################################
    # write_grid_search_to_files(2, 2, "gs", true)
    ##################################################################################################################


    ############## WRITE ONE SAMPLING TO FILE ######################################################################
    # write_to_file("gs")
    ##################################################################################################################
end

main()
