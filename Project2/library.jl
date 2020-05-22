module library

# Library with different methods for Markov Chain Monte Carlo sampling and optmization
# for simluating two electrons in a potential trap.
# Code based on the lecture notes written by Morten Hjorth-Jensen for FYS441: Computational Physics II
# and the code provided for the course at:
# https://github.com/CompPhysics/ComputationalPhysics2/blob/gh-pages/doc/Programs/BoltzmannMachines/MLcpp/src/

using LinearAlgebra
using Random
using Profile

"""
    NQS

Struct for collecting the variables and parameters in the Boltzman machine.
"""
struct NQS

    num_particles::Int64
    num_dims::Int64

    b::Array{Float64, 2}
    a::Array{Float64, 2}

    w::Array{Float64, 2}

    x::Array{Float64, 2}
    h::Array{Float64, 2}

    sigma_squared::Float64
    interacting::Bool

end

"""
    computeInteractionTerm()
"""
function computeInteractionTerm(nqs::NQS)

    interaction_term = 0

    for i = 1:nqs.num_dims:nqs.num_particles*nqs.num_dims

        for j = i+nqs.num_dims:nqs.num_dims:nqs.num_particles*nqs.num_dims

            # println(i, " ", j)
            r_ij = 0

            for k = 0:nqs.num_dims-1

                r_ij += (nqs.x[i + k] - nqs.x[j + k])^2

            end

            r_ij = sqrt(r_ij)

            interaction_term += 1.0/r_ij

        end

    end

    return interaction_term



end

"""
    computeLocalEnergy(nqs::NQS)

Computes the local energy of the system with the given
hamiltonian and wavefunction.
"""
function computeLocalEnergy(nqs::NQS, precalc)

    #Calculates the term that will be used multiple times to increase speed.
    # precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    local_energy::Float64 = 0

    for m = 1:num_visible

        ln_psi_derivative = -(1.0/nqs.sigma_squared)*(nqs.x[m] - nqs.a[m])
        ln_psi_double_derivative = -1.0/nqs.sigma_squared

        for n=1:num_hidden

            ln_psi_derivative += (1.0/nqs.sigma_squared)*nqs.w[m,n]/(exp(-precalc[n]) + 1.0)
            ln_psi_double_derivative += (1.0/nqs.sigma_squared^2)*(exp(precalc[n])/((exp(precalc[n])+1)^2))*(nqs.w[m,n]^2)

        end

        local_energy += -(ln_psi_derivative)^2 - ln_psi_double_derivative + nqs.x[m]^2

    end

    local_energy *= 0.5

    if nqs.interacting
        local_energy += computeInteractionTerm(nqs)
    end

    return local_energy

end

function computeLocalEnergyGibbs(nqs::NQS, precalc::Array{Float64, 2})

    #Calculates the term that will be used multiple times to increase speed.
    # precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    local_energy::Float64 = 0

    for m = 1:num_visible

        ln_psi_derivative = -(0.5/nqs.sigma_squared)*(nqs.x[m] - nqs.a[m])
        ln_psi_double_derivative = -0.5/nqs.sigma_squared

        for n=1:num_hidden

            ln_psi_derivative += (0.5/nqs.sigma_squared)*nqs.w[m,n]/(exp(-precalc[n]) + 1.0)
            ln_psi_double_derivative += (0.5/(nqs.sigma_squared)^2)*(exp(precalc[n])/((exp(precalc[n])+1)^2))*(nqs.w[m,n]^2)

        end

        local_energy += -(ln_psi_derivative)^2 - ln_psi_double_derivative + nqs.x[m]^2

    end

    local_energy *= 0.5

    if nqs.interacting
        local_energy += computeInteractionTerm(nqs)
    end

    return local_energy

end

"""
    computePsi(nqs::NQS)

Computes the wavefunction value of the system with the given
hamiltonian.
"""
function computePsi(nqs::NQS, precalc::Array{Float64, 2})

    # precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    #Extracts the number of hidden and visible units.
    # num_visible = size(nqs.x)[1]
    num_visible = nqs.num_particles*nqs.num_dims
    num_hidden = size(nqs.h)[1]

    exp_argument = 0

    for i=1:num_visible
        exp_argument += (nqs.x[i] - nqs.a[i])^2
    end

    exp_argument /= (2*nqs.sigma_squared)

    prod_term = 1.0

    for j=1:num_hidden
        prod_term *= (1.0 + exp(precalc[j]))
    end

    return exp(-exp_argument)*prod_term
end

"""
    computePsiParameterDerivative(nqs::NQS)

Computes the derivatives of psi with respect to the parameters in
the Boltzman machine.
"""
function computePsiParameterDerivative!(nqs::NQS, psi_derivative_a, psi_derivative_b, psi_derivative_w, precalc::Array{Float64, 2})

    num_visible::Int64 = size(nqs.x)[1]
    # num_visible = nqs.num_particles*nqs.num_dims

    num_hidden::Int64 = size(nqs.h)[1]

    #Calculates the term that will be used multiple times to increase speed.
    # precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    # psi_derivative_parameter_a::Array{Float64, 2} = zeros(Float64, num_visible ,1)
    # psi_derivative_parameter_a[:, 1] = (1.0/nqs.sigma_squared)*(nqs.x - nqs.a)
    psi_derivative_parameter_a = (1.0/nqs.sigma_squared)*(nqs.x - nqs.a)

    # println("TEST1")

    # psi_derivative_parameter_b::Array{Float64, 2} = zeros(Float64, num_hidden ,1)
    # psi_derivative_parameter_b[:, 1] = ones(Float64, num_hidden, 1)./(exp.(-precalc) + ones(Float64, num_hidden, 1))
    # psi_derivative_parameter_b = ones(Float64, num_hidden, 1)./(exp.(-precalc) + ones(Float64, num_hidden, 1))

    for n = 1:num_hidden
        psi_derivative_b[n, 1] = 1.0/((exp(-precalc[n]) + 1.0))
    end

    # println("TEST2")

    # print(num_visible)
    # print(nqs.x)
    # print(precalc)
    # psi_derivative_parameter_w::Array{Float64, 2} = zeros(Float64, num_visible, num_hidden)
    # psi_derivative_parameter_w = nqs.x/(nqs.sigma_squared*(exp.(-precalc)+ones(Float64, num_hidden, 1)))
    for n=1:num_hidden

        for m=1:num_visible

            psi_derivative_w[m,n] = nqs.x[m,1]/(nqs.sigma_squared*(exp(-precalc[n]) + 1.0))

        end

    end
    # println("TEST3")

    # return psi_derivative_parameter_a, psi_derivative_parameter_b, psi_derivative_parameter_w

end

"""
    setUpSystemRandomUniform(M::Int64, N::Int64, sig_sq::Float64 = 0.5, inter::Bool = false)

Sets up the system parameteres randomly from a uniform distribution.
"""
function setUpSystemRandomUniform(num_particles::Int64, num_dims::Int64, M::Int64, N::Int64, sig_sq::Float64 = 0.5, inter::Bool = false)

    b = rand(Float64, N, 1) .-0.5
    a = rand(Float64, M, 1) .-0.5

    w = rand(Float64, M, N) .-0.5

    x = rand(Float64, M, 1) .-0.5
    # h = rand(Float64, N, 1)
    h = rand(0:1, N, 1)

    interacting = inter

    return NQS(num_particles, num_dims, b, a, w, x, h, sig_sq, interacting)

end

"""
    setUpSystemRandomNorm(M::Int64, N::Int64, sig_sq::Float64 = 0.5, inter::Bool = false)

Sets up the system parameteres randomly from a uniform distribution.
"""
# function setUpSystemRandomNorm(M::Int64, N::Int64, sig_sq::Float64 = 0.5, inter::Bool = false)
#
#     b = randn(Float64, N, 1)
#     a = randn(Float64, M, 1)
#
#     w = randn(Float64, M, N)
#
#     x = randn(Float64, M, 1)
#     h = randn(Float64, N, 1)
#
#     interacting = inter
#
#     return NQS(b, a, w, x, h, sig_sq, interacting)
#
# end

function optimizationStep(nqs::NQS, grad_a::Array{Float64, 2}, grad_b::Array{Float64, 2}, grad_w::Array{Float64, 2}, learning_rate::Float64)

    nqs.a[:] = nqs.a - learning_rate*grad_a
    nqs.b[:] = nqs.b - learning_rate*grad_b
    nqs.w[:,:] = nqs.w - learning_rate*grad_w

end

function metropolisStepBruteForce(nqs::NQS, step_length::Float64, precalc)

    #Extracts the number of hidden and visible units.
    num_visible::Float64 = size(nqs.x)[1]
    num_hidden::Float64 = size(nqs.h)[1]

    coordinate::Int64 = rand(1:num_visible)

    old_wavefunction_value = computePsi(nqs, precalc)

    #Update the coordinate:
    old_coordinate = nqs.x[coordinate]
    nqs.x[coordinate] += (rand(Float64) - 0.5)*step_length


    precalc = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    new_wavefunction_value = computePsi(nqs, precalc)

    U = rand(Float64)

    if U < (new_wavefunction_value^2)/(old_wavefunction_value^2)
        return
    else
        nqs.x[coordinate] = old_coordinate
        return
    end

end

function runMetorpolisBruteForce(nqs::NQS, num_mc_iterations::Int64, burn_in::Float64, step_length::Float64)

    local_energy_sum::Float64 = 0.0

    local_energy_psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    local_energy_psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    local_energy_psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    psi_derivative_a::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    psi_derivative_b::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    psi_derivative_w::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    for i = 1:num_mc_iterations

        if i > burn_in*num_mc_iterations

            # println("TEST1")

            metropolisStepBruteForce(nqs, step_length, precalc)
            # println("TEST2")

            precalc = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

            local_energy = computeLocalEnergy(nqs, precalc)
            # println("TEST3")

            computePsiParameterDerivative!(nqs, psi_derivative_a, psi_derivative_b, psi_derivative_w, precalc)
            # println("TEST4")

            local_energy_sum += local_energy

            local_energy_psi_derivative_a_sum += local_energy*psi_derivative_a
            local_energy_psi_derivative_b_sum += local_energy*psi_derivative_b
            local_energy_psi_derivative_w_sum += local_energy*psi_derivative_w

            psi_derivative_a_sum += psi_derivative_a
            psi_derivative_b_sum += psi_derivative_b
            psi_derivative_w_sum += psi_derivative_w

        end

    end

    samples = num_mc_iterations - burn_in*num_mc_iterations

    mc_local_energy = local_energy_sum/samples

    mc_local_energy_psi_derivative_a = local_energy_psi_derivative_a_sum/samples
    mc_local_energy_psi_derivative_b = local_energy_psi_derivative_b_sum/samples
    mc_local_energy_psi_derivative_w = local_energy_psi_derivative_w_sum/samples

    mc_psi_derivative_a = psi_derivative_a_sum/samples
    mc_psi_derivative_b = psi_derivative_b_sum/samples
    mc_psi_derivative_w = psi_derivative_w_sum/samples

    local_energy_derivative_a = 2*(mc_local_energy_psi_derivative_a - mc_local_energy*mc_psi_derivative_a)
    local_energy_derivative_b = 2*(mc_local_energy_psi_derivative_b - mc_local_energy*mc_psi_derivative_b)
    local_energy_derivative_w = 2*(mc_local_energy_psi_derivative_w - mc_local_energy*mc_psi_derivative_w)

    return mc_local_energy, local_energy_derivative_a, local_energy_derivative_b, local_energy_derivative_w

end

function runOptimizationBruteForce(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, mc_burn_in::Float64, mc_step_length::Float64, learning_rate::Float64)

    local_energies::Array{Float64, 2} = zeros(Float64, (num_iterations, 1))

    for k = 1:num_iterations
        local_energy,_grad_a,  _grad_b, _grad_w = runMetorpolisBruteForce(nqs, num_mc_iterations, mc_burn_in, mc_step_length)
        optimizationStep(nqs, _grad_a, _grad_b, _grad_w, learning_rate)
        local_energies[k] = local_energy
    end

    return local_energies

end

function computeDriftForce(nqs::NQS, particle_number::Int64, precalc::Array{Float64, 2})

    m = particle_number

    #Calculates the term that will be used multiple times to increase speed.
    # precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    drift_force = -(1.0/nqs.sigma_squared)*(nqs.x[m]- nqs.a[m])

    for n=1:num_hidden
        drift_force += (1.0/nqs.sigma_squared)*nqs.w[m,n]/(exp(-precalc[n]) + 1.0)
    end

    drift_force*=2.0

    return drift_force

end

function metropolisStepImportanceSampling(nqs::NQS, time_step::Float64, D::Float64, precalc)

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    coordinate::Int64 = rand(1:num_visible)

    old_wavefunction_value = computePsi(nqs, precalc)

    current_drift_force = computeDriftForce(nqs, coordinate, precalc)

    #Update the trial coordinate and store the current position:
    old_coordinate = nqs.x[coordinate]
    nqs.x[coordinate] += D*current_drift_force*time_step + randn(Float64)*sqrt(time_step)


    precalc = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    new_drift_force = computeDriftForce(nqs, coordinate, precalc)

    new_wavefunction_value = computePsi(nqs, precalc)

    greens_function_argument = (old_coordinate - nqs.x[coordinate] - D*time_step*new_drift_force)^2 - (nqs.x[coordinate] - old_coordinate - D*time_step*current_drift_force)^2

    greens_function_argument /= (4.0*D*time_step)

    greens_function = exp(-greens_function_argument)

    U = rand(Float64)

    if U < greens_function*(new_wavefunction_value^2)/(old_wavefunction_value^2)
        return
    else
        nqs.x[coordinate] = old_coordinate
        return
    end

end

function runMetropolisImportanceSampling(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, burn_in::Float64, importance_time_step::Float64, D::Float64)

    local_energy_sum::Float64 = 0.0

    local_energy_psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    local_energy_psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    local_energy_psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    psi_derivative_a::Array{Float64, 2} = zeros(Float64, size(nqs.a))
    psi_derivative_b::Array{Float64, 2} = zeros(Float64, size(nqs.b))
    psi_derivative_w::Array{Float64, 2} = zeros(Float64, size(nqs.w))

    precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    for i = 1:num_mc_iterations

        if i > burn_in*num_mc_iterations

            # println("TEST1")

            metropolisStepImportanceSampling(nqs, importance_time_step, D, precalc)

            # println("TEST2")
            precalc = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

            local_energy = computeLocalEnergy(nqs, precalc)

            computePsiParameterDerivative!(nqs, psi_derivative_a, psi_derivative_b, psi_derivative_w, precalc)

            local_energy_sum += local_energy

            local_energy_psi_derivative_a_sum += local_energy.*psi_derivative_a
            local_energy_psi_derivative_b_sum += local_energy.*psi_derivative_b
            local_energy_psi_derivative_w_sum += local_energy.*psi_derivative_w

            psi_derivative_a_sum += psi_derivative_a
            psi_derivative_b_sum += psi_derivative_b
            psi_derivative_w_sum += psi_derivative_w

            if i % 100 == 0
                # println("iter = ", i, "E_l = ", computeLocalEnergy(nqs))
            end
        end

    end

    samples = num_mc_iterations - burn_in*num_mc_iterations

    mc_local_energy = local_energy_sum/samples

    mc_local_energy_psi_derivative_a = local_energy_psi_derivative_a_sum/samples
    mc_local_energy_psi_derivative_b = local_energy_psi_derivative_b_sum/samples
    mc_local_energy_psi_derivative_w = local_energy_psi_derivative_w_sum/samples

    mc_psi_derivative_a = psi_derivative_a_sum/samples
    mc_psi_derivative_b = psi_derivative_b_sum/samples
    mc_psi_derivative_w = psi_derivative_w_sum/samples

    local_energy_derivative_a = 2*(mc_local_energy_psi_derivative_a - mc_local_energy*mc_psi_derivative_a)
    local_energy_derivative_b = 2*(mc_local_energy_psi_derivative_b - mc_local_energy*mc_psi_derivative_b)
    local_energy_derivative_w = 2*(mc_local_energy_psi_derivative_w - mc_local_energy*mc_psi_derivative_w)

    # println(mc_local_energy)

    return mc_local_energy, local_energy_derivative_a, local_energy_derivative_b, local_energy_derivative_w
end

function runOptimizationImportanceSampling(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, mc_burn_in::Float64, importance_time_step::Float64, D::Float64, learning_rate::Float64)

    local_energies::Array{Float64, 2} = zeros(Float64, (num_iterations, 1))

    for k = 1:num_iterations
        local_energy, _grad_a, _grad_b, _grad_w = runMetropolisImportanceSampling(nqs, num_iterations, num_mc_iterations, mc_burn_in, importance_time_step, D)
        optimizationStep(nqs, _grad_a, _grad_b, _grad_w, learning_rate)
        local_energies[k] = local_energy
        # println("Iteration = ", k, "    E = ", local_energy, nqs.h)
    end

    return local_energies
end

function metropolisStepGibbsSampling(nqs::NQS, precalc)

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    #Calculates the term that will be used multiple times to increase speed.

    mean = nqs.a + nqs.w*nqs.h
    variance = nqs.sigma_squared

    for i = 1:num_visible
        nqs.x[i] = sqrt(variance)*randn(Float64) + mean[i]
    end

    U = rand(Float64)

    # temp::Array{Float64,2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    for j = 1:num_hidden
        # temp::Float64 = nqs.b[j] + (1.0/nqs.sigma_squared)*dot(nqs.x,nqs.w[:,j])

        nqs.h[j] = U < 1.0/(exp(-precalc[j]) + 1.0)
    end

end

function runMetropolisGibbsSampling(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, burn_in::Float64)

        local_energy_sum::Float64 = 0.0

        local_energy_psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
        local_energy_psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
        local_energy_psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

        psi_derivative_a_sum::Array{Float64, 2} = zeros(Float64, size(nqs.a))
        psi_derivative_b_sum::Array{Float64, 2} = zeros(Float64, size(nqs.b))
        psi_derivative_w_sum::Array{Float64, 2} = zeros(Float64, size(nqs.w))

        psi_derivative_a::Array{Float64, 2} = zeros(Float64, size(nqs.a))
        psi_derivative_b::Array{Float64, 2} = zeros(Float64, size(nqs.b))
        psi_derivative_w::Array{Float64, 2} = zeros(Float64, size(nqs.w))

        precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))


        for i = 1:num_mc_iterations

            if i > burn_in*num_mc_iterations

                metropolisStepGibbsSampling(nqs, precalc)

                precalc = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))


                local_energy = computeLocalEnergyGibbs(nqs, precalc)

                computePsiParameterDerivative!(nqs, psi_derivative_a, psi_derivative_b, psi_derivative_w, precalc)

                psi_derivative_a *= 0.5
                psi_derivative_b *= 0.5
                psi_derivative_w *= 0.5

                local_energy_sum += local_energy

                local_energy_psi_derivative_a_sum += local_energy.*psi_derivative_a
                local_energy_psi_derivative_b_sum += local_energy.*psi_derivative_b
                local_energy_psi_derivative_w_sum += local_energy.*psi_derivative_w

                psi_derivative_a_sum += psi_derivative_a
                psi_derivative_b_sum += psi_derivative_b
                psi_derivative_w_sum += psi_derivative_w

                if i % 100 == 0
                    # println("iter = ", i, "E_l = ", computeLocalEnergy(nqs))
                end

            end

        end

        samples = num_mc_iterations - burn_in*num_mc_iterations

        mc_local_energy = local_energy_sum/samples

        mc_local_energy_psi_derivative_a = local_energy_psi_derivative_a_sum/samples
        mc_local_energy_psi_derivative_b = local_energy_psi_derivative_b_sum/samples
        mc_local_energy_psi_derivative_w = local_energy_psi_derivative_w_sum/samples

        mc_psi_derivative_a = psi_derivative_a_sum/samples
        mc_psi_derivative_b = psi_derivative_b_sum/samples
        mc_psi_derivative_w = psi_derivative_w_sum/samples

        local_energy_derivative_a = 2*(mc_local_energy_psi_derivative_a - mc_local_energy*mc_psi_derivative_a)
        local_energy_derivative_b = 2*(mc_local_energy_psi_derivative_b - mc_local_energy*mc_psi_derivative_b)
        local_energy_derivative_w = 2*(mc_local_energy_psi_derivative_w - mc_local_energy*mc_psi_derivative_w)

        return mc_local_energy, local_energy_derivative_a, local_energy_derivative_b, local_energy_derivative_w
end

function runOptimizationGibbsSampling(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, mc_burn_in::Float64, learning_rate::Float64)

    local_energies::Array{Float64, 2} = zeros(Float64, (num_iterations, 1))

    for k = 1:num_iterations
        local_energy, _grad_a, _grad_b, _grad_w = runMetropolisGibbsSampling(nqs::NQS, num_iterations::Int64, num_mc_iterations::Int64, mc_burn_in::Float64)
        optimizationStep(nqs, _grad_a, _grad_b, _grad_w, learning_rate)
        local_energies[k] = local_energy
        # println(local_energy)
    end

    return local_energies

end

"""
    grid_search_to_files()

Function for doing a grid search and writing the results from each run to files.
"""
function write_grid_search_to_files(num_particles::Int64, num_dims::Int64, method::String, interacting::Bool)
    learning_rates::Array{Float64,1} = [1.0, 0.1, 0.01]
    hidden_nodes::Array{Int64,1} = [2, 3, 4]

    len_learning_rates = length(learning_rates)
    len_hidden_nodes = length(hidden_nodes)

    M::Int64 = num_particles*num_dims          # Number of visible nodes

    mc_burn_in = 0.2

    mc_step_length = 0.5
    importance_sampling_step_length = 0.5

    num_mc_iterations = 100000
    num_optimization_steps = 200

    D = 0.5

    for i = 1:len_learning_rates
        for j = 1:len_hidden_nodes

            num_hidden_nodes = hidden_nodes[j]
            learning_rate = learning_rates[i]

            if method == "bf"
                sigma_squared = 1.0
                system = setUpSystemRandomUniform(num_particles, num_dims, M, num_hidden_nodes, sigma_squared, interacting)
                start = time()
                local_energies = runOptimizationBruteForce(system, num_optimization_steps, num_mc_iterations, mc_burn_in, mc_step_length, learning_rate)
                runtime = time() - start
                println(runtime)
                filename = string("output/interacting_", interacting, "/" , method, "/num_particles_",num_particles, "_num_dims_",num_dims , "_lr_", string(learning_rate), "_hidden_", string(num_hidden_nodes), "_sigma_", string(sigma_squared), "_mc_step_length_", string(mc_step_length), "_num_mc_iterations_", num_mc_iterations, ".txt")
                open(filename, "w") do io
                    println(io, "sigma=", sigma_squared, " time=", runtime, " mc_step_length=", mc_step_length, " num_mc_iterations=", num_mc_iterations)
                    # println(io, "TEST = ", local_energies[1])
                    for e in local_energies
                        println(io, e)
                    end
                end

            elseif method == "is"
                sigma_squared = 1.0
                system = setUpSystemRandomUniform(num_particles, num_dims, M, num_hidden_nodes, sigma_squared, interacting)
                start = time()
                local_energies = runOptimizationImportanceSampling(system, num_optimization_steps, num_mc_iterations, mc_burn_in, importance_sampling_step_length, D, learning_rate)
                runtime = time() - start
                println(runtime)
                filename = string("output/interacting_", interacting, "/" , method, "/num_particles_",num_particles, "_num_dims_",num_dims , "_lr_", string(learning_rate), "_hidden_", string(num_hidden_nodes),"_sigma_", string(sigma_squared), "_importance_step_length_", string(importance_sampling_step_length), "_num_mc_iterations_", num_mc_iterations, ".txt")
                open(filename, "w") do io
                    println(io, "sigma=", sigma_squared, " time=", runtime, " importance_sampling_step_length=", importance_sampling_step_length, " num_mc_iterations=", num_mc_iterations)
                    # println(io, "TEST = ", local_energies[1])
                    for e in local_energies
                        println(io, e)
                    end
                end


            elseif method == "gs"
                sigma_squared = 0.5
                system = setUpSystemRandomUniform(num_particles, num_dims, M, num_hidden_nodes, sigma_squared, interacting)
                start = time()
                local_energies = runOptimizationGibbsSampling(system, num_optimization_steps, num_mc_iterations, mc_burn_in, learning_rate)
                runtime = time() - start
                println(runtime)
                filename = string("output/interacting_", interacting, "/" , method, "/num_particles_",num_particles, "_num_dims_",num_dims , "_lr_", string(learning_rate), "_hidden_", string(num_hidden_nodes),"_sigma_", string(sigma_squared), "_num_mc_iterations_", num_mc_iterations, ".txt")
                open(filename, "w") do io
                    println(io, "sigma=", sigma_squared, " time=", runtime, " num_mc_iterations=", num_mc_iterations)
                    # println(io, "TEST = ", local_energies[1])
                    for e in local_energies
                        println(io, e)
                    end
                end
            end


        end
    end

end

function main()

    # SET UP THE SYSTEM AND BOLTZMAN MACHINE:
    num_particles = 1                          # Number of particles
    num_dims = 2                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 2                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    mc_burn_in = 0.2                           # Fraction of steps before sampling
    num_mc_cycles = 20000                    # Number of steps in the MC algorithm
    num_optimization_steps = 100                # Number of optimization steps

    brute_force_step_length = 1.0              # Step-length in the Brute-force Metropolis
    importance_sampling_step_length = 0.05     # Time-step in the importance sampling algorithm
    learning_rate = 0.05                        # Learning rate in the optimization algorithm
    D = 0.5                                    # Diffusion constant for importance sampling

    # system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    # system_gibbs = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared_gibbs, interacting)

    # @time runOptimizationBruteForce(system, num_optimization_steps, num_mc_cycles, mc_burn_in, brute_force_step_length, learning_rate)
    # @time runOptimizationImportanceSampling(system, num_optimization_steps, num_mc_cycles, mc_burn_in, importance_sampling_step_length, D, learning_rate)
    # @time runOptimizationGibbsSampling(system_gibbs, num_optimization_steps, num_mc_cycles, mc_burn_in, learning_rate)

    write_grid_search_to_files(num_particles, num_dims, "is", interacting)
end

main()

end
