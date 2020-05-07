module Test

using LinearAlgebra

struct Hamiltonian
    localEnergy::Array{Float64, 1}
end

struct NQS
    b::Array{Float64, 2}
    a::Array{Float64, 2}
    w::Array{Float64, 2}
    x::Array{Float64, 2}
    h::Array{Float64, 2}
    sigma_squared::Float64
    interacting::Bool

end

"""
    computeLocalEnergy(nqs::NQS)

Computes the local energy of the system with the given
hamiltonian and wavefunction.
"""
function computeLocalEnergy(nqs::NQS)
    #Calculates the term that will be used multiple times to increase speed.
    precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    local_energy::Float64 = 0

    for m = 1:num_visible

        ln_psi_derivative = -(1.0/nqs.sigma_squared)*(nqs.x[m] - nqs.a[m])
        ln_psi_double_derivative = -1.0/nqs.sigma_squared

        for n=1:num_hidden

            ln_psi_derivative += (1.0/nqs.sigma_squared)*nqs.w[m,n]/(exp(-precalc[n]) + 1.0)
            ln_psi_double_derivative += (1.0/nqs.sigma_squared^2)*exp(precalc[n])/((exp(precalc[n])+1)^2)

        end

        local_energy += -(ln_psi_derivative)^2 - ln_psi_double_derivative + nqs.x[m]^2

    end

    local_energy *= 0.5

    if nqs.interacting
        local_energy += 100.
    end

    return local_energy

end

function computePsiParameterDerivative(nqs::NQS)

    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    #Calculates the term that will be used multiple times to increase speed.
    precalc::Array{Float64, 2} = nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w))

    psi_derivative_parameter::Array{Float64, 2} = zeros(Float64, num_visible + num_hidden + num_visible*num_hidden ,1)

    psi_derivative_parameter[1:num_visible, 1] = (1.0/nqs.sigma_squared)*(nqs.x - nqs.a)
    psi_derivative_parameter[num_visible + 1: num_visible + num_hidden + 1, 1] =

    println(psi_derivative_parameter)
end

function setUpSystemRandom(M::Int64, N::Int64, sig_sq::Float64 = 0.5, inter::Bool = false)
    b = rand(Float64, N, 1)
    a = rand(Float64, M, 1)
    w = rand(Float64, M, N)
    x = rand(Float64, M, 1)
    h = rand(Float64, N, 1)
    interacting = inter
    return NQS(b, a, w, x, h, sig_sq, interacting)
end

function main()
    M::Int64 = 10
    N::Int64 = 10
    sigma_squared = 0.5
    interacting = false

    test = setUpSystemRandom(M, N, sigma_squared, interacting)

    print(computeLocalEnergy(test))
    computePsiParameterDerivative(test)
end

main()

end
