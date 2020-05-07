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

end

"""
    computeLocalEnergy(nqs::NQS)

Computes the local energy of the system with the given
hamiltonian and wavefunction. 
"""
function computeLocalEnergy(nqs::NQS)
    #Calculates the term that will be used multiple times to increase speed.
    precalc::Array{Float64, 2} = exp.(nqs.b + transpose((1.0/nqs.sigma_squared)*(transpose(nqs.x)* nqs.w)))

    #Extracts the number of hidden and visible units.
    num_visible = size(nqs.x)[1]
    num_hidden = size(nqs.h)[1]

    for m = 1:num_visible
        ln_psi_derivative =
    end

    return precalc

end
#
# function computeDoubleDerivative(nqs::NQS, particles::Array{Float64, 2})
#
# end
#
# function setLocalEnergy(locE::Float64, ham::Hamiltonian)
#     ham.localEnergy[1] = locE
# end


function main()
    M::Int64 = 10
    N::Int64 = 10
    b = rand(Float64, N, 1)
    a = rand(Float64, M, 1)
    w = rand(Float64, M, N)
    x = rand(Float64, M, 1)
    h = rand(Float64, N, 1)
    print(transpose(b)* a)
    test = NQS(b, a, w, x, h, 0.5)

    print(computeLocalEnergy(test))
end

main()

end
