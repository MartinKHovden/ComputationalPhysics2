include("library.jl")
using .library
using Test
using Random

@testset "test setUpSystemRandomUniform 1 particle in 2 dims" begin
    num_particles = 1                          # Number of particles
    num_dims = 2                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)

    @test system.num_particles == num_particles
    @test system.num_dims == num_dims
    @test system.sigma_squared == sigma_squared
    @test system.interacting == interacting
    @test size(system.a)[1] == M
    @test size(system.b)[1] == N
    @test size(system.x)[1] == M
    @test size(system.h)[1] == N
    @test size(system.a)[2] == 1
    @test size(system.b)[2] == 1
    @test size(system.x)[2] == 1
    @test size(system.h)[2] == 1
    @test size(system.w) == (M,N)

end

@testset "test setUpSystemRandomUniform 10 particles in 3 dims" begin
    num_particles = 10                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)

    @test system.num_particles == num_particles
    @test system.num_dims == num_dims
    @test system.sigma_squared == sigma_squared
    @test system.interacting == interacting
    @test size(system.a)[1] == M
    @test size(system.b)[1] == N
    @test size(system.x)[1] == M
    @test size(system.h)[1] == N
    @test size(system.a)[2] == 1
    @test size(system.b)[2] == 1
    @test size(system.x)[2] == 1
    @test size(system.h)[2] == 1
    @test size(system.w) == (M,N)


end

@testset "test computeInteractionTerm" begin
    Random.seed!(123)
    num_particles = 2                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)

    system.x[:, 1] = randn(M)*0.01
    @test abs(computeInteractionTerm(system) - 29.7155) <= 0.01

    system.x[:,1] = [1;0;0; 0;0;0]
    @test computeInteractionTerm(system) ≈ 1

    system.x[:,1] = [0;1;0; 0;0;0]
    @test computeInteractionTerm(system) ≈ 1

    system.x[:,1] = [1;1;0; 1;1;1]
    @test computeInteractionTerm(system) ≈ 1

    system.x[:,1] = [0;0;0; 0;0;0]
    @test computeInteractionTerm(system) == Inf

    system.x[:,1] = [1;0;0; 1;0;0]
    @test computeInteractionTerm(system) == Inf

    system.x[:,1] = [1;2.3;0; 1;2.3;0]
    @test computeInteractionTerm(system) == Inf

    for i = 1:10
        system.x[:,1] = [i/10.0, i/10.0, i/10.0, i/10.0, i/10.0, i/10.0]
        @test computeInteractionTerm(system) == Inf
    end

    for i = 1:10
        system.x[:,1] = [i/10.0+1, i/10.0+2, i/10.0+3, i/10.0+1, i/10.0+2, i/10.0+3]
        @test computeInteractionTerm(system) == Inf
    end

    system.x[:,1] = [1;0;0; 0.5;0;0]
    @test computeInteractionTerm(system) ≈ 2

    system.x[:,1] = [0.2;0;0; 0.7;0;0]
    @test computeInteractionTerm(system) ≈ 2




end

@testset "test computeLocalEnergy" begin
    Random.seed!(123)
    num_particles = 2                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    precalc::Array{Float64, 2} = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computeLocalEnergy(system, precalc)) ≈ 2.2086045478681564

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    precalc = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computeLocalEnergy(system, precalc)) ≈ 1.2550477488718554

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, true)
    precalc = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computeLocalEnergy(system, precalc)) ≈ 3.086798612593536

end

@testset "test computePsi" begin
    Random.seed!(123)
    num_particles = 2                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    precalc::Array{Float64, 2} = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computePsi(system, precalc)) ≈ 11.225542310264295

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    precalc = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computePsi(system, precalc)) ≈ 4.487620906057123

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, true)
    precalc = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test (computePsi(system, precalc)) ≈ 14.980332190247944
end

@testset "test optimizationStep" begin
    Random.seed!(123)
    num_particles = 2                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    system.w[:,:] = ones(M, N)
    system.a[:,1] = ones(M,1)
    system.b[:,1] = ones(N,1)

    grad_w = 0.5*ones(M,N)
    grad_a = 0.5*ones(M,1)
    grad_b = 0.5*ones(N,1)

    optimizationStep(system, grad_a, grad_b, grad_w, 2.0)

    @test system.w ≈ zeros(M,N)
    @test system.a ≈ zeros(M,1)
    @test system.b ≈ zeros(N,1)

    system.w[:,:] = ones(M, N)
    system.a[:,1] = ones(M,1)
    system.b[:,1] = ones(N,1)

    grad_w = 0.5*ones(M,N)
    grad_a = 0.5*ones(M,1)
    grad_b = 0.5*ones(N,1)

    optimizationStep(system, grad_a, grad_b, grad_w, 1.0)

    @test system.w ≈ 0.5*ones(M,N)
    @test system.a ≈ 0.5*ones(M,1)
    @test system.b ≈ 0.5*ones(N,1)

end

@testset "test computeDriftForce" begin
    Random.seed!(123)
    num_particles = 2                          # Number of particles
    num_dims = 3                               # Number of dimensions
    M::Int64 = num_particles*num_dims          # Number of visible nodes
    N::Int64 = 4                               # Number of hidden nodes
    sigma_squared = 1.0                        # RBM variance
    sigma_squared_gibbs = 0.5
    interacting = false                        # Interacting system?

    system = setUpSystemRandomUniform(num_particles, num_dims, M, N, sigma_squared, interacting)
    precalc::Array{Float64, 2} = system.b + transpose((1.0/system.sigma_squared)*(transpose(system.x)* system.w))

    @test computeDriftForce(system, 1, precalc) ≈ -1.293323006602911
    @test computeDriftForce(system, 2, precalc) ≈  1.8481682201086667

end
