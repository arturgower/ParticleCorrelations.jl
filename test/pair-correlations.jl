using ParticleCorrelations
using LinearAlgebra, Statistics, Test
using CSV

@testset "Percus-Yevick" begin

    # Load reference data
    stdir = if intersect(readdir(),["test"]) |> isempty
        ""
    else "test/"
    end

    # The data assume R = 1, where R is the minimum distance between the centres of any two spheres.
    file = CSV.File("$(stdir)data/P-Y_f=0.25.txt"; header = 1)
    distances = file.r
    g_reference = file.field 

    # NOTE: increasing the exclusion_distance would increase the effective volume fraction needed for PY. The saved data assumes and effective volume fraction of 25%

    φ = 0.25
    r = 0.9
    γ = 1.5

    # generate data from function in package
    s1 = Specie(
        HardMedium{3}(),
        Sphere(r),
        volume_fraction = φ / γ^3,
        separation_ratio = γ  # <---- changes the effective volume fraction
    );

    R = 2*outer_radius(s1) * s1.separation_ratio

    pairtype = PercusYevick(3; rtol=1e-3, maxevals = Int(2e5))

    # Need to scale the distance by R
    py = DiscretePairCorrelation(s1, pairtype, R .* distances)

    i = findfirst(distances .> 1.0)

    @test norm(py.g[i:end] - g_reference[i:end]) / norm(g_reference[i:end]) < 0.02
    @test abs(py.g[i+1] - g_reference[i+1]) / norm(g_reference[i+1]) < 0.01

    # calculate the structure factor
    # Need to change the mesh of r
    py = DiscretePairCorrelation(s1, pairtype)

    py_factor = structure_factor(py) 

    @test true
end

@testset "structure-factor" begin

    # use a pair correlation with a known exact formula for the structure factor
    α = 1.5; 
    β = 3.0;

    number_density = 0.5
    dim = 3

    f(r) = 1 + (exp(im * β * r) + exp(-im * β * r)) * exp(- α * r) / r

    rs = 0.001:0.001:10.0;

    pair = DiscretePairCorrelation(dim, rs, f.(rs); number_density = 0.5)

    @test pair.g[end] - 1.0 < 1e-8

    sfactor = structure_factor(pair)

    ks = sfactor.k

    # using Plots 

    # plot(ks, sfactor.S, xlims = (0.,20))

    S(k) = 1 + 4π * number_density * real(1 / (k^2 + (α - im * β)^2) + 1 / (k^2 + (α + im * β)^2) )

    @test norm(S.(ks) - sfactor.S) / norm(sfactor.S) < 1e-4
    # plot!(ks, S.(ks))


end  