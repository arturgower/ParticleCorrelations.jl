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
end