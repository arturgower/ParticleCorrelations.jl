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


@testset "disordered particulate" begin

    # choose the spatial dimension
    dim = 2

    # choose the medium for the particles. Currently only one type
    medium = HardMedium{dim}()

    # choose the shapes of the particles
    radius = 0.5
    particle_shapes = [Sphere(dim,radius)]

    # choose a region to place the particles within
    dimensions = repeat([90.0],dim)
    region_shape = Box(dimensions)

    # create a uniform random arrangement of particles using Sequential Addition 
    particles = random_particles(medium, particle_shapes;
        seed = 3,
        volume_fraction = 0.1,
        region_shape = region_shape
    )

    correlation_length = 5.0
    ks = 1.5:0.1:10.0

    sfactor = structure_factor(particles, ks; correlation_length = correlation_length)

    # using MultipleScattering
    # plot(time_to_frequency(sfactor.S, sfactor.k) |> real)

    # ω = t_to_ω(sfactor.k)

    # gs = time_to_frequency(sfactor.S, sfactor.k);
    # gs[20:end] .= 0.0

    # plot(ω_to_t(ω), frequency_to_time(gs, ω))
    # plot!(sfactor.k, sfactor.S)

    # using Plots
    # plot(sfactor.k, sfactor.S)

    rs = (2radius):0.001:correlation_length
    pair = pair_correlation(particles, rs)

    # If you have the Plots package
    # plot(pair.r, pair.g)

    sfactor2 = structure_factor(pair, ks)

    # plot(sfactor.k, sfactor.S)
    # plot!(sfactor2.k, sfactor2.S, linestyle = :dash)

    @test norm(sfactor.S - sfactor2.S) / norm(sfactor.S) < 0.1

    # choose the spatial dimension
    dim = 3

    # choose the medium for the particles. Currently only one type
    medium = HardMedium{dim}()

    # choose the shapes of the particles
    particle_shapes = [Sphere(dim,radius)]

    # choose a region to place the particles within
    dimensions = repeat([25.0],dim)
    region_shape = Box(dimensions)

    # create a uniform random arrangement of particles using Sequential Addition 
    particles = random_particles(medium, particle_shapes;
        seed = 3,
        volume_fraction = 0.1,
        region_shape = region_shape
    )

    correlation_length = 5.0
    ks = 2.0:0.1:10.0

    sfactor = structure_factor(particles, ks; correlation_length = correlation_length)

    rs = (2radius):0.01:correlation_length
    pair = pair_correlation(particles, rs)

    # If you have the Plots package
    # plot(pair.r, pair.g)

    sfactor2 = structure_factor(pair, ks)

    # plot(sfactor.k, sfactor.S)
    # plot!(sfactor2.k, sfactor2.S, linestyle = :dash)

    @test norm(sfactor.S - sfactor2.S) / norm(sfactor.S) < 0.05

end  