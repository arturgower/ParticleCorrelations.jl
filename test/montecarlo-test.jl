using ParticleCorrelations
using LinearAlgebra, Statistics, Test
using CSV

@testset "Monte-Carlo" begin

    dim = 3;

    # choose the medium for the particles. Currently only one type
    medium = HardMedium{dim}()

    # choose the particle radius
    radius = 0.5

    # Choose the species, which represents a collection of one type of particle
    s = Specie(
        medium,
        Sphere(dim, radius),
        volume_fraction = 0.15,
        separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
    );

    R = 2*outer_radius(s) * s.separation_ratio

    pairtype = PercusYevick(3; rtol=1e-3, maxevals = Int(2e5))
    # Need to scale the distance by R
    rs = (1.05:0.1:4.05) .* R

    pc = pair_correlation(s, pairtype,rs)

    # If you have the Plots package
        # plot(pc.r, pc.g)
    
    # The Monte-Carlo method results in a slightly lower volume fraciton, as particles concetrate near the boundary which are then ignored.
    s_mc = Specie(
        medium,
        Sphere(dim, radius),
        volume_fraction = 0.15 * 1.15,
        # volume_fraction = 0.15 * 1.4,
        separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
    );

    # decreasing the meshsize, or increasing the iterations or numberofparticles, does not seem to improve the match with PY.
    pairtype_mc = MonteCarloPairCorrelation(dim; 
        rtol = 1e-3, 
        maxlength = 20, 
        meshsize = 0.2, 
        # iterations = 200, 
        iterations = 100, 
        # numberofparticles = 3000
        numberofparticles = 3000
    )

    pc_mc = pair_correlation(s_mc, pairtype_mc, rs)
        # scatter!(pc_mc.r, pc_mc.g, 
        #     ylims = (0.9,1.55), 
        #     xlims = (0.9,3.0)
        # )

    inds = findall(1.0 .< pc_mc.r .< 2.5)
    pc_mc.r[inds]
    pc.r[inds]

    # The agreement between PY and MC is not as good as shown in 
    #  Figure 8.3.1 from [1]. Though figure 8.3.7 which uses the sequential addition method for the MC, like we do, does show that MC underestimates PY.
    # [1] Kong, Jin Au, Leung Tsang, Kung-Hau Ding, and Chi On Ao. Scattering of electromagnetic waves: numerical simulations. John Wiley & Sons, 2004.
    @test mean(abs.(pc_mc.g[inds] - pc.g[inds])) .< 0.03


    # if we use a far higher volume fraction for the MC then it matches PY better.
    s_mc = Specie(
        medium,
        Sphere(dim, radius),
        volume_fraction = 0.15 * 1.38,
        separation_ratio = 1.0 
    );
    pc_mc = pair_correlation(s_mc, pairtype_mc, rs)

    # plot(pc.r, pc.g)

    # scatter!(pc_mc.r, pc_mc.g, 
    #     ylims = (0.9,1.55), 
    #     xlims = (0.9,3.0)
    # )

    @test mean(abs.(pc_mc.g[inds] - pc.g[inds])) < 0.015
end