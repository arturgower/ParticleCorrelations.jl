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

    ks = 0.2:0.2:20.0 
    sfactor = structure_factor(pair, ks)

    S(k) = 1 + 4π * number_density * real(1 / (k^2 + (α - im * β)^2) + 1 / (k^2 + (α + im * β)^2) )

    @test norm(S.(ks) - sfactor.S) / norm(sfactor.S) < 1e-5

end  