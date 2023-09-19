# Here we calculate a particle configuration from a structure factor

using ParticleCorrelations
using Test, Statistics

using Optim
using Plots


a = 1.0
rs = 0.0:0.01:8.0

g = 1.0 .+ exp.(-(rs .- 4a).^2)

dpair = DiscretePairCorrelation(2, rs, g)

number_density(dpair)


plot(rs,g)

