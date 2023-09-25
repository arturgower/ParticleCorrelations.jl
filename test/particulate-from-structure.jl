# Here we calculate a particle configuration from a structure factor

using ParticleCorrelations
using Test, Statistics, LinearAlgebra

using Plots

# let us first aim to recover a MonteCarlo pair-correlation from a crystalline arrangement

# choose the spatial dimension
Dim = 2

# choose the medium for the particles. Currently only one type
medium = HardMedium{Dim}()

# choose the shapes of the particles

pairtype = MonteCarloPairCorrelation(Dim; 
    rtol = 1e-3, 
    maxlength = 100, 
    iterations = 10, 
    numberofparticles = 3000
)

# choose the medium for the particles. Currently only one type
medium = HardMedium{Dim}()

# Choose the species, which represents a collection of one type of particle
radius = 0.5

specie = Specie(
    medium,
    Sphere(Dim, radius),
    volume_fraction = 0.15,
    separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
);

specie.particle.medium

rs = 0.2:0.4:8.0
pair = pair_correlation(specie, pairtype, rs)


# If you have the Plots package
plot(pair.r, pair.g)

ks = 1.0:0.5:50.0
ks = 1.0:0.3:40.0
ks2 = 1.0:0.2:160.0
sfactor = structure_factor(pair, ks)
sfactor2 = structure_factor(pair, ks2)

plot(sfactor.k, sfactor.S)
plot!(sfactor2.k, sfactor2.S)

using Optim

# target_S = sfactor
# Dim = 2
# numberofparticles = 200
# correlation_length = 4

res = optimise_particulate(sfactor, specie;
    numberofparticles = 300,
    optimoptions = Optim.Options(g_tol = 0.0, f_tol = 1e-5, iterations = 20)
)

show(res)

x = res.minimizer
points = Iterators.partition(x,Dim) |> collect
x1 = [p[1] for p in points]
x2 = [p[2] for p in points]

scatter(x1,x2)
# scatter(x1,x2, xlims =(-30.0,30.0), ylims = (-30.0,30.0))

outer_box = Box(points);
inner_box = Box(outer_box.origin, outer_box.dimensions .- 8.0);

result_S = DiscreteStructureFactor(points, ks; 
    inner_box = inner_box
)
norm(result_S.S - sfactor.S)

plot(result_S.k, result_S.S)
plot!(sfactor.k, sfactor.S, linestyle = :dash)

result_S = DiscreteStructureFactor(points, ks2; 
    inner_box = inner_box
)

plot(result_S.k, result_S.S)
plot!(sfactor2.k, sfactor2.S, linestyle = :dash)


result_pair = DiscretePairCorrelation(points, pair.r)
plot(pair.r, pair.g)
plot!(result_pair.r, result_pair.g)


points = origin.(particles)
x = vcat(points...)

points = Iterators.partition(x,Dim) |> collect
points1 = filter(p -> p ∈ reg1, points)
J1 = length(points1)  

S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)
dS = S.S - target_S.S # add quadradture weight w here if wanted


function objective(x)
    points = Iterators.partition(x,Dim) |> collect
    S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)

    return sum(abs2.(S.S - target_S.S))
end

ks = abs.(target_S.k)



G = vcat(map(points) do p
    sum(
        if p1 == p
            zeros(Float64, Dim)
        else     
            Rji = p - p1
            nRji = norm(Rji)
            
            Gji = (p ∈ reg1 ? 2 : 1) * sum(diffbesselj.(0, ks .* nRji) .* dS) / J1
            2 * Gji * Rji / nRji
        end    
    for p1 in points1)
end...)


DfDrs = Iterators.partition(G,Dim) |> collect

plot(particles)
plot!(reg1)

x1 = [p[1] for p in points];
x2 = [p[2] for p in points];

u = [p[1] for p in DfDrs];
v = [p[2] for p in DfDrs];

scale = 6
quiver!(x1,x2,quiver=(scale .* u,scale .* v))
plot!()

i_max = 40
xs = map(1:i_max) do i
    x - G .* 2i/i_max 
end    

fs = objective.(xs)

plot(fs)

## An arteficial pair-correlation like hyper-uniform disorder.
a = 1.0
rs = 0.2:0.4:8.0

g = 1.0 .+ exp.(-(rs .- 4a).^2)

dpair = DiscretePairCorrelation(2, rs, g)

number_density(dpair)

plot(rs,g)

