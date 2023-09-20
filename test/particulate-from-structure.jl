# Here we calculate a particle configuration from a structure factor

using ParticleCorrelations
using Test, Statistics

using Plots

# let us first aim to recover a MonteCarlo pair-correlation from a crystalline arrangement

# choose the spatial dimension
dim = 2

# choose the medium for the particles. Currently only one type
medium = HardMedium{dim}()

# choose the shapes of the particles

pairtype = MonteCarloPairCorrelation(dim; 
    rtol = 1e-3, 
    maxlength = 100, 
    iterations = 10, 
    numberofparticles = 3000
)

# choose the medium for the particles. Currently only one type
medium = HardMedium{dim}()

# Choose the species, which represents a collection of one type of particle
radius = 0.5

s = Specie(
    medium,
    Sphere(dim, radius),
    volume_fraction = 0.15,
    separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
);

rs = 0.2:0.4:8.0
pair = pair_correlation(s, pairtype, rs)


# If you have the Plots package
plot(pair.r, pair.g)

ks = 1.0:0.5:20.0
sfactor = structure_factor(pair, ks)

plot(sfactor.k, sfactor.S)


# Define two regions R1 and R2 to place the particles
reg1 = Box([[-40.0, -40.0],[40.0, 40.0]]);
reg2 = Box([[-42.0, -42.0],[42.0, 42.0]]);

Dim = 2 

# box = reg2
function periodic_particles(box::Box{T,Dim}, specie::Specie{Dim}) where {T, Dim}
    
    Id = Matrix(I, Dim, Dim)
    vs = [Id[:,i] for i = 1:Dim]

    n = number_density(specie)

    cell_length = (1 / n)^(1/Dim) / 2

    d1 = origin(box) - box.dimensions ./ (2)
    d2 = origin(box) + box.dimensions ./ (2)
    
end    

using Optim

function g!(G, S)
    G[1] = -2.0 * (1.0 - S[1]) 
end

g = 1.0 .+ exp.(-(rs .- 4a).^2)

dpair = DiscretePairCorrelation(2, rs, g)

number_density(dpair)


plot(rs,g)

