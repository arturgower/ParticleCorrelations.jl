# Here we calculate a particle configuration from a structure factor

using ParticleCorrelations
using Test, Statistics, LinearAlgebra

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

specie = Specie(
    medium,
    Sphere(dim, radius),
    volume_fraction = 0.15,
    separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
);

specie.particle.medium

rs = 0.2:0.4:8.0
pair = pair_correlation(specie, pairtype, rs)


# If you have the Plots package
plot(pair.r, pair.g)

ks = 1.0:0.5:20.0
sfactor = structure_factor(pair, ks)

plot(sfactor.k, sfactor.S)

using Optim


target_S = sfactor
Dim = 2
numberofparticles = 200
correlation_length = 4

function optimise_particulate(target_S::DiscreteStructureFactor{Dim}, specie::Specie{Dim};
        correlation_length::Float64 = 4.0, 
        numberofparticles::Int = 200) where Dim

    # Define two regions R1 and R2 to place the particles
    cell_number = (numberofparticles)^(1/Dim) |> round

    cell_volume = 1 / number_density(specie);
    cell_length = cell_volume ^ (1/Dim)
    box_length = cell_length * cell_number

    reg1 = Box([
        [-box_length/2, -box_length/2],
        [box_length/2, box_length/2]
    ]);

    cell_number_boundary = ceil(correlation_length/(2 * cell_length))
    boundary_length = cell_number_boundary * cell_length

    reg2 = Box([
        [-box_length/2 - boundary_length, -box_length/2 - boundary_length],
        [box_length/2 + boundary_length, box_length/2 + boundary_length]
    ]);

    particles = periodic_particles(reg2,specie; random_perturbation = true)
    particles1 = filter(p -> p ⊆ reg1, particles)  

    S = structure_factor(particles,ks; inner_box = reg1)

    points = origin.(particles)
    x0 = vcat(points...)

    function objective(x)
        points = Iterators.partition(x,Dim) |> collect
        S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)

        return sum(abs2.(S.S - target_S.S))
    end

    ks = abs.(target_S.k)

    function g!(G, x)
        points = Iterators.partition(x,Dim) |> collect
        points1 = filter(p -> p ∈ reg1, points)
        J1 = length(points1)  

        S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)
        dS = S.S - target_S.S # add quadradture weight w here if wanted

        G[:] = vcat(map(points) do p
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
    end    
    
    optimize(objective, g!, x0, LBFGS())
end



points = Iterators.partition(x,Dim) |> collect
points1 = filter(p -> p ∈ reg1, points)
J1 = length(points1)  

S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)
dS = S.S - target_S.S # add quadradture weight w here if wanted

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

x = [p[1] for p in points];
y = [p[2] for p in points];

u = [p[1] for p in DfDrs];
v = [p[2] for p in DfDrs];

scale = 6
quiver!(x,y,quiver=(scale .* u,scale .* v))
plot!()

g = 1.0 .+ exp.(-(rs .- 4a).^2)

dpair = DiscretePairCorrelation(2, rs, g)

number_density(dpair)


plot(rs,g)

