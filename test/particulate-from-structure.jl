# Here we calculate a particle configuration from a structure factor

using ParticleCorrelations
using Test, Statistics, LinearAlgebra

using Plots, Optim

# let us first aim to recover a MonteCarlo pair-correlation from a crystalline arrangement

# choose the spatial dimension
Dim = 2

# choose the medium for the particles. Currently only one type
medium = HardMedium{Dim}()

# choose the shapes of the particles

pairtype = MonteCarloPairCorrelation(Dim; 
    rtol = 1e-3, 
    maxlength = 100, 
    iterations = 800, 
    numberofparticles = 4000
)

# pairtype = MonteCarloPairCorrelation(Dim; 
#     rtol = 1e-3, 
#     maxlength = 100, 
#     iterations = 1, 
#     numberofparticles = 10000
# )

# choose the medium for the particles. Currently only one type
medium = HardMedium{Dim}()

# Choose the species, which represents a collection of one type of particle
radius = 0.5

specie = Specie(
    medium,
    MultipleScattering.Sphere(Dim, radius),
    volume_fraction = 0.15,
    separation_ratio = 1.0 # minimal distance from this particle = r * (separation_ratio - 1.0) 
);
s = specie

specie.particle.medium

dr = 0.05
rs = (2radius + dr/2):dr:4.0
distances = rs
T = Float64

pair = pair_correlation(specie, pairtype, rs)

# check the restriction from prob. dist.
σ = trapezoidal_scheme(pair.r)

sum(σ .* pair.g .* pair.r)
sum(dr .* pair.g .* pair.r)

pair.r[end]^2 / 2 - 1 / (2pi * pair.number_density)

(pair.r[end] + dr/2)^2 / 2 - 1 / (2pi * pair.number_density)

# If you have the Plots package

h = 220
gr(size = (h * 1.6, h), linewidth = 2.0)

plot(pair.r, pair.g, 
    title = "Percus-Yevick",
    ylab = "pair-correlation", xlab = "z (distance)",
    lab = ""
)
# savefig("percus-yevick-pair.pdf")


using Interpolations

extrap = LinearInterpolation(pair.r[2:end], pair.g[2:end], extrapolation_bc = Line())


# itp = interpolate(pair.g[2:end], BSpline(Cubic(Line(OnGrid()))))
# itp = interpolate(pair.g[2:end], BSpline(Linear()))

# sitp = scale(itp, rs[2:end])

rs2 = (2radius):0.01:rs[end]
g2 = extrap.(rs2)

plot!(rs2, g2, linestyle = :dash, xlims = (0.95,3.5))

σ2 = trapezoidal_scheme(rs2)

sum(σ2 .* g2 .* rs2)


# rs = 0.2:0.4:8.0

D1 = 40; 
k1 = 4pi /D1;
k2 = 2pi /radius;

ks = k1:k1:(2k2)
ks2 = k1:(k1/2):(3k2)

ks = 0.1:0.3:40.0
ks2 = 0.1:0.2:50.0

sfactor = structure_factor(pair, ks)
sfactor2 = structure_factor(pair, ks2)

h = 230
pyplot(size = (h*1.6,h))
plot(sfactor.k, sfactor.S,
    title = "Percus-Yevick",
    ylab = "structure factor", xlab = "k (wavenumber)",
    lab = ""
)
# savefig("percus-yevick-S.pdf")


plot!(sfactor2.k, sfactor2.S, linestyle = :dash)

# target_S = sfactor
Dim = 2
numberofparticles = 600
correlation_length = 4
target_S = sfactor

res1, res = optimise_particulate(sfactor, specie;
    numberofparticles = numberofparticles,
    method = LBFGS(),
    optimoptions = Optim.Options(g_tol = 0.0, f_tol = 1e-10, iterations = 15)
    # optimoptions = Optim.Options(method = LBFGS(), g_tol = 0.0, f_tol = 1e-10, iterations = 20)
)

show(res)

x = res1.minimizer
points1 = Iterators.partition(x,Dim) |> collect
particles1 = [Particle(medium,congruent(specie.particle.shape,p)) for p in points1]

pyplot(size = (420,420), linewidth = 1.0)
plot(particles1)

x = res.minimizer
points = Iterators.partition(x,Dim) |> collect
particles = [Particle(medium,congruent(specie.particle.shape,p)) for p in points]

plot(particles)


outer_box = Box(points);
inner_box = Box(outer_box.origin, outer_box.dimensions .- 8.0);

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


inner_box = reg1


pyplot(size = (420,420), linewidth = 1.0)
plot(particles, title="Predicted particle configuration",  
    xguide = "", yguide =""
    # xlims = (-40,40), ylims = (-30,30)
)
plot!(reg1,  xlab = "", ylab = "",axis = false) 
# savefig("predicted_particles.pdf")   
# scatter(x1,x2, xlims =(-30.0,30.0), ylims = (-30.0,30.0))


particles0 = periodic_particles(reg2,specie; random_perturbation = true)

pyplot(size = (420,420), linewidth = 1.0)
plot(particles0, title="Initial particle configuration",  
    xguide = "", yguide ="", 
    # xlims = (-40,40), ylims = (-30,30)
    )
plot!(reg1,  xlab = "", ylab = "",axis = false) 
# savefig("initial_particles.pdf") 


result1_S = DiscreteStructureFactor(points1, ks; 
    inner_box = inner_box
)
sum(abs2.(result1_S.S - sfactor.S)) / sum(abs2.(sfactor.S))

result_S = DiscreteStructureFactor(points, ks; 
    inner_box = inner_box
)
sum(abs2.(result_S.S - sfactor.S)) / sum(abs2.(sfactor.S))

h = 220
gr(size = (h * 1.6, h), linewidth = 2.0)
plot(result1_S.k, result1_S.S)
plot(sfactor.k, sfactor.S, label = "target")
plot!(result_S.k, result_S.S, 
    linestyle = :dash, label = "fitted")
plot!(ylab = "Structure factor", xlab = "k (wavenumber)")
# savefig("fitter-structure-factor.pdf")

# plot!(xlims = (0.0,20.0), ylims = (0.0,2.0))

result_S = DiscreteStructureFactor(points, ks2; 
    inner_box = inner_box
)

plot(result_S.k, result_S.S)
plot!(sfactor2.k, sfactor2.S, linestyle = :dash, 
    ylims = (0,1.4), xlims = (0.01,20))

rs = pair.r
rs = 0.2:0.4:12.0
rs = 0.2:0.35:10.0
result_pair = DiscretePairCorrelation(points, rs)
result_pair1 = DiscretePairCorrelation(points1, rs)

plot(pair.r, pair.g, label = "target")
plot!(rs, result_pair1.g, linestyle = :dash, label = "fitted")
plot!(ylab = "Pair correlation", xlab = "z (distance)")
# savefig("fit-pair-correlation.pdf")

plot!(rs, result_pair.g, linestyle = :dash)


# Check gradient and steps

box = Box([[-20.0,-20.0],[20.0,20.0]])
reg1 = Box([[-17.0,-17.0],[17.0,17.0]])

particles = periodic_particles(box, specie; random_perturbation = true);
points = origin.(particles);
x1 = [p[1] for p in points]
x2 = [p[2] for p in points]

scatter(x1,x2)

points1 = filter(p -> p ∈ reg1, points)
J1 = length(points1)  

S = DiscreteStructureFactor(points, sfactor.k; inner_box = reg1)
dS = S.S - sfactor.S # add quadradture weight w here if wanted

function objective(x)
    points = Iterators.partition(x,Dim) |> collect
    S = DiscreteStructureFactor(points, sfactor.k; inner_box = reg1)

    return sum(abs2.(S.S - sfactor.S))
end

ks = abs.(sfactor.k)


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

scale = 2.0
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

