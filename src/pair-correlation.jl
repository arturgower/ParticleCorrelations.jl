## A the top of the type tree we have:

"""
    PairCorrelationType

A type used to specify what type of pair correlation is to be used. This is like a tag, or a option, to specify which pair correlation is wanted.
"""
abstract type PairCorrelationType end

"""
    PairCorrelation

A type used to store a pair-correlation. This represents the calculated pair-correlation.
"""
abstract type PairCorrelation end


struct PercusYevick{Dim} <: PairCorrelationType
    "The size of each element of the mesh relative to the particle radius"
    meshsize::Float64
    "Relative tolerance used when calculating the Percus-Yevick approximation"
    rtol::Float64
    "Maximum number of quadture evaluations when calculating the Percus-Yevick approximation"
    maxevals::Int
    "Maximum number of points for the pair correlation"
    maxlength::Int
end

PercusYevick(Dim; rtol::AbstractFloat = 1e-3, meshsize::AbstractFloat = 0.2, maxevals::Int = Int(2e4), maxlength::Int = 50) = PercusYevick{Dim}(meshsize, rtol, maxevals, maxlength)

"""
    MonteCarloPairCorrelation{Dim} <: PairCorrelationType

Currently only used to create pair-correlations for particles that are uniformly randomly placed, except they can not overlap.
"""
struct MonteCarloPairCorrelation{Dim} <: PairCorrelationType
    "The max nubmer elements for the mesh"
    maxlength::Int
    "The size of each element of the mesh relative to the particle radius"
    meshsize::Float64
    "Avoid sampling particles that are boundary_distance * particle radius close to the boundary. This is because particles are more concentrated near the boundary. "
    boundary_distance::Float64
    "Number of partilces configurations to take into account"
    iterations::Int
    "Maximum number of points for the pair correlation"
    numberofparticles::Int
    "Relative tolerance"
    rtol::Float64
end

function MonteCarloPairCorrelation(Dim;
        iterations::Int = 1, meshsize::AbstractFloat = 0.2,
        numberofparticles::Number = 1e4,
        boundary_distance::AbstractFloat = 2.0,
        maxlength::Int = 50,
        rtol::AbstractFloat = 1e-3
    )
    MonteCarloPairCorrelation{Dim}(maxlength, meshsize, boundary_distance, iterations, Int(round(numberofparticles)), rtol)
end

struct HoleCorrection <: PairCorrelationType
end

"""
    DiscretePairCorrelation

Represents the pair correlation between two types of particles. The particles could be the same or different species. 

This struct has the field `r`, which represent the radial distance, and the field `g` is the value of the pair correlation. For radial distances `< r[1]` we assume that the pair correlation is zero. For radial distances `> r[end]` we assume the pair correlation is 1. 
"""
struct DiscretePairCorrelation{Dim} <: PairCorrelation
    "distance between particles centres"
    r::Vector{Float64}
    "the value of the pair correlation"
    g::Vector{Float64}
    "Minimal distance between particle centres'"
    minimal_distance::Float64
    "the average number of particles divided by the volume containing the centre of the particles"
    number_density::Float64

    function DiscretePairCorrelation{Dim}(r::AbstractVector, g::AbstractVector, minimal_distance::Float64, number_density::Float64 = -1.0;
            tol::AbstractFloat = 1e-3, rescale = false 
        ) where Dim
        if !isempty(g) && size(g) != size(r)
            @error "the size of vector of distances `r` (currently $(size(r))) should be the same as the size of the pair-correlation variation `g` (currently $(size(g)))."
        end
        if !isempty(g) && abs(g[end]) > tol
            @warn "For the pair-correlation to be accurate, we expect it to be long enough (in terms of the distance `r`) such that the particle positions become uncorrelatd. They become uncorrelated when `g[end]` tends to zero."
        end
        
        # Predict the number density from the restriction that is due to the link with a probability distribution, and assuming homogeneous distribution of particles. 
        if !isempty(r)
            σs = trapezoidal_scheme(r)

            sumg = sum(g .* r .^ (Dim - 1) .* σs)

            number_density_predict = 1 / (2*(Dim - 1) * π * (r[end] ^ Dim / Dim - sumg))

            if number_density == -1.0
                number_density = number_density_predict
            end

            number_density_error = abs(number_density / number_density_predict - 1)
        
            # This formula seems to sensitive to numerical errors.
            # if number_density_error > tol
            #     println("The number density that was specified was $(number_density), whereas the calculated number density was $(number_density_predict), a reltive error of $(number_density_error)")
            # end
        
            if rescale
                println("rescaling this pair correlation according to the number density provided.")
                # this rescaling is due to a restriction imposed by the connection to probability distributions

                b = r[end] ^ Dim / Dim - 1 / (2*(Dim - 1) * π * number_density)
                α = b / sumg

                if α < 0
                    @error "rescaling has failed and has given a negative scale for the pair-correlation"
                end     

                g = α .* g
            end
        end

        new{Dim}(r,g,minimal_distance,number_density)
    end
end

function number_density(discrete_pair::DiscretePairCorrelation{Dim}) where Dim

    σs = trapezoidal_scheme(discrete_pair.r)

    sumg = sum(discrete_pair.g .* discrete_pair.r .^ (Dim - 1) .* σs)

    term = discrete_pair.r[end] ^ Dim / Dim - sumg;

    number_density_predict = 1 / (2*(Dim - 1) * π * term)
end

function translate_pair_correlation(discrete_pair::DiscretePairCorrelation{Dim}, number_density::Float64; dg::Function = r -> exp(- 6r / discrete_pair.r[end])) where Dim

    R =  discrete_pair.r[end]
    σs = trapezoidal_scheme(discrete_pair.r)

    dps = dg.(discrete_pair.r)

    sumg = sum(discrete_pair.g .* discrete_pair.r .^ (Dim - 1) .* σs)
    sumdp = sum(dps .* discrete_pair.r .^ (Dim - 1) .* σs)

    α = (R^Dim / Dim - (1 / (2*(Dim - 1) * pi * number_density)) - sumg) / sumdp

    scaled_discrete_pair = DiscretePairCorrelation(Dim, discrete_pair.r, discrete_pair.g + α .* dps; number_density = number_density)

    return scaled_discrete_pair
end

function DiscretePairCorrelation(Dim::Int,r::AbstractVector, g::AbstractVector;
    number_density::AbstractFloat = -1.0,
    minimal_distance::AbstractFloat = r[1],
    tol::AbstractFloat = 1e-3
) 
    DiscretePairCorrelation{Dim}(r,g,minimal_distance,number_density; tol = tol)
end



# Have no pair correlation, then only hole correction will be used.
function DiscretePairCorrelation(s1::Specie{Dim}, s2::Specie{Dim}) where Dim
    T = typeof(outer_radius(s1))
    
    return DiscretePairCorrelation(Dim,T[],T[]; minimal_distance = s1.separation_ratio * outer_radius(s1) + s2.separation_ratio * outer_radius(s2))
end

pair_correlation(Dim::Int, r::AbstractVector{T},g::AbstractVector{T}) where T <: AbstractFloat = DiscretePairCorrelation(Dim,r,g)

"""
    pair_correlation(s::Specie, pairtype::PairCorrelationType)

Generates a DiscretePairCorrelation for the specie `s` by using the type of paircorrelation pairtype provided. The distances where the pair correlation is sampled is calculated from the properties of pairtype.
"""
pair_correlation(s::Specie, pairtype::PT; kws...) where PT <: PairCorrelationType = DiscretePairCorrelation(s, pairtype; kws...)

"""
    pair_correlation(s::Specie, pairtype::PairCorrelationType, distances::AbstractVector)

Generates a DiscretePairCorrelation for the specie `s` by using the type of paircorrelation pairtype provided, and for the radial distances provided.
"""
pair_correlation(s::Specie, pairtype::PairCorrelationType, distances::AbstractVector) = DiscretePairCorrelation(s, pairtype, distances)

"""
    pair_correlation(s1::Specie, s2::Specie, pairtype::PairCorrelationType; kws...)

currently provides an approximation for the pair-correlation for two different types of particles `s1` and `s2`.    
"""
pair_correlation(s1::Specie, s2::Specie, pairtype::PairCorrelationType; kws...) = DiscretePairCorrelation(s1, s2, pairtype; kws...)


"""
    pair_correlation(particles::Vector{AbstractParticle}, distances::AbstractVector)

Calculates the isotropic pair correlation from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.
"""
pair_correlation(particles::Vector{p} where p <: AbstractParticle, distances::AbstractVector; kws...) = DiscretePairCorrelation(particles, distances; kws...) 

"""
    pair_correlation(particle_centres::Vector, distances::AbstractVector)

Calculates the isotropic pair correlation from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.
"""
pair_correlation(particle_centres::Vector{v} where v <: AbstractVector, distances::AbstractVector; kws...) = DiscretePairCorrelation(particle_centres, distances; kws...)


function DiscretePairCorrelation(s::Specie{Dim}, pairtype::PT; 
        distances::AbstractVector{T} = Float64[]
    ) where {Dim,  T<:AbstractFloat, PT <: PairCorrelationType}

    r1 = exclusion_distance(s)

    R = 2r1
    numdensity = number_density(s)

    automatic_dist = if isempty(distances)
        dr = pairtype.meshsize * r1

        # For MonteCarlo calculations the mesh points are at the centre of the bins.
        d1 = if PT <: MonteCarloPairCorrelation
            R+dr/2
        else R
        end

        maximum_distance = d1 .+ dr .* pairtype.maxlength
        
        # the distances should be in the centre of the mesh element.
        distances = d1:dr:maximum_distance

        true
    else false
    end

    d = DiscretePairCorrelation(s, pairtype, distances);
    g = d.g

    if automatic_dist
        i = findfirst(reverse(abs.(d.g .- 1)) .> pairtype.rtol)
        if isnothing(i)
            g = typeof(d.g)[]
            distances = typeof(d.g)[]
        elseif i > 1
            g = d.g[1:end-i+2]
            distances = distances[1:end-i+2]
        end
    end

    return DiscretePairCorrelation(Dim, distances, g; number_density = d.number_density)
end


function DiscretePairCorrelation(s1::Specie{Dim}, s2::Specie{Dim}, pairtype::PairCorrelationType; kws...) where Dim

    if outer_radius(s1) != outer_radius(s2)
        @warn "Calculating any pair-correlation has only been implemented for particles of the same size. Will use a crude approximation which combines both particles into one average particle"
    end

    a = (outer_radius(s1) + outer_radius(s2)) / 2
    sep_ratio = (s1.separation_ratio + s2.separation_ratio) / 2
    numdensity = number_density(s1) + number_density(s2)


    # replace both particles by an average particle
    sm = Specie(
        HardMedium{Dim}(), radius1;
        number_density = numdensity,
        separation_ratio = sep_ratio
    );

    return DiscretePairCorrelation(sm, pairtype; kws...)
end

function DiscretePairCorrelation(s::Specie{3, P} where P<:AbstractParticle{3}, pairtype::PercusYevick{3}, distances::AbstractVector{T}) where T

    R = 2 * outer_radius(s) * s.separation_ratio
    numdensity = number_density(s)

    rtol = pairtype.rtol;
    maxevals = pairtype.maxevals;

    # Note that f is the volume fraction of the species including the exclusion volume around each particle
    f = numdensity * π * R^3 / 6
    α = - (1 + 2f)^2 / (1 - f)^4
    β = 6f * (1 + f/2)^2 / (1 - f)^4
    δ = - f * (1 + 2f)^2 / (2 * (1 - f)^4)

    function F(x)
        A = 24δ / x^6 - 2β / x^4
        B = (α + 2β + 4δ) / x^3 - 24δ / x^5
        C = - (α + β + δ) / x^2 + (2β + 12δ) / x^4 - 24δ / x^6

        return A + B * sin(x) + C * cos(x)
    end

    G(x) = (α + 2β + 4δ) * sin(x) / x^2 - (α + β + δ) * cos(x) / x

    ker_fun(r) = x -> (x * F(x) / (1 - 24f * F(x)) - G(x)) * sin(r * x / R)

    ## Below we estimate the domain of integration needed

    ker = ker_fun(R)

    maxker = abs.(ker.(distances ./ R)) |> maximum

    d = distances[end] / R
    max_xs = LinRange(d, 20 * d, 200)

    imax = findfirst(abs.(ker.(max_xs) ./ maxker) .< rtol^2)
    max_x = isnothing(imax) ? max_xs[end] : max_xs[imax]

    (I,E) = hquadrature(ker, eps(T), max_x;
        rtol = rtol, maxevals = maxevals
    )
    if abs(E/I) > rtol
        @warn "For rtol = $rtol, the estimated relative error when calculating the integral of the Percus-Yevick is: ", abs(E/I)
    end

    dp = map(distances) do r
        if r < R return zero(T) end

        ker = ker_fun(r)
        int_ker = hquadrature(ker, 10eps(T), max_x;
            rtol = rtol, maxevals = maxevals
        )[1]

        η = r / R
        9f * (1 + f) / (2 * η * (1 - f)^3) + 2 / (η*pi) * int_ker
    end

    g = dp .+ 1.0;
    g[distances .< R] .= zero(T)

    return DiscretePairCorrelation(3, distances, g; number_density = numdensity, tol = pairtype.rtol)
end

function DiscretePairCorrelation(s::Specie{Dim}, pairtype::MonteCarloPairCorrelation{Dim}, distances::AbstractVector{T}) where {T, Dim}

    numdensity = number_density(s)
    a = outer_radius(s)
    R = a * s.separation_ratio * pairtype.boundary_distance

    vol = pairtype.numberofparticles / numdensity
    l = (vol)^(1/Dim)

    # NOTE: when calling random_particles, the region specified will completely contain the whole of all particles. However, to better align with the theory, we need to use the region containing the particle centres, which leeds to the correction below.
    la = l + 2a

    zs = zeros(T,Dim)
    region_shape = Box(zs .+ la)
    region_shape_numdensity = Box(zs .+ l)

    large_region_shape = Box(zs .+ la  .+ 2R)
    large_region_shape_numdensity = Box(zs .+ l .+ 2R)

    large_N = Int(round(numdensity * volume(large_region_shape_numdensity)))

    dpcs = map(1:pairtype.iterations) do i
        ps = random_particles(s.particle.medium, s.particle.shape, large_region_shape, large_N;
            separation_ratio = s.separation_ratio
        );

        # using the cookie cutter method to keep only the particles slightly away from the boundary where particles concentrate.
        # ps = filter(p -> p ⊆ region_shape, ps);
        particle_centres = origin.(filter(p -> p ⊆ region_shape, ps));

        DiscretePairCorrelation(particle_centres, distances;
            region_particle_centres = region_shape_numdensity
        )
    end

    nums = [d.number_density for d in dpcs]
    achieved_number_density = mean(nums)
    std_number_density = std(nums)

    if abs(achieved_number_density / number_density(s) - 1) > 0.01
        @warn "The requested volume fraction of the pair correlation was $(s.volume_fraction). The achieved volume fraction was $(achieved_number_density * volume(s)) with std $( std(nums) *  volume(s))"
    end

    return DiscretePairCorrelation(Dim, distances, mean(d.g for d in dpcs); number_density = achieved_number_density)

end

function hole_correction_pair_correlation(x1::AbstractVector{T},s1::Specie, x2::AbstractVector{T},s2::Specie) where T <: Number
    overlapping = norm(x1 - x2) > outer_radius(s1) * s1.separation_ratio + outer_radius(s2) * s2.separation_ratio

    return  overlapping ? one(T) : zero(T)
end


function DiscretePairCorrelation(particles::Vector{p} where p <: AbstractParticle, distances::AbstractVector{T}; kws...
) where {T}

    particle_centres = [origin(p) for p in particles]

    return DiscretePairCorrelation(particle_centres, distances; kws...)
end

function DiscretePairCorrelation(particle_centres::Vector{v} where v <: AbstractVector{T}, distances::AbstractVector{T};
        dz::T = distances[2] - distances[1], 
        maximum_distance::T = distances[end] + dz/2, # the max possible distance between any two particles
        minimum_distance::T = distances[1] - dz/2, # the min possible distance between any two particles. 
        region_particle_centres::Shape = Box(zeros(length(particle_centres[1])))
    ) where {T}

    Dim = length(particle_centres[1])
    p2s = particle_centres

    if norm(region_particle_centres.dimensions) == 0.0
        ind = CartesianIndices(p2s[1])
        xs = [p[i] for p in p2s, i in ind]

        xmin = minimum(xs; dims = 1)
        xmax = maximum(xs; dims = 1)

        c = (xmin + xmax)[:] ./ 2.0;
        dims = (xmax - xmin)[:];
        outer_box = Box(c,dims)

    else
        dims = region_particle_centres.dimensions
        c = origin(region_particle_centres)
        outer_box = region_particle_centres
    end
    inner_box = Box(c,dims .- 2 * maximum_distance)

    p1s = filter(x -> x ∈ inner_box, p2s)

    if length(p1s) < 10
        @error "There are only $(length(p1s)) particles in the feasible region. This is not enough to calculate the pair-correlation. To increase this number, and get a more accurate result, try: 1) increasing the number of particles or 2) using a shorter distance for the pair-correlation"
    end

    N = length(distances)
    bins = zeros(N);

    for p1 in p1s, p2 in p2s
        dist = norm(p1 - p2)

        if minimum_distance < dist < maximum_distance
            n = 1 + Int(round(
                -1/2 + N * (dist - minimum_distance) / (maximum_distance - minimum_distance)
            ))
            bins[n] += 1.0
        end
    end

    J1 = length(p1s)
    J2 = length(p2s)

    numdensity = J2 / volume(outer_box)
    # numdensity1 = J1 / volume(inner_box)

    # scaling = (1 / ((2 * (Dim - 1)) * pi * dz)) * J2 / ((J2 - 1) * J1 * numdensity)
    # scaling = (1 / ((2 * (Dim - 1)) * pi * dz)) / (J1 * numdensity)
    scaling = (3 / ((Dim + 1) * pi)) / (J1 * numdensity)

    # g = scaling .* bins ./ (distances .^(Dim - 1))
    g = scaling .* bins ./ ((distances .+ dz/2) .^ Dim - (distances .- dz/2) .^ Dim) 

    return DiscretePairCorrelation(Dim,distances,g; number_density = numdensity)

end

function smooth_pair_corr_distance(pair_corr_distance::Function, a12::T; smoothing::T = T(0), max_distance::T = T(20*a12),
        polynomial_order::Int = 15, mesh_size::Int = 10*polynomial_order + 4
    ) where T

    if smoothing > T(1.0)
        @warn "smoothing should be in [0,1], setting smoothing = 1"
        smoothing = T(1.0)
    elseif smoothing < T(0.0)
        @warn "smoothing should be in [0,1], setting smoothing = 0"
        smoothing < T(0.0)
    end

    zs = collect(LinRange(zero(T), max_distance, mesh_size))

    inds = findall(abs.(zs .- a12) .< a12*smoothing/T(2) )
    deleteat!(zs,inds)
    data = pair_corr_distance.(zs)

    # P = Legendre()
    # ls = 0:polynomial_order |> collect

    # Pmat = P[2zs ./ max_distance .- T(1.0), ls .+ 1];

    Pmats = [collectPl(2z / max_distance - T(1.0), lmax = polynomial_order) for z in zs];
    Pmat = hcat(Pmats...) |> transpose;

    # projector_mat = inv(transpose(Pmat) * (Pmat)) * transpose(Pmat);
    # pls = projector_mat * data
    # Pmat * pls ~ data
    pls = Pmat \ data

    return function (z)
        # Ps = P[2z / max_distance - T(1.0), ls .+ 1]
        Ps = collectPl(2z / max_distance - T(1.0), lmax = polynomial_order) |> collect
        return sum(Ps .* pls)
    end
end

# function constant_number_density_fun(material::Material)
#
#     @warn "the constant_number_density doesn't currently verify if the whole particle is inside the shape"
#
#     prob = function (x1::AbstractVector, s1::Specie)
#         if x1 ∈ sh
#             number_density(s1)
#         else
#             zero(typeof(V))
#         end
#     end
#
#     return prob
# end
"""
    gls_pair_radial_fun(pair_corr_distance::Function, a12::T; polynomial_order::Int, mesh_size::Int)

Return a function `gls_fun`. For any radial distances `r1` and `r2` we have `gls = gls_fun(r1,r2)` where `gls` is an array of the Legendre coefficients for the pair correlation.
    
Using mathematics, we have that such that ``g(r_1,r_2,\\cos \\theta_{12}) = \\sum_{\\ell_1 =0} \\frac{2\\ell_1 + 1}{4\\pi} gls[\\ells+1] P_{\\ell_1}(\\cos \\theta_{12}) ``, where ``g(r_1,r_2,\\cos \\theta_{12})`` is the radially symmetric pair-correlation, so it depends only on the radial distances ``r_1`` and ``r_2``, and the angle between two position vectors ``\\theta_{12}``.

The function `gls_fun` is calculated from the function `pair_corr_distance`, where `pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12)))` gives the pair correlation.
"""
function gls_pair_radial_fun(pair_corr_distance::Union{Function,AbstractArray}, a12::T;
            polynomial_order::Int = 15,
            mesh_size::Int = 10*polynomial_order + 4,
            sigma_approximation = true
        ) where T

    # P = Legendre()
    ls = 0:polynomial_order |> collect

    if sigma_approximation
        sigmas = [one(T); sin.(pi .* ls[2:end] ./ (polynomial_order+1)) ./ (pi .* ls[2:end] ./ (polynomial_order+1))]
    else
        sigmas = ones(T,polynomial_order+1)
    end
    S = diagm(sigmas)

    us = LinRange(-1.0, 1.0, mesh_size)

    # Pmat = P[us, ls .+ 1];
    Pmats = [collectPl(u, lmax = polynomial_order) for u in us];
    Pmat = hcat(Pmats...) |> transpose;

    Pmat = [Pmat[i] * T(2i[2] - 1) / (4pi) for i in CartesianIndices(Pmat)];

    # Pmat = S * Pmat;

    # projector_mat =  S * inv(transpose(Pmat) * (Pmat)) * transpose(Pmat);

    return function (r1,r2)
        if (r1 + r2 < a12)
            return zeros(typeof(a12), polynomial_order + 1)
        else
            data = pair_corr_distance.(sqrt.(r1^2 .+ r2^2 .- 2r1 .* r2 .* us))
            # pls = projector_mat * data
            # data ~ Pmat * pls
            pls = Pmat \ data

            return S * pls
        end
    end
end

"""
    pair_radial_fun(pair_corr::Function, a12::T; polynomial_order::Int, mesh_size::Int)

Return a function `pair_radial` such that `pair_radial(r1,r2, cos(θ12))` gives the pair correlation particles at the radial distances r1 and r2, with and angle of `θ12` between them.

The function `pair_radial` is calculated from the function `pair_corr_distance`, where `pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12)))` gives the pair correlation.
"""
function pair_radial_fun(pair_corr_distance::Function, a12::T; polynomial_order::Int = 15, kws...) where T

    gls_fun = gls_pair_radial_fun(pair_corr_distance, a12; polynomial_order = polynomial_order, kws...)
    # P = Legendre()

    return function (r1,r2,u)
        # Pus = P[u, 1:(polynomial_order + 1) |> collect] .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)

        Pus = collectPl(u, lmax = polynomial_order) |> collect
        Pus = Pus .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)

        return sum(Pus .* gls_fun(r1,r2))
    end
end

function pair_radial_to_pair_corr(pair_radial::Function)
    function (x1,s1,x2,s2)
        if norm(x1) < 1e-12 || norm(x2) < 1e-12
            pair_radial(norm(x1),norm(x2),0.0)
        else
            pair_radial(norm(x1),norm(x2),dot(x1,x2) / (norm(x1)*norm(x2)))
        end
    end
end
