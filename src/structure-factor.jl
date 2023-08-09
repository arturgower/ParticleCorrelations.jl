"""
    DiscreteStructureFactor

Represents the structure factor between two types of particles. The particles could be the same or different species.
"""
struct DiscreteStructureFactor{Dim}
    "the wavenumber"
    k::Vector{Float64}
    "the vector of the radial structure factor"
    S::Vector{Float64}
    "the minimal distance between the two species"
    minimal_distance::Float64
    "the average number of particles divided by the volume containing the centre of the particles"
    number_density::Float64

    function DiscreteStructureFactor{Dim}(k::AbstractVector, S::AbstractVector, minimal_distance::Float64, number_density::Float64;
            tol::AbstractFloat = 1e-3 
        ) where Dim
        if !isempty(S) && size(S) != size(k)
            @error "the size of vector of wavenumbers `k` (currently $(size(k))) should be the same as the size of the structure factor vector `S` (currently $(size(S)))."
        end
        new{Dim}(k,S,minimal_distance,number_density)
    end
end

function DiscreteStructureFactor(Dim::Int, k::AbstractVector, S::AbstractVector;
    number_density::AbstractFloat = 0.0,
    minimal_distance::AbstractFloat = 0.0
) 
    DiscreteStructureFactor{Dim}(k,S,minimal_distance,number_density)
end

structure_factor(pair::DiscretePairCorrelation, ks::AbstractVector = Float64[]; kws...) = DiscreteStructureFactor(pair, ks; kws...)

"""
    structure_factor(particles::Vector{AbstractParticle}, distances::AbstractVector)

Calculates the isotropic structure_factor from one configuration of particles.
"""

structure_factor(particles::Vector{p} where p <: AbstractParticle, wavenumbers::AbstractVector = Float64[]; kws...) = DiscreteStructureFactor(particles, wavenumbers; kws...) 

"""
structure_factor(particle_centres::Vector, distances::AbstractVector)

Calculates the isotropic structure_factor from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.
"""
structure_factor(particle_centres::Vector{v} where v <: AbstractVector, wavenumbers::AbstractVector  = Float64[]; kws...) = DiscreteStructureFactor(particle_centres, wavenumbers; kws...)


function DiscreteStructureFactor(particles::Vector{p} where p <: AbstractParticle{Dim}, ks::AbstractVector{T}; kws...
    ) where {T, Dim}
    
    particle_centres = [origin(p) for p in particles]

    return DiscreteStructureFactor(particle_centres, ks; kws...)
end

function DiscreteStructureFactor(points::Vector{v} where v <: AbstractVector{T}, ks::AbstractVector{T} = Float64[]; correlation_length::T = 0.0) where T <: AbstractFloat
    
    box = Box(points)
    Dim = length(box.dimensions)  
    
    dims = box.dimensions .- 2 * correlation_length
    inner_box = Box(box.origin, dims)

    inner_points = filter(p -> p ∈ inner_box, points) 

    Rijs = [ 
        norm(p1 - p2)        
    for p1 in inner_points, p2 in points][:]
    
    filter!(Rij -> Rij > 0, Rijs)

    J1 = length(inner_points)

    dS = if Dim == 3
        [sum(sin.(k .* Rijs) ./ (k .* Rijs)) for k in ks] ./ J1
    elseif Dim == 2    
        [sum(besselj0.(abs(k) .* Rijs)) for k in ks] ./ J1
    else @error "not implemented for Dim = $Dim"    
    end  

    number_density = J1 / volume(inner_box)

    return DiscreteStructureFactor(Dim, ks, T(1) .+ dS;
        number_density = number_density,
        minimal_distance = minimum(Rijs)    
    ) 
end   

function DiscreteStructureFactor(pair::DiscretePairCorrelation{2}, ks::AbstractVector = Float64[])

    rs = pair.r
    drs = rs[2:end] - circshift(rs,1)[2:end]
    dr = mean(drs)

    if isempty(ks)
        # estimate maxk from max dr 
        maxk = 2π / dr 

        # estimate dk 
        rmax = min(dr * length(drs), rs[end])
        dk = 2π / rmax 
    
        # structure factor not correct for k = 0
        ks = dk:dk:maxk
        ks = ks[1:end-1] # the last term is the same as k = 0 when using exp(iu * dr * k)
    end    

    σs = trapezoidal_scheme(rs)

    M = [besselj(0,k * r) for k in ks, r in rs]

    f = σs .* (pair.g .- 1.0) .* rs 
    s = M * f
    s = 1.0 .+ (2π * pair.number_density) .* s

    # Need to add to the region of the integral where the pair correlation is zero
    if rs[1] >= dr 
        rs = 0.0:dr:rs[1]
        σs = trapezoidal_scheme(rs)

        M = [besselj(0, k * r) for k in ks, r in rs]

        w = - rs .* σs  
        s0 = M * w
        s0 = (2π * pair.number_density) .* s0

        s = s + s0
    end   

    return DiscreteStructureFactor(2, ks, s; 
        number_density = pair.number_density, 
        minimal_distance = pair.minimal_distance
    )
end

function DiscreteStructureFactor(pair::DiscretePairCorrelation{3}, ks::AbstractVector = Float64[])

    rs = pair.r
    drs = rs[2:end] - circshift(rs,1)[2:end]
    dr = mean(drs)

    if isempty(ks)
        # estimate maxk from max dr 
        maxk = 2π / dr 

        # estimate dk 
        rmax = min(dr * length(drs), rs[end])
        dk = 2π / rmax 
    
        # structure factor not correct for k = 0
        ks = dk:dk:maxk
        ks = ks[1:end-1] # the last term is the same as k = 0 when using exp(iu * dr * k)
    end    
    
    σs = trapezoidal_scheme(rs)

    M = [sin(k * r) / (k * r) for k in ks, r in rs]

    w = (pair.g .- 1.0) .* rs .^2 .* σs  
    s = M * w
    s = 1 .+ (4π * pair.number_density) .* s

    # Need to add to the integral the pair of the pair correlation that is missing
    if rs[1] >= dr 
        rs = 0.0:dr:rs[1]
        σs = trapezoidal_scheme(rs)

        M = [sin(k * r) / (k) for k in ks, r in rs]

        w = - rs .* σs  
        s0 = M * w
        s0 = (4π * pair.number_density) .* s0

        s = s + s0
    end    

    return DiscreteStructureFactor(3, ks, s; 
        number_density = pair.number_density, 
        minimal_distance = pair.minimal_distance
    )
end