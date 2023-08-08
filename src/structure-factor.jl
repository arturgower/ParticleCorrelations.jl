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

structure_factor(pair::DiscretePairCorrelation; kws...) = DiscreteStructureFactor(pair; kws...)

function DiscreteStructureFactor(pair::DiscretePairCorrelation{2}; dk::AbstractFloat = 0.0, maxk::AbstractFloat = 0.0)

    rs = pair.r
    drs = rs[2:end] - circshift(rs,1)[2:end]
    dr = mean(drs)

    # estimate maxk from max dr 
    if maxk == 0.0
        maxk = 2π / dr 
    end    

    # estimate dk 
    if dk == 0.0
       rmax = min(dr * length(drs), rs[end])
       dk = 2π / rmax 
    end
    
    # structure factor not correct for k = 0
    ks = dk:dk:maxk
    ks = ks[1:end-1] # the last term is the same as k = 0 when using exp(iu * dr * k)

    σs = trapezoidal_scheme(rs)

    M = [besselj(0,k * r) for k in ks, r in rs]

    f = σs .* (pair.g .- 1.0) .* rs 
    S = M * f
    S = 1.0 .+ (2π * pair.number_density) .* S

    # Need to add to the region of the integral where the pair correlation is zero
    if rs[1] >= dr 
        rs = 0.0:dr:rs[1]
        σs = trapezoidal_scheme(rs)

        M = [besselj(0, k * r) for k in ks, r in rs]

        w = - rs .* σs  
        s0 = M * w
        s0 = 1 .+ (2π * pair.number_density) .* s0

        s = s + s0
    end   

    return DiscreteStructureFactor(2, ks, S; 
        number_density = pair.number_density, 
        minimal_distance = pair.minimal_distance
    )
end

function DiscreteStructureFactor(pair::DiscretePairCorrelation{3}; dk::AbstractFloat = 0.0, maxk::AbstractFloat = 0.0)

    rs = pair.r
    drs = rs[2:end] - circshift(rs,1)[2:end]
    dr = mean(drs)

    # estimate maxk from max dr 
    if maxk == 0.0
        maxk = 2π / dr 
    end    

    # estimate dk 
    if dk == 0.0
       rmax = min(dr * length(drs), rs[end])
       dk = 2π / rmax 
    end
    
    # structure factor not correct for k = 0
    ks = dk:dk:maxk
    ks = ks[1:end-1] # the last term is the same as k = 0 when using exp(iu * dr * k)

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