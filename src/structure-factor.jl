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
    DiscreteStructureFactor{Dim}(r,g,minimal_distance,number_density)
end

# function DiscreteStructureFactor(pair::DiscretePairCorrelation{2})

#     DiscreteStructureFactor{2}(r,g,minimal_distance,number_density)
# end