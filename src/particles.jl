"""
    HardMedium

Represents a solid type of material. Solid particles can not overlap. Currently the package has no soft medium, but this can easily be added.
"""
struct HardMedium{Dim} <: PhysicalMedium{Dim,1} end

name(e::HardMedium{Dim}) where {Dim} = "$(Dim)D HardMedium"


""" Specie

Represents a set of particles which are all the same. The type of particle is given by `Specie.particle` and the volume fraction this specie occupies is given by `Specie.volume_fraction`.

We can use `Specie.numberofparticles` to specify the number of particles, otherwise for an infinite `Specie.numberofparticles = Inf`.

The minimum distance between any two particles will equal `outer_radius(Specie) * Specie.separation_ratio`.
"""
struct Specie{Dim,P<:AbstractParticle{Dim}}
    particle::P
    volume_fraction::Float64
    separation_ratio::Float64
end

# Convenience constructor which does not require explicit types/parameters
function Specie(p::AbstractParticle{Dim};
        number_density::AbstractFloat = 0.0,
        volume_fraction::AbstractFloat = number_density * volume(p),
        separation_ratio::AbstractFloat = 1.0,
        exclusion_distance::AbstractFloat = separation_ratio
    ) where Dim

    if number_density == 0.0 && volume_fraction == 0.0
            @warn println("zero volume fraction or number density was chosen.")
    end

    Specie{Dim,typeof(p)}(p, volume_fraction, exclusion_distance)
end

function Specie(medium::P,s::S; kws...) where {Dim,P<:PhysicalMedium{Dim},S<:Shape{Dim}}
    Specie(Particle(medium, s); kws...)
end

function Specie(medium::P, radius::AbstractFloat; kws...) where P<:PhysicalMedium
    Specie(Particle(medium, radius); kws...)
end

# Shorthand for all Vectors of species
Species{Dim,P} = Vector{S} where S<:Specie{Dim,P}


"Returns the volume fraction of the specie."
volume_fraction(s::Specie) = s.volume_fraction
volume_fraction(ss::Species) = sum(volume_fraction.(ss))

volume(s::Specie) = volume(s.particle)


"""
    number_density(s::Specie)

Gives the number of particles per unit volume. Note this is given exactly by `N / V` where `V` is the volume of the region containing the origins of all particles. For consistency, [`volume_fraction`](@ref) is given by `N * volume(s) / V`.
"""
number_density(s::Specie) = s.volume_fraction / volume(s)
# number_density(s::Specie{2}) where {T} = s.volume_fraction / (outer_radius(s.particle)^2 * pi)
# number_density(s::Specie{3}) where {T} = s.volume_fraction / (T(4/3) * outer_radius(s.particle)^3 * pi)
number_density(ss::Species) = sum(number_density.(ss))

outer_radius(s::Specie) = outer_radius(s.particle)
exclusion_distance(s::Specie) = outer_radius(s) * s.separation_ratio