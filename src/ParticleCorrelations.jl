module ParticleCorrelations

export Specie, Species, number_density, volume_fraction, exclusion_distance
export HardMedium

## Pair correlation
export PairCorrelationType, PairCorrelation
export PercusYevick, MonteCarloPairCorrelation, HoleCorrection, DiscretePairCorrelation
export hole_correction_pair_correlation, gls_pair_radial_fun, pair_radial_fun
export calculate_pair_correlation, smooth_pair_corr_distance, pair_radial_to_pair_corr

using Statistics
using LinearAlgebra
using Reexport

@reexport using MultipleScattering
# using MultipleScattering: PhysicalMedium, AbstractParticle, Particle, Shape, Box, Sphere, Circle, outer_radius, random_particles

import MultipleScattering: outer_radius, volume, random_particles

using HCubature: hcubature, hquadrature
using ClassicalOrthogonalPolynomials: Legendre

include("particles.jl")
include("pair-correlations.jl")

end
