module ParticleCorrelations

export Specie, Species, volume_fraction, exclusion_distance
export HardMedium
export number_density, translate_pair_correlation
export periodic_particles, trapezoidal_scheme, optimise_particulate

## Pair correlation
export DiscretePairCorrelation, DiscreteStructureFactor
export PairCorrelationType, PairCorrelation
export PercusYevick, MonteCarloPairCorrelation, HoleCorrection

export pair_correlation, structure_factor

export hole_correction_pair_correlation, gls_pair_radial_fun, pair_radial_fun
export calculate_pair_correlation, smooth_pair_corr_distance, pair_radial_to_pair_corr

using Statistics
using SpecialFunctions
using LinearAlgebra
using Reexport
using Optim

@reexport using MultipleScattering
# using MultipleScattering: PhysicalMedium, AbstractParticle, Particle, Shape, Box, Sphere, Circle, outer_radius, random_particles

import MultipleScattering: outer_radius, volume, random_particles

using HCubature: hcubature, hquadrature
# using ClassicalOrthogonalPolynomials: Legendre
using LegendrePolynomials

Sphere = MultipleScattering.Sphere

include("particles.jl")
include("pair-correlation.jl")
include("structure-factor.jl")
include("numerical.jl")
include("inverse-problem.jl")


end
