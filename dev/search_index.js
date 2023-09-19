var documenterSearchIndex = {"docs":
[{"location":"library/library/#Base","page":"Library","title":"Base","text":"","category":"section"},{"location":"library/library/","page":"Library","title":"Library","text":"CurrentModule = ParticleCorrelations","category":"page"},{"location":"library/library/#All-types","page":"Library","title":"All types","text":"","category":"section"},{"location":"library/library/","page":"Library","title":"Library","text":"Modules = [ParticleCorrelations]\nOrder   = [:function, :constant, :type]","category":"page"},{"location":"library/library/#ParticleCorrelations.gls_pair_radial_fun-Union{Tuple{T}, Tuple{Union{Function, AbstractArray}, T}} where T","page":"Library","title":"ParticleCorrelations.gls_pair_radial_fun","text":"gls_pair_radial_fun(pair_corr_distance::Function, a12::T; polynomial_order::Int, mesh_size::Int)\n\nReturn a function gls_fun. For any radial distances r1 and r2 we have gls = gls_fun(r1,r2) where gls is an array of the Legendre coefficients for the pair correlation.\n\nUsing mathematics, we have that such that g(r_1r_2cos theta_12) = sum_ell_1 =0 frac2ell_1 + 14pi glsells+1 P_ell_1(cos theta_12), where g(r_1r_2cos theta_12) is the radially symmetric pair-correlation, so it depends only on the radial distances r_1 and r_2, and the angle between two position vectors theta_12.\n\nThe function gls_fun is calculated from the function pair_corr_distance, where pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12))) gives the pair correlation.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.number_density-Tuple{Specie}","page":"Library","title":"ParticleCorrelations.number_density","text":"number_density(s::Specie)\n\nGives the number of particles per unit volume. Note this is given exactly by N / V where V is the volume of the region containing the origins of all particles. For consistency, volume_fraction is given by N * volume(s) / V.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_correlation-Tuple{Specie, PairCorrelationType, AbstractVector{T} where T}","page":"Library","title":"ParticleCorrelations.pair_correlation","text":"pair_correlation(s::Specie, pairtype::PairCorrelationType, distances::AbstractVector)\n\nGenerates a DiscretePairCorrelation for the specie s by using the type of paircorrelation pairtype provided, and for the radial distances provided.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_correlation-Tuple{Specie, Specie, PairCorrelationType}","page":"Library","title":"ParticleCorrelations.pair_correlation","text":"pair_correlation(s1::Specie, s2::Specie, pairtype::PairCorrelationType; kws...)\n\ncurrently provides an approximation for the pair-correlation for two different types of particles s1 and s2.    \n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_correlation-Tuple{Vector{p} where p<:AbstractParticle, AbstractVector{T} where T}","page":"Library","title":"ParticleCorrelations.pair_correlation","text":"pair_correlation(particles::Vector{AbstractParticle}, distances::AbstractVector)\n\nCalculates the isotropic pair correlation from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_correlation-Tuple{Vector{v} where v<:(AbstractVector{T} where T), AbstractVector{T} where T}","page":"Library","title":"ParticleCorrelations.pair_correlation","text":"pair_correlation(particle_centres::Vector, distances::AbstractVector)\n\nCalculates the isotropic pair correlation from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_correlation-Union{Tuple{PT}, Tuple{Specie, PT}} where PT<:PairCorrelationType","page":"Library","title":"ParticleCorrelations.pair_correlation","text":"pair_correlation(s::Specie, pairtype::PairCorrelationType)\n\nGenerates a DiscretePairCorrelation for the specie s by using the type of paircorrelation pairtype provided. The distances where the pair correlation is sampled is calculated from the properties of pairtype.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.pair_radial_fun-Union{Tuple{T}, Tuple{Function, T}} where T","page":"Library","title":"ParticleCorrelations.pair_radial_fun","text":"pair_radial_fun(pair_corr::Function, a12::T; polynomial_order::Int, mesh_size::Int)\n\nReturn a function pair_radial such that pair_radial(r1,r2, cos(θ12)) gives the pair correlation particles at the radial distances r1 and r2, with and angle of θ12 between them.\n\nThe function pair_radial is calculated from the function pair_corr_distance, where pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12))) gives the pair correlation.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.structure_factor","page":"Library","title":"ParticleCorrelations.structure_factor","text":"structurefactor(particlecentres::Vector, distances::AbstractVector)\n\nCalculates the isotropic structure_factor from one configuration of particles. To use many configurations of particles, call this function for each, then take the average of the pair-correlation.\n\n\n\n\n\n","category":"function"},{"location":"library/library/#ParticleCorrelations.volume_fraction-Tuple{Specie}","page":"Library","title":"ParticleCorrelations.volume_fraction","text":"Returns the volume fraction of the specie.\n\n\n\n\n\n","category":"method"},{"location":"library/library/#ParticleCorrelations.DiscretePairCorrelation","page":"Library","title":"ParticleCorrelations.DiscretePairCorrelation","text":"DiscretePairCorrelation\n\nRepresents the pair correlation between two types of particles. The particles could be the same or different species. \n\nThis struct has the field r, which represent the radial distance, and the field g is the value of the pair correlation. For radial distances < r[1] we assume that the pair correlation is zero. For radial distances > r[end] we assume the pair correlation is 1. \n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.DiscreteStructureFactor","page":"Library","title":"ParticleCorrelations.DiscreteStructureFactor","text":"DiscreteStructureFactor\n\nRepresents the structure factor between two types of particles. The particles could be the same or different species.\n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.HardMedium","page":"Library","title":"ParticleCorrelations.HardMedium","text":"HardMedium\n\nRepresents a solid type of material. Solid particles can not overlap. Currently the package has no soft medium, but this can easily be added.\n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.MonteCarloPairCorrelation","page":"Library","title":"ParticleCorrelations.MonteCarloPairCorrelation","text":"MonteCarloPairCorrelation{Dim} <: PairCorrelationType\n\nCurrently only used to create pair-correlations for particles that are uniformly randomly placed, except they can not overlap.\n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.PairCorrelation","page":"Library","title":"ParticleCorrelations.PairCorrelation","text":"PairCorrelation\n\nA type used to store a pair-correlation. This represents the calculated pair-correlation.\n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.PairCorrelationType","page":"Library","title":"ParticleCorrelations.PairCorrelationType","text":"PairCorrelationType\n\nA type used to specify what type of pair correlation is to be used. This is like a tag, or a option, to specify which pair correlation is wanted.\n\n\n\n\n\n","category":"type"},{"location":"library/library/#ParticleCorrelations.Specie","page":"Library","title":"ParticleCorrelations.Specie","text":"Specie\n\nRepresents a set of particles which are all the same. The type of particle is given by Specie.particle and the volume fraction this specie occupies is given by Specie.volume_fraction.\n\nWe can use Specie.numberofparticles to specify the number of particles, otherwise for an infinite Specie.numberofparticles = Inf.\n\nThe minimum distance between any two particles will equal outer_radius(Specie) * Specie.separation_ratio.\n\n\n\n\n\n","category":"type"},{"location":"#ParticleCorrelations.jl","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"","category":"section"},{"location":"","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"A Julia library to calculate pair correlations (or Structure factors) for disordered particulates, and calculate particle configurations from pair correlations.","category":"page"},{"location":"#Installation","page":"ParticleCorrelations.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"Install Julia v1.0 or later, then run","category":"page"},{"location":"","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"using Pkg # or enter the package mode by pressing ]\nPkg.add(\"ParticleCorrelations\")","category":"page"},{"location":"","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"Press backspace to exit the package mode.","category":"page"},{"location":"#Contents","page":"ParticleCorrelations.jl","title":"Contents","text":"","category":"section"},{"location":"","page":"ParticleCorrelations.jl","title":"ParticleCorrelations.jl","text":"Depth = 2","category":"page"}]
}