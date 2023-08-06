using Documenter
using ParticleCorrelations

makedocs(
    sitename = "ParticleCorrelations",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [ParticleCorrelations],
    pages = [
        "Library" => "library/library.md"
    ]

)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/arturgower/ParticleCorrelations.jl.git",
)
