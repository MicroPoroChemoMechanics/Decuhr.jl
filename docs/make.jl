using Documenter
using Documenter.Remotes
using Pkg

# Activate the package so DECUHR is loadable from the docs environment
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))

using DECUHR

makedocs(
    sitename = "DECUHR.jl",
    modules  = [DECUHR],
    authors  = "Jean-François Barthélémy",
    remotes  = Dict(
        joinpath(@__DIR__, "..") => (Remotes.GitHub("MicMacTools", "DECUHR.jl"), "main"),
    ),
    format   = Documenter.HTML(
        prettyurls       = get(ENV, "CI", nothing) == "true",
        canonical        = nothing,
        assets           = String[],
        ansicolor        = true,
    ),
    pages = [
        "Home"           => "index.md",
        "Algorithm"      => "algorithm.md",
        "Examples"       => "examples.md",
        "API Reference"  => "api.md",
    ],
    checkdocs = :exports,   # warn only on exported names
    warnonly  = true,
)

deploydocs(
    repo      = "github.com/MicMacTools/DECUHR.jl.git",
    devbranch = "main",
)
