using Documenter
using WavesDelta

makedocs(;
    sitename = "WavesDelta.jl",
    modules  = [WavesDelta],
    remotes  = nothing,
    warnonly = [:missing_docs, :cross_references],
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home"                     => "index.md",
        "Theory"                   => "theory.md",
        "Julia ↔ MATLAB"           => "julia_matlab.md",
        "Validation"               => "validation.md",
        "API Reference"            => "api.md",
    ],
)
