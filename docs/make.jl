using Documenter
using ForcedInterfacialWaves

makedocs(;
    sitename = "ForcedInterfacialWaves.jl",
    modules  = [ForcedInterfacialWaves],
    remotes  = nothing,
    warnonly = [:missing_docs, :cross_references],
    format   = Documenter.HTML(;
        prettyurls = true,
        mathengine = Documenter.KaTeX(),
        assets = ["assets/custom.css"],
    ),
    pages = [
        "Home"                     => "index.md",
        "Theory"                   => "theory.md",
        "Julia ↔ MATLAB"           => "julia_matlab.md",
        "Validation"               => "validation.md",
        "API Reference"            => "api.md",
    ],
)
