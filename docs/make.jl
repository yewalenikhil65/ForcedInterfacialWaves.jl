using Documenter
using ForcedInterfacialWaves
using Plots  # Preload before @example evaluation (avoids Julia 1.12 world-age warnings).

makedocs(;
    sitename = "ForcedInterfacialWaves.jl",
    modules  = [ForcedInterfacialWaves],
    remotes  = nothing,
    doctest  = true,
    # :missing_docs only flags internal helper types/functions (e.g. CGProfileIntegrand,
    # cg_integrate_scalar) that are intentionally undocumented implementation details.
    # Everything else (cross-references, doctests, @example execution) must fail the build.
    warnonly = [:missing_docs],
    format   = Documenter.HTML(;
        prettyurls = true,
        repolink = "https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl",
        mathengine = Documenter.KaTeX(Dict(
            # The source derivation defines \vp as \phi (supplementary_vinod_jfm.tex:69).
            # KaTeX must receive that macro or it leaves every affected equation unrendered.
            :macros => Dict("\\vp" => "\\phi"),
        )),
        assets = ["assets/custom.css", "assets/custom.js"],
    ),
    pages = [
        "Home"                     => "index.md",
        "Theory"                   => [
            "Overview"                 => "theory/overview.md",
            "Steady-State"               => "theory/steady_state.md",
            "Pure Gravity"              => "theory/pure_gravity.md",
            "Capillary–Gravity"         => "theory/capillary_gravity.md",
        ],
        "API Reference"            => "api.md",
    ],
)
