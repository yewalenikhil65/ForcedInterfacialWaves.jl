using Documenter

makedocs(;
    sitename = "Interfacial Waves from Pressure Forcing",
    remotes  = nothing,
    format   = Documenter.HTML(;
        prettyurls = true,
        repolink = "https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl",
        inventory_version = "0.1.0",
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
            "Basilisk setup for forced Capillary-Gravity waves"        => "theory/basilisk_capillary_gravity.md",
        ],
        "Validation"               => [
            "Nonlinear-regime : Capillary-Gravity waves evolution"  => "validation/cfd_conformal_mapping.md",
        ],
    ],
)
