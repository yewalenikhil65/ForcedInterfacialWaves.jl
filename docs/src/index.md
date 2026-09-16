# ForcedInterfacialWaves.jl

*Julia/MATLAB implementation of the initial value problem (IVP) for pressure-forced interfacial waves in a two-fluid system.*


## Dependencies

The theory pages use the following Julia packages. Install them from the Julia REPL:

```julia
using Pkg

# Core numerical dependencies (used in all code blocks)
Pkg.add("QuadGK")           # adaptive Gauss–Kronrod quadrature
Pkg.add("FresnelIntegrals") # closed-form Fresnel C and S integrals

# Plotting (used in all profile and figure blocks)
Pkg.add("Plots")
Pkg.add("LaTeXStrings")

```

Requires Julia ≥ 1.9. All packages above are registered in the General registry and install with a single `Pkg.add` call.



## Contents

```@contents
Pages = [
    "theory/overview.md",
    "theory/steady_state.md",
    "theory/pure_gravity.md",
    "theory/capillary_gravity.md",
]
Depth = 2
```
