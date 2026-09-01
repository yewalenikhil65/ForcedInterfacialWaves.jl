module ForcedInterfacialWaves

using QuadGK: quadgk, quadgk!, alloc_segbuf
using FresnelIntegrals: fresnelc, fresnels
using CommonSolve: solve

include("parameters.jl")
include("capillary_gravity.jl")
include("pure_gravity.jl")
include("solve.jl")

# ─── Abstract type ────────────────────────────────────────────────────────────
export AbstractWaveParams

# ─── Parameters ───────────────────────────────────────────────────────────────
export CapillaryGravityParams, PureGravityParams
export compute_cg_parameters, compute_gravity_parameters

# ─── CommonSolve problem types ────────────────────────────────────────────────
export ForcedGravityProblem, ForcedGCProblem, WaveSolution
export IVP, steady
export solve

# ─── Unicode convenience ──────────────────────────────────────────────────────

"""
    ∫(f, a, b; kw...) → value

Unicode alias for numerical integration. Returns the integral value only
(discards the error estimate). Wraps `QuadGK.quadgk`.

# Example

```julia
using ForcedInterfacialWaves: ∫
result = ∫(sin, 0, π)   # ≈ 2.0
```
"""
∫(f, a, b; kw...) = first(quadgk(f, a, b; kw...))
∫(f, a, b, c...; kw...) = first(quadgk(f, a, b, c...; kw...))

export ∫

# ─── Capillary-gravity (α > 0) ───────────────────────────────────────────────
export dispersion_chi, dispersion_omega, χ
export cg_integrand_steady, cg_integrand_I3, cg_integrand_I4
export cg_I3, cg_I4, 𝕀₃, 𝕀₄
export cg_combined_integrand, cg_partial_integrals
export cg_Gx_integral
export ivp_surface_elevation, steady_surface_elevation
export cg_profile_integrand!, cg_profile_partial_integrals
export compute_cg_ivp_profile, compute_cg_steady_profile
export compute_cg_I3_profile, compute_cg_I4_profile
export compute_cg_profile, make_cg_xgrid

# ─── Pure gravity (α = 0) ────────────────────────────────────────────────────
export T₀, T₁, T₂, T₃, T₄
export fresnel_C, fresnel_S
export gravity_T0
export gravity_T1_left, gravity_T2_left, gravity_T3_left, gravity_T4_left
export gravity_T1_right, gravity_T2_right, gravity_T3_right, gravity_T4_right
export gravity_analytical_left, gravity_analytical_right
export gravity_combined_integrand, gravity_numerical_cpv
export compute_gravity_profile, make_gravity_xgrid

end # module
