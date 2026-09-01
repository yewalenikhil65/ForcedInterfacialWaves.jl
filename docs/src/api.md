# API Reference

## Abstract type

```@docs
AbstractWaveParams
```

## Primary workflow: `solve`

The public entry point for both regimes. Construct a problem, then call `solve`:

```@docs
ForcedGravityProblem
ForcedGCProblem
solve
WaveSolution
```

### IVP decomposition options

```@docs
IVP
steady
```

For capillary–gravity problems, `IVP()` and `IVP(asym_cancel=false)` return the
symmetric decomposition from equations (4.5a–d). Use
`IVP(asym_cancel=true)` when the transient is referenced to the asymmetric
classical radiation profile produced by the long-time cancellation.


## Parameters

`CapillaryGravityParams` and `PureGravityParams` are **mutable**. Every physical field is
accessible under both an ASCII and a Unicode name (`p.alpha ≡ p.α`, `p.k_l ≡ p.kₗ`, etc.),
including for destructuring (`(; α, β, kₗ, kₛ) = p`).

!!! warning "Mutation does not recompute dependent fields"
    Setting a field in place (e.g. `p.alpha = 0.2`) does **not** re-derive quantities that
    depend on it (`k_l`, `k_s`, `F0`, ...). For any physically consistent change, construct a
    fresh parameter set with [`compute_cg_parameters`](@ref)/[`compute_gravity_parameters`](@ref)
    keywords instead of mutating an existing instance.

```@docs
CapillaryGravityParams
PureGravityParams
compute_cg_parameters
compute_gravity_parameters
```

## Pure gravity: unified `T₀`–`T₄` terms

Branch automatically on the side of the wavefront (`x` vs. `t`); no separate left/right calls
are needed.

```@docs
T₀
T₁
T₂
T₃
T₄
```

### Fresnel integrals

```@docs
fresnel_C
fresnel_S
```

### Numerical CPV (independent cross-check)

```@docs
gravity_combined_integrand
gravity_numerical_cpv
```

### Profiles

```@docs
compute_gravity_profile
make_gravity_xgrid
```

## Capillary–gravity: dispersion and integrands

```@docs
dispersion_chi
cg_integrand_steady
cg_integrand_I3
cg_integrand_I4
cg_combined_integrand
cg_I3
cg_I4
```

`𝕀₃` and `𝕀₄` are Unicode aliases for `cg_I3` and `cg_I4`. `cg_I3` uses the
standalone pole-split CPV convention; `cg_I4` uses the algebraically regularized
non-singular integrand. These component functions are diagnostic APIs; the
production `solve` path integrates the combined integrand.

### Quadrature and profiles

```@docs
cg_partial_integrals
ivp_surface_elevation
cg_Gx_integral
steady_surface_elevation
compute_cg_ivp_profile
compute_cg_I3_profile
compute_cg_I4_profile
compute_cg_steady_profile
compute_cg_profile
make_cg_xgrid
```

## Unicode utilities

```@docs
∫
```

## Backward compatibility

Older, pre-`solve` low-level pure-gravity entry points (`gravity_T0`,
`gravity_T1_left`/`gravity_T1_right` through `gravity_T4_left`/`gravity_T4_right`,
`gravity_analytical_left`, `gravity_analytical_right`) remain in the package as thin,
undocumented wrappers around [`T₀`](@ref)–[`T₄`](@ref) for compatibility with earlier scripts.
`dispersion_omega` is likewise kept as an undocumented `const` alias of
[`dispersion_chi`](@ref). None of these are part of the recommended API surface — new code
should use [`solve`](@ref)`(`[`ForcedGravityProblem`](@ref)`(...))`, the unified `T₀`–`T₄`
functions, and [`dispersion_chi`](@ref) directly.
