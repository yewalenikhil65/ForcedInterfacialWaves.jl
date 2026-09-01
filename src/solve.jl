"""
    CommonSolve.jl integration for ForcedInterfacialWaves.

    Provides a `solve(prob)` interface following the SciML convention:
    - `ForcedGravityProblem`  — pure-gravity IVP (α = 0)
    - `ForcedGCProblem`       — capillary–gravity IVP (α > 0)

    Both single-point and profile (vector) problems are supported via
    the same type — pass a scalar `x` or a vector `x_grid`.
"""

import CommonSolve: solve

# ═══════════════════════════════════════════════════════════════════════════════
# Solve method objects
# ═══════════════════════════════════════════════════════════════════════════════

"""
    IVP(; asym_cancel=false) → IVP

Marker method for the full time-dependent initial-value problem. For a
capillary–gravity problem, `asym_cancel=false` returns the documented
equation (4.5b) symmetric steady/transient decomposition. Setting
`asym_cancel=true` uses the asymmetric classical radiation steady reference,
so the returned transient is the remainder relative to that reference.
The option has no effect on the pure-gravity IVP dispatch.
"""
struct IVP
    asym_cancel::Bool
end

IVP(; asym_cancel::Bool=false) = IVP(asym_cancel)

"""Steady-state method configuration."""
struct SteadyMethod
    rayleigh_dissipation::Bool
end

"""
    steady(; rayleigh_dissipation=false) → SteadyMethod

Construct a steady-state solve method. Enabling Rayleigh dissipation selects
the radiation-based asymmetric steady response; disabling it returns the
symmetric Fourier steady response.
"""
steady(; rayleigh_dissipation::Bool=false) = SteadyMethod(rayleigh_dissipation)

function requireTime(t)
    t === nothing &&
        throw(ArgumentError("method=IVP() requires a nondimensional time t; construct the problem with t"))
    return t
end

"""Return `(η_s_local, η_s_farfield, η_s_total)` for pure gravity."""
function steadyComponents(x::Real, p::PureGravityParams,
                          rayleigh_dissipation::Bool)
    xa = abs(x)
    β² = p.beta^2
    localIntegral, _ = quadgk(y -> exp(-y * xa) * y / (β² + y^2),
                              0.0, Inf;
                              atol=p.atol_transformed, rtol=p.rtol_transformed)
    η_local = p.F0 * localIntegral / (π * (1.0 + p.rho_r))
    η_far_symmetric = -p.F0 * sin(p.beta * xa) / (1.0 + p.rho_r)
    η_farfield = rayleigh_dissipation ?
        η_far_symmetric - p.F0 * sin(p.beta * x) / (1.0 + p.rho_r) :
        η_far_symmetric
    return η_local, η_farfield, η_local + η_farfield
end

"""Return `(η_s_local, η_s_farfield, η_s_total)` for capillary gravity."""
function steadyComponents(x::Real, p::CapillaryGravityParams,
                          rayleigh_dissipation::Bool)
    G_x = cg_Gx_integral(x, p)
    η_local = p.F0 * G_x / (π * p.alpha)
    denominator = p.alpha * (p.k_l - p.k_s)
    η_farfield = if rayleigh_dissipation
        if x > 0.0
            -2.0 * p.F0 * sin(p.k_s * x) / denominator
        else
            -2.0 * p.F0 * sin(p.k_l * x) / denominator
        end
    else
        p.F0 / denominator *
        (-sin(p.k_s * abs(x)) + sin(p.k_l * abs(x)))
    end
    return η_local, η_farfield, η_local + η_farfield
end

function steadyComponents(x::AbstractVector{<:Real}, p::PureGravityParams,
                          rayleigh_dissipation::Bool)
    localIntegral = zeros(Float64, length(x))
    β² = p.beta^2
    localIntegrand = (values, y) -> begin
        factor = y / (β² + y^2)
        @inbounds @simd for ix in eachindex(x, values)
            values[ix] = exp(-y * abs(x[ix])) * factor
        end
        values
    end
    quadgk!(localIntegrand, localIntegral, 0.0, Inf;
            atol=p.atol_transformed, rtol=p.rtol_transformed,
            norm=values -> maximum(abs, values))

    localScale = p.F0 / (π * (1.0 + p.rho_r))
    η_local = localScale .* localIntegral
    η_far_symmetric = @. -p.F0 * sin(p.beta * abs(x)) / (1.0 + p.rho_r)
    η_farfield = if rayleigh_dissipation
        @. η_far_symmetric - p.F0 * sin(p.beta * x) / (1.0 + p.rho_r)
    else
        η_far_symmetric
    end
    return η_local, η_farfield, η_local .+ η_farfield
end

function steadyComponents(x::AbstractVector{<:Real}, p::CapillaryGravityParams,
                          rayleigh_dissipation::Bool)
    # This profile helper performs one threaded vector-valued quadrature per
    # chunk, instead of one scalar cg_Gx_integral call for every grid point.
    η_rayleigh = compute_cg_steady_profile(x, p)
    denominator = p.alpha * (p.k_l - p.k_s)
    η_far_rayleigh = Vector{Float64}(undef, length(x))
    @inbounds @simd for ix in eachindex(x, η_far_rayleigh)
        xv = x[ix]
        η_far_rayleigh[ix] = xv > 0.0 ?
            -2.0 * p.F0 * sin(p.k_s * xv) / denominator :
            -2.0 * p.F0 * sin(p.k_l * xv) / denominator
    end
    η_local = η_rayleigh .- η_far_rayleigh
    η_farfield = if rayleigh_dissipation
        η_far_rayleigh
    else
        @. p.F0 / denominator *
           (-sin(p.k_s * abs(x)) + sin(p.k_l * abs(x)))
    end
    return η_local, η_farfield, η_local .+ η_farfield
end

"""
    WaveSolution

Result returned by `solve`. Fields:

- `η` — total surface displacement (scalar or vector)
- `η_steady` — nominal time-independent component used in the IVP decomposition
  (the symmetric Fourier/PV term for capillary–gravity, equation (4.5b))
- `η_s_local` — localized part of the steady component
- `η_s_farfield` — far-field wave part of the steady component
- `η_transient` — IVP remainder `η - η_steady`
- `x` — spatial coordinate(s) used
- `t` — time used (`nothing` for a time-independent steady solve)
"""
struct WaveSolution{X,E}
    η::E
    η_steady::E
    η_s_local::E
    η_s_farfield::E
    η_transient::Union{E,Nothing}
    x::X
    t::Union{Nothing,Float64}
end

# ═══════════════════════════════════════════════════════════════════════════════
# Pure-gravity problem  (α = 0)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ForcedGravityProblem(params, x, t)
    ForcedGravityProblem(params, x)  # steady-only problem
    ForcedGravityProblem(; params, x, t=nothing)

Define a forced pure-gravity interfacial wave problem. The time `t` is
required for `method=IVP()` and may be omitted for `method=steady(rayleigh_dissipation=false)`.

- `params::PureGravityParams` — physical parameters (from [`compute_gravity_parameters`](@ref))
- `x` — spatial point (`Real`) or spatial grid (`AbstractVector{<:Real}`)
- `t::Real` — nondimensional time for an IVP; `nothing` for a steady-only problem

# Example

```julia
pg = compute_gravity_parameters()
prob_ivp = ForcedGravityProblem(pg, -2.0, 1.0 / pg.t_c)
prob_steady = ForcedGravityProblem(pg, make_gravity_xgrid(pg))
sol = solve(prob_steady; method=steady(rayleigh_dissipation=false))
```
"""
struct ForcedGravityProblem{X}
    params::PureGravityParams
    x::X
    t::Union{Nothing,Float64}
end

function ForcedGravityProblem(params::PureGravityParams, x, t::Real)
    return ForcedGravityProblem{typeof(x)}(params, x, Float64(t))
end

function ForcedGravityProblem(params::PureGravityParams, x)
    return ForcedGravityProblem{typeof(x)}(params, x, nothing)
end

function ForcedGravityProblem(; params::PureGravityParams, x, t=nothing)
    return t === nothing ? ForcedGravityProblem(params, x) : ForcedGravityProblem(params, x, t)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Capillary–gravity problem  (α > 0)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ForcedGCProblem(params, x, t)
    ForcedGCProblem(params, x)  # steady-only problem
    ForcedGCProblem(; params, x, t=nothing)

Define a forced capillary–gravity interfacial wave problem. The time `t` is
required for `method=IVP()` and may be omitted for `method=steady(rayleigh_dissipation=false)`.

- `params::CapillaryGravityParams` — physical parameters (from [`compute_cg_parameters`](@ref))
- `x` — spatial point (`Real`) or spatial grid (`AbstractVector{<:Real}`)
- `t::Real` — nondimensional time for an IVP; `nothing` for a steady-only problem

# Example

```julia
p = compute_cg_parameters()
prob_ivp = ForcedGCProblem(p, 3.0, 110.0)
prob_steady = ForcedGCProblem(p, make_cg_xgrid(p))
sol = solve(prob_steady; method=steady(rayleigh_dissipation=false))
```
"""
struct ForcedGCProblem{X}
    params::CapillaryGravityParams
    x::X
    t::Union{Nothing,Float64}
end

function ForcedGCProblem(params::CapillaryGravityParams, x, t::Real)
    return ForcedGCProblem{typeof(x)}(params, x, Float64(t))
end

function ForcedGCProblem(params::CapillaryGravityParams, x)
    return ForcedGCProblem{typeof(x)}(params, x, nothing)
end

function ForcedGCProblem(; params::CapillaryGravityParams, x, t=nothing)
    return t === nothing ? ForcedGCProblem(params, x) : ForcedGCProblem(params, x, t)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Typed solve dispatches
# ═══════════════════════════════════════════════════════════════════════════════

"""Solve using the selected typed method object."""
function solve(prob::ForcedGravityProblem{<:Real}; method=IVP())
    return solve(prob, method)
end

function solve(prob::ForcedGravityProblem{<:AbstractVector}; method=IVP())
    return solve(prob, method)
end

function solve(prob::ForcedGCProblem{<:Real}; method=IVP())
    return solve(prob, method)
end

function solve(prob::ForcedGCProblem{<:AbstractVector}; method=IVP())
    return solve(prob, method)
end

# Legacy Symbol spellings remain accepted while the public API uses objects.
function solve(prob::Union{ForcedGravityProblem,ForcedGCProblem}, method::Symbol)
    if method === :IVP
        return solve(prob, IVP())
    elseif method === :steady
        return solve(prob, steady(rayleigh_dissipation=true))
    end
    throw(ArgumentError("method must be IVP(), steady(...), :IVP, or :steady; got $(repr(method))"))
end

function solve(prob::Union{ForcedGravityProblem,ForcedGCProblem}, method)
    throw(ArgumentError("method must be IVP(), steady(...), :IVP, or :steady; got $(repr(method))"))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Pure gravity: scalar and profile methods
# ═══════════════════════════════════════════════════════════════════════════════

"""Pure-gravity full initial-value solution at one point."""
function solve(prob::ForcedGravityProblem{<:Real}, ::IVP)
    p = prob.params
    x = prob.x
    t = requireTime(prob.t)
    η_s_local, η_s_farfield, η_s = steadyComponents(x, p, false)
    a = t - x

    if abs(a) > p.front_band
        η_transient = p.F0 * (T₁(x, t, p) + T₂(x, t, p) +
                              T₃(x, t, p) + T₄(x, t, p))
        return WaveSolution(η_s + η_transient, η_s, η_s_local,
                            η_s_farfield, η_transient, x, t)
    end

    η_cpv = gravity_numerical_cpv(x, t, p)
    η_transient = η_cpv - η_s
    return WaveSolution(η_cpv, η_s, η_s_local, η_s_farfield,
                        η_transient, x, t)
end

"""Pure-gravity steady solution at one point."""
function solve(prob::ForcedGravityProblem{<:Real}, method::SteadyMethod)
    p = prob.params
    x = prob.x
    η_s_local, η_s_farfield, η_s =
        steadyComponents(x, p, method.rayleigh_dissipation)
    return WaveSolution(η_s, η_s, η_s_local, η_s_farfield,
                        zero(η_s), x, prob.t)
end

"""Pure-gravity full initial-value profile."""
function solve(prob::ForcedGravityProblem{<:AbstractVector}, ::IVP)
    p = prob.params
    x = prob.x
    t = requireTime(prob.t)
    η_ivp, η_steady, η_transient, _ = compute_gravity_profile(x, t, p)
    η_s_local, η_s_farfield, _ = steadyComponents(x, p, false)
    return WaveSolution(η_ivp, η_steady, η_s_local, η_s_farfield,
                        η_transient, x, t)
end

"""Pure-gravity steady profile."""
function solve(prob::ForcedGravityProblem{<:AbstractVector}, method::SteadyMethod)
    p = prob.params
    x = prob.x
    η_s_local, η_s_farfield, η_s =
        steadyComponents(x, p, method.rayleigh_dissipation)
    return WaveSolution(η_s, η_s, η_s_local, η_s_farfield,
                        zeros(Float64, length(x)), x, prob.t)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Capillary–gravity: scalar and profile methods
# ═══════════════════════════════════════════════════════════════════════════════

"""Capillary–gravity full initial-value solution at one point."""
function solve(prob::ForcedGCProblem{<:Real}, method::IVP)
    p = prob.params
    x = prob.x
    t = requireTime(prob.t)
    # asym_cancel=false is the symmetric nominal decomposition in (4.5b).
    # asym_cancel=true selects the asymmetric classical radiation reference.
    η_s_local, η_s_farfield, η_s =
        steadyComponents(x, p, method.asym_cancel)
    η_ivp = ivp_surface_elevation(x, t, p)
    η_transient = η_ivp - η_s
    return WaveSolution(η_ivp, η_s, η_s_local, η_s_farfield,
                        η_transient, x, t)
end

"""Capillary–gravity steady solution at one point."""
function solve(prob::ForcedGCProblem{<:Real}, method::SteadyMethod)
    p = prob.params
    x = prob.x
    η_s_local, η_s_farfield, η_s =
        steadyComponents(x, p, method.rayleigh_dissipation)
    return WaveSolution(η_s, η_s, η_s_local, η_s_farfield,
                        zero(η_s), x, prob.t)
end

"""Capillary–gravity full initial-value profile."""
function solve(prob::ForcedGCProblem{<:AbstractVector}, method::IVP)
    p = prob.params
    x = prob.x
    t = requireTime(prob.t)
    η_ivp, η_steady = compute_cg_profile(x, t, p;
                                         rayleigh_dissipation=method.asym_cancel)
    η_transient = η_ivp .- η_steady
    η_s_local, η_s_farfield, _ =
        steadyComponents(x, p, method.asym_cancel)
    return WaveSolution(η_ivp, η_steady, η_s_local, η_s_farfield,
                        η_transient, x, t)
end

"""Capillary–gravity steady profile."""
function solve(prob::ForcedGCProblem{<:AbstractVector}, method::SteadyMethod)
    p = prob.params
    x = prob.x
    η_s_local, η_s_farfield, η_s =
        steadyComponents(x, p, method.rayleigh_dissipation)
    return WaveSolution(η_s, η_s, η_s_local, η_s_farfield,
                        zeros(Float64, length(x)), x, prob.t)
end
