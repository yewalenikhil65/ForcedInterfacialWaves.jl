"""
    Capillary-gravity two-fluid IVP solution (α > 0), Figure 10.

    The dispersion function is χ(k) = √(β k + γ_ρ α k³).
    The combined integrand 𝕀(x,t) sums the steady and two transient integrands
    so that the removable singularities at k_s, k_l cancel before quadrature.

    Individual integrands and partial integrals are exposed for open validation.
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Dispersion function
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dispersion_chi(k, p) → χ(k)

Two-fluid capillary–gravity dispersion function:

```math
\\chi(k) = \\sqrt{\\beta k + \\gamma_\\rho \\alpha k^3}
```
"""
@inline function dispersion_chi(k::Real, p::CapillaryGravityParams)
    return sqrt(p.beta * k + p.gamma_rho * p.alpha * k^3)
end

const dispersion_omega = dispersion_chi

# ═══════════════════════════════════════════════════════════════════════════════
# Individual integrands  (for validation: each can be evaluated at any k)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cg_integrand_steady(k, x, p)

The ``\\eta_s`` integrand (eq. 2.22 of supplementary):

```math
\\frac{2\\cos(kx)}{\\alpha(k - k_l)(k - k_s)}
```
"""
@inline function cg_integrand_steady(k::Real, x::Real, p::CapillaryGravityParams)
    return 2.0 * cos(k * x) / (p.alpha * (k - p.k_l) * (k - p.k_s))
end

"""
    cg_integrand_I3(k, x, t, p)

The ``\\mathbb{I}_3`` integrand (eq. 2.40 of supplementary):

```math
-\\frac{(1+\\rho_r)(k+\\chi)}{1+\\alpha k^2-\\rho_r}\\cdot
\\frac{\\cos[k(t-x)-t\\chi]}{\\alpha(k-k_l)(k-k_s)}
```
"""
@inline function cg_integrand_I3(k::Real, x::Real, t::Real, p::CapillaryGravityParams)
    χ = dispersion_chi(k, p)
    denom_poles = p.alpha * (k - p.k_l) * (k - p.k_s)
    denom_disp  = 1.0 + p.alpha * k^2 - p.rho_r
    phase = k * (t - x) - t * χ
    return -(1.0 + p.rho_r) * (k + χ) / denom_disp * cos(phase) / denom_poles
end

"""
    cg_integrand_I4(k, x, t, p)

The ``\\mathbb{I}_4`` integrand (eq. 2.41 of supplementary):

```math
-\\frac{(1+\\rho_r)(k-\\chi)}{1+\\alpha k^2-\\rho_r}\\cdot
\\frac{\\cos[k(t-x)+t\\chi]}{\\alpha(k-k_l)(k-k_s)}
```
"""
@inline function cg_integrand_I4(k::Real, x::Real, t::Real, p::CapillaryGravityParams)
    χ = dispersion_chi(k, p)
    denom_poles = p.alpha * (k - p.k_l) * (k - p.k_s)
    denom_disp  = 1.0 + p.alpha * k^2 - p.rho_r
    phase = k * (t - x) + t * χ
    return -(1.0 + p.rho_r) * (k - χ) / denom_disp * cos(phase) / denom_poles
end

"""
    cg_combined_integrand(k, x, t, p)

Combined integrand ``\\mathbb{I}(k; x, t)`` with singularities cancelling in the sum:

```math
\\mathbb{I} = \\frac{2\\cos(kx)}{\\alpha(k-k_l)(k-k_s)}
  - \\frac{(1+\\rho_r)(k+\\chi)}{1+\\alpha k^2 - \\rho_r}\\,
    \\frac{\\cos[k(t-x)-t\\chi]}{\\alpha(k-k_l)(k-k_s)}
  - \\frac{(1+\\rho_r)(k-\\chi)}{1+\\alpha k^2 - \\rho_r}\\,
    \\frac{\\cos[k(t-x)+t\\chi]}{\\alpha(k-k_l)(k-k_s)}
```

Uses `@fastmath`; computes ``\\chi(k)`` once and shares common subexpressions.
"""
@inline function cg_combined_integrand(k::Real, x::Real, t::Real, p::CapillaryGravityParams)
    @fastmath begin
        χ = sqrt(p.beta * k + p.gamma_rho * p.alpha * k^3)
        inv_denom_poles = 1.0 / (p.alpha * (k - p.k_l) * (k - p.k_s))
        inv_denom_disp  = (1.0 + p.rho_r) / (1.0 + p.alpha * k^2 - p.rho_r)

        kt_x = k * (t - x)
        tχ = t * χ

        I_s = 2.0 * cos(k * x) * inv_denom_poles
        I_3 = -inv_denom_disp * (k + χ) * cos(kt_x - tχ) * inv_denom_poles
        I_4 = -inv_denom_disp * (k - χ) * cos(kt_x + tχ) * inv_denom_poles

        return I_s + I_3 + I_4
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Individual partial integrals  (for open validation of each sub-interval)
# ═══════════════════════════════════════════════════════════════════════════════


"""Immutable callable scalar integrand used by concurrent reference quadratures."""
struct CGScalarIntegrand{TX<:Real,TT<:Real,P<:CapillaryGravityParams}
    x::TX
    t::TT
    p::P
end

@inline (integrand::CGScalarIntegrand)(k::Real) =
    cg_combined_integrand(k, integrand.x, integrand.t, integrand.p)
"""
    cg_partial_integrals(x, t, p) → (I1, I2, I3)

Returns the three partial integrals of the combined integrand ``\\mathbb{I}``
(eq. 2.17 of supplementary), split around the removable poles at ``k_s`` and ``k_l``:

```math
I_1 = \\int_0^{k_s-\\varepsilon} \\mathbb{I}\\,dk, \\quad
I_2 = \\int_{k_s+\\varepsilon}^{k_l-\\varepsilon} \\mathbb{I}\\,dk, \\quad
I_3 = \\int_{k_l+\\varepsilon}^{k_{\\max}} \\mathbb{I}\\,dk
```

These can be validated individually against MATLAB's `integral_1`, `integral_2`, `integral_3`.

For performance-critical profiles, use [`compute_cg_ivp_profile`](@ref) which evaluates all
spatial points simultaneously via in-place vector-valued quadrature.
"""
function cg_partial_integrals(x::Real, t::Real, p::CapillaryGravityParams)
    integrand = CGScalarIntegrand(x, t, p)
    ε = p.epsilon_pv

    I1, _ = quadgk(integrand, 0.0, p.k_s - ε; atol=p.atol, rtol=p.rtol)
    I2, _ = quadgk(integrand, p.k_s + ε, p.k_l - ε; atol=p.atol, rtol=p.rtol)
    I3, _ = quadgk(integrand, p.k_l + ε, p.k_max;
                   atol=p.atol, rtol=p.rtol, order=15)

    return I1, I2, I3
end

# ═══════════════════════════════════════════════════════════════════════════════
# Full IVP solution
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ivp_surface_elevation(x, t, p)

Evaluate ``\\eta(x,t)`` at a single spatial point (eq. 2.17 of supplementary):

```math
\\eta(x,t) = -\\frac{F_0}{2\\pi}\\left(I_1 + I_2 + I_3\\right)
```

where ``I_1, I_2, I_3`` are the partial integrals over `[0, k_s−ε]`, `[k_s+ε, k_l−ε]`, `[k_l+ε, ∞)`.
"""
function ivp_surface_elevation(x::Real, t::Real, p::CapillaryGravityParams)
    I1, I2, I3 = cg_partial_integrals(x, t, p)
    return -p.F0 / (2.0π) * (I1 + I2 + I3)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Steady solution (individual G(x) integral also exposed)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cg_Gx_integral(x, p) → G(x)

Auxiliary integral for the steady solution:

```math
G(x) = \\frac{1}{k_l - k_s}\\int_0^{k_{\\max}}
\\left[\\frac{\\cos(kx)}{k+k_s} - \\frac{\\cos(kx)}{k+k_l}\\right]dk
```

Exposed for independent validation against MATLAB.
"""
function cg_Gx_integral(x::Real, p::CapillaryGravityParams)
    integrand = k -> cos(k * x) / (k + p.k_s) - cos(k * x) / (k + p.k_l)
    G_val, _ = quadgk(integrand, 0.0, p.k_max_steady; atol=p.atol, rtol=p.rtol)
    return G_val / (p.k_l - p.k_s)
end

"""
    steady_surface_elevation(x, p)

Classical steady capillary–gravity solution (obtained via residue calculus):

```math
\\frac{\\eta_{\\text{steady}}(x)}{F_0} = \\begin{cases}
  -\\dfrac{2}{\\alpha(k_l-k_s)}\\sin(k_s x) + \\dfrac{G(x)}{\\pi\\alpha}, & x > 0 \\\\[6pt]
  -\\dfrac{2}{\\alpha(k_l-k_s)}\\sin(k_l x) + \\dfrac{G(x)}{\\pi\\alpha}, & x < 0
\\end{cases}
```

where ``G(x)`` is computed by [`cg_Gx_integral`](@ref).
"""
function steady_surface_elevation(x::Real, p::CapillaryGravityParams)
    G_x = cg_Gx_integral(x, p)

    if x > 0.0
        return p.F0 * (-2.0 / (p.alpha * (p.k_l - p.k_s)) * sin(p.k_s * x) +
                        G_x / (π * p.alpha))
    else
        return p.F0 * (-2.0 / (p.alpha * (p.k_l - p.k_s)) * sin(p.k_l * x) +
                        G_x / (π * p.alpha))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Profile computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cg_profile_integrand!(values, k, x_grid, t, p)

Evaluate the combined capillary-gravity integrand for every `x_grid` point in
place. Algebraically,

    𝕀(k;x,t) = C(k,t) cos(kx) + S(k,t) sin(kx),

so χ(k), both transient amplitudes, and both time phases are computed only once
per quadrature node. Each spatial point then needs one `sincos(kx)` pair.
"""
function cg_profile_integrand!(values::AbstractVector, k::Real,
                               x_grid::AbstractVector{<:Real}, t::Real,
                               p::CapillaryGravityParams)
    χ = dispersion_chi(k, p)
    pole_factor = inv(p.alpha * (k - p.k_l) * (k - p.k_s))
    transient_factor = -(1.0 + p.rho_r) * pole_factor /
                       (1.0 + p.alpha * k^2 - p.rho_r)

    amplitude_minus = transient_factor * (k + χ)
    amplitude_plus  = transient_factor * (k - χ)
    phase_minus = t * (k - χ)
    phase_plus  = t * (k + χ)

    sin_minus, cos_minus = sincos(phase_minus)
    sin_plus,  cos_plus  = sincos(phase_plus)

    cosine_coefficient = 2.0 * pole_factor +
                         amplitude_minus * cos_minus +
                         amplitude_plus * cos_plus
    sine_coefficient = amplitude_minus * sin_minus +
                       amplitude_plus * sin_plus

    @inbounds @simd for ix in eachindex(x_grid, values)
        sin_kx, cos_kx = sincos(k * x_grid[ix])
        values[ix] = cosine_coefficient * cos_kx + sine_coefficient * sin_kx
    end
    return values
end

"""Maximum absolute component norm used to enforce pointwise profile tolerances."""

"""Callable in-place profile integrand with uniform-grid recurrence metadata."""
struct CGProfileIntegrand{X<:AbstractVector,T<:Real,P<:CapillaryGravityParams}
    x_grid::X
    t::T
    p::P
    uniform::Bool
    spacing::Float64
end

function CGProfileIntegrand(x_grid::AbstractVector{<:Real}, t::Real,
                            p::CapillaryGravityParams)
    if length(x_grid) < 2
        return CGProfileIntegrand(x_grid, t, p, false, 0.0)
    end
    spacing = Float64(x_grid[firstindex(x_grid) + 1] - x_grid[firstindex(x_grid)])
    uniform = all(ix -> isapprox(x_grid[ix] - x_grid[ix - 1], spacing;
                                  rtol=1e-12, atol=1e-14),
                  (firstindex(x_grid) + 2):lastindex(x_grid))
    return CGProfileIntegrand(x_grid, t, p, uniform, spacing)
end

@inline function cg_profile_coefficients(k::Real, t::Real,
                                         p::CapillaryGravityParams)
    χ = dispersion_chi(k, p)
    pole_factor = inv(p.alpha * (k - p.k_l) * (k - p.k_s))
    transient_factor = -(1.0 + p.rho_r) * pole_factor /
                       (1.0 + p.alpha * k^2 - p.rho_r)
    amplitude_minus = transient_factor * (k + χ)
    amplitude_plus  = transient_factor * (k - χ)
    sin_minus, cos_minus = sincos(t * (k - χ))
    sin_plus,  cos_plus  = sincos(t * (k + χ))
    cosine_coefficient = 2.0 * pole_factor +
                         amplitude_minus * cos_minus +
                         amplitude_plus * cos_plus
    sine_coefficient = amplitude_minus * sin_minus +
                       amplitude_plus * sin_plus
    return cosine_coefficient, sine_coefficient
end

function cg_uniform_profile_integrand!(values::AbstractVector, k::Real,
                                       integrand::CGProfileIntegrand)
    cosine_coefficient, sine_coefficient =
        cg_profile_coefficients(k, integrand.t, integrand.p)
    sin_step, cos_step = sincos(k * integrand.spacing)
    block_length = 64
    first = firstindex(integrand.x_grid)
    final = lastindex(integrand.x_grid)

    while first <= final
        block_end = min(first + block_length - 1, final)
        sin_kx, cos_kx = sincos(k * integrand.x_grid[first])
        @inbounds for ix in first:block_end
            values[ix] = cosine_coefficient * cos_kx +
                         sine_coefficient * sin_kx
            sin_kx, cos_kx = (sin_kx * cos_step + cos_kx * sin_step,
                              cos_kx * cos_step - sin_kx * sin_step)
        end
        first = block_end + 1
    end
    return values
end

function (integrand::CGProfileIntegrand)(values::AbstractVector, k::Real)
    if integrand.uniform
        return cg_uniform_profile_integrand!(values, k, integrand)
    end
    return cg_profile_integrand!(values, k, integrand.x_grid,
                                 integrand.t, integrand.p)
end
@inline cg_profile_norm(values::AbstractVector) = maximum(abs, values)

"""
    cg_profile_partial_integrals(x_grid, t, p; atol=p.atol, rtol=p.rtol)

Compute `(I1, I2, I3)` for every spatial point with three shared adaptive,
in-place vector quadratures. `quadgk!` avoids allocating a vector at every
integrand evaluation, and the maximum norm makes `atol`/`rtol` apply to the
worst point in the profile rather than to an L² aggregate.
"""
function cg_profile_partial_integrals(x_grid::AbstractVector{<:Real}, t::Real,
                                      p::CapillaryGravityParams;
                                      atol::Real=p.atol, rtol::Real=p.rtol)
    length(x_grid) > 0 || throw(ArgumentError("x_grid must not be empty"))

    I1 = zeros(Float64, length(x_grid))
    I2 = similar(I1)
    I3 = similar(I1)
    integrand! = CGProfileIntegrand(x_grid, t, p)
    ε = p.epsilon_pv

    quadgk!(integrand!, I1, 0.0, p.k_s - ε;
            atol=atol, rtol=rtol, norm=cg_profile_norm)
    quadgk!(integrand!, I2, p.k_s + ε, p.k_l - ε;
            atol=atol, rtol=rtol, norm=cg_profile_norm)
    quadgk!(integrand!, I3, p.k_l + ε, p.k_max;
            atol=atol, rtol=rtol, order=15, norm=cg_profile_norm)

    return I1, I2, I3
end

"""Integrate one scalar capillary-gravity profile point at supplied tolerances."""
function cg_integrate_scalar(x::Real, t::Real, p::CapillaryGravityParams,
                             atol::Real, rtol::Real)
    f = CGScalarIntegrand(x, t, p)
    ε = p.epsilon_pv
    I1, _ = quadgk(f, 0.0, p.k_s - ε; atol=atol, rtol=rtol)
    I2, _ = quadgk(f, p.k_s + ε, p.k_l - ε; atol=atol, rtol=rtol)
    I3, _ = quadgk(f, p.k_l + ε, p.k_max;
                   atol=atol, rtol=rtol, order=15)
    return -p.F0 / (2.0π) * (I1 + I2 + I3)
end

"""Integrate one vector-valued spatial chunk and return its IVP profile."""
function cg_integrate_profile_chunk(x_chunk::AbstractVector{<:Real}, t::Real,
                                    p::CapillaryGravityParams,
                                    atol::Real, rtol::Real)
    I1, I2, I3 = cg_profile_partial_integrals(x_chunk, t, p;
                                               atol=atol, rtol=rtol)
    prefactor = -p.F0 / (2.0π)
    eta = similar(I1)
    @inbounds @simd for ix in eachindex(eta)
        eta[ix] = prefactor * (I1[ix] + I2[ix] + I3[ix])
    end
    return eta
end

"""
    compute_cg_ivp_profile(x_grid, t, p; method=:threaded_vector,
                           atol=p.atol, rtol=p.rtol)

Available algorithms:
- `:threaded_vector` (production default): one in-place vector quadrature per
  thread-local spatial chunk; shares adaptation and uses all Julia threads.
- `:vector`: one in-place vector quadrature for the complete profile.
- `:threaded_scalar`: independent per-x quadrature retained as the reference.
"""
function compute_cg_ivp_profile(x_grid::AbstractVector{<:Real}, t::Real,
                                p::CapillaryGravityParams;
                                method::Symbol=:threaded_vector,
                                atol::Real=p.atol, rtol::Real=p.rtol)
    if method === :threaded_vector
        N = length(x_grid)
        eta_ivp = Vector{Float64}(undef, N)
        chunk_count = min(Threads.nthreads(), N)
        chunk_length = cld(N, chunk_count)

        Threads.@threads :static for chunk_index in 1:chunk_count
            first_index = (chunk_index - 1) * chunk_length + 1
            last_index = min(chunk_index * chunk_length, N)
            if first_index <= last_index
                indices = first_index:last_index
                chunk_eta = cg_integrate_profile_chunk(@view(x_grid[indices]), t, p,
                                                       atol, rtol)
                copyto!(@view(eta_ivp[indices]), chunk_eta)
            end
        end
        return eta_ivp
    elseif method === :vector
        return cg_integrate_profile_chunk(x_grid, t, p, atol, rtol)
    elseif method === :threaded_scalar
        eta_ivp = Vector{Float64}(undef, length(x_grid))
        Threads.@threads :static for ix in eachindex(x_grid)
            eta_ivp[ix] = cg_integrate_scalar(x_grid[ix], t, p, atol, rtol)
        end
        return eta_ivp
    end
    throw(ArgumentError("method must be :threaded_vector, :vector, or :threaded_scalar"))
end

"""Callable in-place steady-profile integrand with uniform-grid metadata."""
struct CGSteadyProfileIntegrand{X<:AbstractVector,P<:CapillaryGravityParams}
    x_grid::X
    p::P
    uniform::Bool
    spacing::Float64
end

function CGSteadyProfileIntegrand(x_grid::AbstractVector{<:Real},
                                  p::CapillaryGravityParams)
    if length(x_grid) < 2
        return CGSteadyProfileIntegrand(x_grid, p, false, 0.0)
    end
    spacing = Float64(x_grid[firstindex(x_grid) + 1] - x_grid[firstindex(x_grid)])
    uniform = all(ix -> isapprox(x_grid[ix] - x_grid[ix - 1], spacing;
                                  rtol=1e-12, atol=1e-14),
                  (firstindex(x_grid) + 2):lastindex(x_grid))
    return CGSteadyProfileIntegrand(x_grid, p, uniform, spacing)
end

function (integrand::CGSteadyProfileIntegrand)(values::AbstractVector, k::Real)
    coefficient = (inv(k + integrand.p.k_s) - inv(k + integrand.p.k_l)) /
                  (integrand.p.k_l - integrand.p.k_s)
    if !integrand.uniform
        @inbounds @simd for ix in eachindex(integrand.x_grid, values)
            values[ix] = coefficient * cos(k * integrand.x_grid[ix])
        end
        return values
    end

    sin_step, cos_step = sincos(k * integrand.spacing)
    block_length = 64
    first = firstindex(integrand.x_grid)
    final = lastindex(integrand.x_grid)
    while first <= final
        block_end = min(first + block_length - 1, final)
        sin_kx, cos_kx = sincos(k * integrand.x_grid[first])
        @inbounds for ix in first:block_end
            values[ix] = coefficient * cos_kx
            sin_kx, cos_kx = (sin_kx * cos_step + cos_kx * sin_step,
                              cos_kx * cos_step - sin_kx * sin_step)
        end
        first = block_end + 1
    end
    return values
end

"""Integrate `G(x)` for one spatial chunk."""
function cg_integrate_steady_chunk(x_chunk::AbstractVector{<:Real},
                                   p::CapillaryGravityParams,
                                   atol::Real, rtol::Real)
    G = zeros(Float64, length(x_chunk))
    integrand! = CGSteadyProfileIntegrand(x_chunk, p)
    quadgk!(integrand!, G, 0.0, p.k_max_steady;
            atol=atol, rtol=rtol, norm=cg_profile_norm)
    return G
end

"""
    compute_cg_steady_profile(x_grid, p; atol=p.atol, rtol=p.rtol)

Compute the steady profile with thread-local in-place vector quadratures.
"""
function compute_cg_steady_profile(x_grid::AbstractVector{<:Real},
                                   p::CapillaryGravityParams;
                                   atol::Real=p.atol, rtol::Real=p.rtol)
    N = length(x_grid)
    N > 0 || throw(ArgumentError("x_grid must not be empty"))
    G = Vector{Float64}(undef, N)
    chunk_count = min(Threads.nthreads(), N)
    chunk_length = cld(N, chunk_count)

    Threads.@threads :static for chunk_index in 1:chunk_count
        first_index = (chunk_index - 1) * chunk_length + 1
        last_index = min(chunk_index * chunk_length, N)
        if first_index <= last_index
            indices = first_index:last_index
            chunk_G = cg_integrate_steady_chunk(@view(x_grid[indices]), p,
                                                atol, rtol)
            copyto!(@view(G[indices]), chunk_G)
        end
    end

    eta_steady = similar(G)
    common = p.F0 / (π * p.alpha)
    wave = -2.0 * p.F0 / (p.alpha * (p.k_l - p.k_s))
    @inbounds @simd for ix in eachindex(x_grid, eta_steady)
        x = x_grid[ix]
        eta_steady[ix] = wave * sin((x > 0 ? p.k_s : p.k_l) * x) + common * G[ix]
    end
    return eta_steady
end

"""
    compute_cg_profile(x_grid, t, p; method=:threaded_vector, compute_steady=true,
                       atol=p.atol, rtol=p.rtol)

Compatibility wrapper returning `(eta_ivp, eta_steady)`. Set
`compute_steady=false` at times where Figure 10 does not plot the steady curve;
the second return value is then `nothing`.
"""
function compute_cg_profile(x_grid::AbstractVector{<:Real}, t::Real,
                            p::CapillaryGravityParams;
                            method::Symbol=:threaded_vector, compute_steady::Bool=true,
                            atol::Real=p.atol, rtol::Real=p.rtol)
    eta_ivp = compute_cg_ivp_profile(x_grid, t, p;
                                     method=method, atol=atol, rtol=rtol)
    eta_steady = compute_steady ?
                 compute_cg_steady_profile(x_grid, p; atol=atol, rtol=rtol) :
                 nothing
    return eta_ivp, eta_steady
end

"""
    make_cg_xgrid(p; Nx=2001)

Nondimensional spatial grid, excluding x ≈ 0.
"""
function make_cg_xgrid(p::CapillaryGravityParams; Nx::Int=2001)
    xg = collect(range(-p.L / (2.0 * p.l_c), p.L / (2.0 * p.l_c); length=Nx))
    filter!(x -> abs(x) > 1.0e-12, xg)
    return xg
end
