"""
Pure-gravity two-fluid IVP (α = 0), manuscript Figure 6.

The numerical core follows the package implementation: every quadrature receives
concrete Float64 arguments and a typed PureGravityParams object. The complete
profile calculation is inside computeGravityProfile, not global scope.
"""

using QuadGK, FresnelIntegrals
using DelimitedFiles
using Plots, LaTeXStrings

"""Concrete physical and quadrature parameters for the pure-gravity IVP."""
struct PureGravityParams
    ρᵣ::Float64
    β::Float64
    sqrtβ::Float64
    F₀::Float64
    εCPV::Float64
    atolCPV::Float64
    rtolCPV::Float64
    atolTransformed::Float64
    rtolTransformed::Float64
    kₘₐₓAnalytical::Float64
    kₘₐₓCPV::Float64
    frontBand::Float64
end

function makeParams()
    U, g = 26.7046, 981.0
    ρₗ, ρᵤ = 1.0, 0.001
    ρᵣ = ρᵤ / ρₗ
    β = (1.0 - ρᵣ) / (1.0 + ρᵣ)
    ℓc = U^2 / g
    F₀ = 0.01 * 72.0 / (ρₗ * U^2 * ℓc)
    return PureGravityParams(ρᵣ, β, sqrt(β), F₀,
                             1e-6, 1e-10, 1e-8, 1e-10, 1e-8,
                             100.0, 100.0, 1.0)
end

# ─── T₀–T₄: scalar, type-stable numerical functions ───────────────────────────
function T₀(x::Float64, p::PureGravityParams)
    xₐ = abs(x)
    β² = p.β^2
    Ilocal, _ = quadgk(y -> exp(-y * xₐ) * y / (β² + y^2), 0.0, Inf;
                       atol=p.atolTransformed, rtol=p.rtolTransformed)
    return (-π * sin(p.β * xₐ) + Ilocal) / (π * (1.0 + p.ρᵣ))
end

@inline function T₁(x::Float64, t::Float64, p::PureGravityParams)
    s = x < t ? -1.0 : 1.0
    return s * sin(p.β * x) / (1.0 + p.ρᵣ)
end

function T₂(x::Float64, t::Float64, p::PureGravityParams)
    a = t - x
    s = a > 0.0 ? 1.0 : -1.0
    tSqrtβ = t * p.sqrtβ
    integrand = v -> begin
        arg = v * tSqrtβ
        numerator = p.sqrtβ * cos(arg) + s * (2.0 * v - p.sqrtβ) * sin(arg)
        denominator = p.β + (2.0 * v - p.sqrtβ)^2
        v^2 * exp(-s * 2.0 * v^2 * a + s * v * tSqrtβ) * numerator / denominator
    end
    I, _ = quadgk(integrand, 0.0, p.kₘₐₓAnalytical;
                  atol=p.atolTransformed, rtol=p.rtolTransformed)
    return -4.0 * I / (π * (1.0 + p.ρᵣ) * p.β)
end

function T₃(x::Float64, t::Float64, p::PureGravityParams)
    a = t - x
    absA = abs(a)
    s = a > 0.0 ? 1.0 : -1.0
    b = 0.5 * t * p.sqrtβ
    X = b * sqrt(2.0 / (π * absA))
    prefactor = (1.0 + t / (2.0 * a)) * sqrt(π / (2.0 * absA)) /
                (π * (1.0 + p.ρᵣ) * p.sqrtβ)
    phase = b^2 / absA
    return prefactor * (cos(phase) * (0.5 - s * fresnelc(X)) +
                        sin(phase) * (0.5 - s * fresnels(X)))
end

function T₄(x::Float64, t::Float64, p::PureGravityParams)
    a = t - x
    tSqrtβ = t * p.sqrtβ
    integrand = v -> cos(v^2 * a + v * tSqrtβ) / (v + p.sqrtβ)
    I, _ = quadgk(integrand, 0.0, p.kₘₐₓAnalytical;
                  atol=p.atolTransformed, rtol=p.rtolTransformed, order=15)
    return -I / (π * (1.0 + p.ρᵣ))
end

# ─── Direct numerical CPV ─────────────────────────────────────────────────────
@inline function combinedIntegrand(k::Float64, x::Float64, t::Float64,
                                   p::PureGravityParams)
    sqrtβk = sqrt(p.β * k)
    phase = k * (t - x)
    tSqrtβk = t * sqrtβk
    return cos(k * x) / (π * (1.0 + p.ρᵣ) * (k - p.β)) -
           k * cos(phase - tSqrtβk) / (2.0 * π * (1.0 - p.ρᵣ) * (k - sqrtβk)) -
           k * cos(phase + tSqrtβk) / (2.0 * π * (1.0 - p.ρᵣ) * (k + sqrtβk))
end

function ηCPV(x::Float64, t::Float64, p::PureGravityParams)
    G = k -> combinedIntegrand(k, x, t, p)
    I₋, _ = quadgk(G, 0.0, p.β - p.εCPV; atol=p.atolCPV, rtol=p.rtolCPV)
    I₊, _ = quadgk(G, p.β + p.εCPV, p.kₘₐₓCPV;
                   atol=p.atolCPV, rtol=p.rtolCPV, order=15)
    return p.F₀ * (I₋ + I₊)
end

"""
    computeGravityProfile(xgrid, t, p)

Compute `(η, ηₛ, ηₜᵣ, ηCPV)` for the spatial profile. Each iteration owns its
QuadGK calls, so static outer threading is safe and matches the package
pure-gravity implementation. Points in `|t-x| ≤ p.frontBand` are `NaN`.
"""
function computeGravityProfile(xgrid::Vector{Float64}, t::Float64,
                               p::PureGravityParams)
    N = length(xgrid)
    η = fill(NaN, N)
    ηₛ = Vector{Float64}(undef, N)
    ηₜᵣ = fill(NaN, N)
    ηCPVvalues = fill(NaN, N)

    Threads.@threads :static for i in eachindex(xgrid)
        x = xgrid[i]
        ηsteady = p.F₀ * T₀(x, p)
        ηₛ[i] = ηsteady

        if abs(t - x) > p.frontBand
            ηtransient = p.F₀ * (T₁(x, t, p) + T₂(x, t, p) +
                                 T₃(x, t, p) + T₄(x, t, p))
            ηₜᵣ[i] = ηtransient
            η[i] = ηsteady + ηtransient
            ηCPVvalues[i] = ηCPV(x, t, p)
        end
    end
    return η, ηₛ, ηₜᵣ, ηCPVvalues
end

function makeXGrid(; Nₓ::Int=2001, L::Float64=24.0)
    xgrid = collect(range(-L / 2.0, L / 2.0; length=Nₓ))
    filter!(x -> abs(x) > 1e-6, xgrid)
    return xgrid
end

function makeFigure(xgrid::Vector{Float64}, η::Vector{Float64},
                    ηₛ::Vector{Float64}, ηₜᵣ::Vector{Float64},
                    ηCPVvalues::Vector{Float64})
    plt = plot(xgrid, η .* 1e3; label=L"\eta", color="blue", linestyle=:dash,
               xlabel=L"x", ylabel=L"\eta \times 10^{3}", legend=:outerright,
               ylims=(-4.8, 8.2), yticks=[-4, 0, 4, 8], xlim=(-12, 12))
    plot!(plt, xgrid, ηCPVvalues .* 1e3; label=L"\eta_{\mathrm{CPV}}", color="red", linestyle=:dashdot)
    plot!(plt, xgrid, ηₛ .* 1e3; label=L"\eta_s", color="black")
    plot!(plt, xgrid, ηₜᵣ .* 1e3; label=L"\eta_{tr}", color="magenta", linestyle=:dot)
    return plt
end

function main(; tFig6::Float64=183.68, Nₓ::Int=2001)
    default(fontfamily="Computer Modern", linewidth=3, framestyle=:box,
            label=nothing, color="blue", grid=false, fg_legend=false,
            background_color_legend=false, guidefontsize=16,
            tickfontsize=14, legendfontsize=14)

    p = makeParams()
    xgrid = makeXGrid(Nₓ=Nₓ)
    η, ηₛ, ηₜᵣ, ηCPVvalues = computeGravityProfile(xgrid, tFig6, p)

    overlayCSV = joinpath(@__DIR__, "julia_fig6_profiles.csv")
    writedlm(overlayCSV, hcat(xgrid, η, ηCPVvalues, ηₛ, ηₜᵣ), ',')
    println("Saved Julia Figure 6 profiles -> ", overlayCSV)

    valid = .!isnan.(η) .& .!isnan.(ηCPVvalues)
    println("Threads: ", Threads.nthreads())
    println("max |η - ηCPV| (excluding front band) = ",
            maximum(abs.(η[valid] .- ηCPVvalues[valid])))

    figure = makeFigure(xgrid, η, ηₛ, ηₜᵣ, ηCPVvalues)
    output = joinpath(@__DIR__, "..", "pure_gravity_fig6.png")
    savefig(figure, output)
    println("Saved Fig. 6 -> ", output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
