using QuadGK          # adaptive Gauss–Kronrod quadrature (∫₀^∞ …)
using Plots, LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=3,
        framestyle=:box, label=nothing, color="blue", grid=false,
        fg_legend=false, background_color_legend=false,
        guidefontsize=16, tickfontsize=14, legendfontsize=14)

# ─── Nondimensional parameters (see the Overview page) ───
U, g, T = 26.7046, 981.0, 72.0
ρₗ, ρᵤ  = 1.0, 0.001
l_c = U^2 / g
α   = T / (ρₗ * U^2 * l_c)
ρᵣ  = ρᵤ / ρₗ
F₀  = 0.01 * T / (ρₗ * U^2 * l_c)

# Steady capillary/gravity roots  k_{l,s}  (roots of  α k² − (1+ρᵣ)k + (1−ρᵣ) = 0)
Δ  = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
kₗ = ((1 + ρᵣ) + √Δ) / (2α)      # long-wavelength (capillary) root
kₛ = ((1 + ρᵣ) - √Δ) / (2α)      # short-wavelength (gravity) root

# ─── Equation (3.11):  ηₛ(x)/F₀ = far-field + local ───

# Far-field: closed-form sine combination
η_farfield(x) = F₀ / (α * (kₗ - kₛ)) * (-sin(kₛ * abs(x)) + sin(kₗ * abs(x)))

# Local: exponentially-decaying integral  (kₗ+kₛ)/(πα) ∫₀^∞ y e^{-|x|y} / [(y²+kₗ²)(y²+kₛ²)] dy
# quadgk returns (value, error); `first` keeps the value.
η_local(x) = F₀ * (kₗ + kₛ) / (π * α) *
             first(quadgk(y -> y * exp(-abs(x) * y) / ((y^2 + kₗ^2) * (y^2 + kₛ^2)),  0, Inf; atol = 1e-10, rtol = 1e-8))

# ─── Evaluate on a grid (skip x = 0, where the local integrand is singular) ───
x = range(-15.0, 15.0; length = 2001)
x = filter(xi -> abs(xi) > 1e-12, collect(x))

η_ff  = η_farfield.(x)
η_loc = η_local.(x)

plot(x, η_loc .* 1e3; label="local", color="red", lw=3,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
     legend=:outerright)
plot!(x, η_ff .* 1e3; label="far-field", color="blue", lw=3)

savefig(joinpath(@__DIR__, "..", "steady_no_rayleigh.png"))
println("Saved Fig 5a(i) -> docs/src/assets/steady_no_rayleigh.png")

# ─── Persist data for MATLAB overlay comparison ───
using DelimitedFiles
writedlm(joinpath(@__DIR__, "julia_ssl_no_rayleigh_local.csv"),    hcat(x, η_loc),  ',')
writedlm(joinpath(@__DIR__, "julia_ssl_no_rayleigh_farfield.csv"), hcat(x, η_ff),   ',')
println("Saved Julia CSV data and steady_no_rayleigh.png")
