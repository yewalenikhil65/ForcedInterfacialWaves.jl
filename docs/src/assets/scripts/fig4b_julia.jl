using QuadGK
using Plots, LaTeXStrings
using DelimitedFiles

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

Δ  = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
kₗ = ((1 + ρᵣ) + √Δ) / (2α)
kₛ = ((1 + ρᵣ) - √Δ) / (2α)

η_local(x) = F₀ * (kₗ + kₛ) / (π * α) *
             first(quadgk(y -> y * exp(-abs(x) * y) / ((y^2 + kₗ^2) * (y^2 + kₛ^2)),  0, Inf; atol = 1e-10, rtol = 1e-8))

x = range(-15.0, 15.0; length = 2001)
x = filter(xi -> abs(xi) > 1e-12, collect(x))

# ─── Equation (3.12):  Rayleigh-dissipation far-field (asymmetric) ───
η_farfield_rayleigh(x) = -2F₀ / (α * (kₗ - kₛ)) * sin((x < 0 ? kₗ : kₛ) * x)

η_ff_r  = η_farfield_rayleigh.(x)
η_loc_r = η_local.(x)          # identical local term as in eqn. (3.11)

plot(x, η_loc_r .* 1e3; label="local", color="red", lw=3,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
     legend=:outerright)
plot!(x, η_ff_r .* 1e3; label="far-field", color="blue", lw=3)

savefig(joinpath(@__DIR__, "..", "steady_rayleigh.png"))
println("Saved Fig 5b(i) -> docs/src/assets/steady_rayleigh.png")

writedlm(joinpath(@__DIR__, "julia_ssl_rayleigh_local.csv"),    hcat(x, η_loc_r), ',')
writedlm(joinpath(@__DIR__, "julia_ssl_rayleigh_farfield.csv"), hcat(x, η_ff_r),  ',')
println("Saved Julia CSV data and steady_rayleigh.png")
