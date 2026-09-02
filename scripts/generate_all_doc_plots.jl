"""
    generate_all_doc_plots.jl

Generates all plot PNGs needed by the documentation pages, replacing
the @example blocks that previously ran during makedocs.

Usage:
    julia -t auto --project=scripts scripts/generate_all_doc_plots.jl
"""

using ForcedInterfacialWaves
using Plots, LaTeXStrings
using DelimitedFiles

gr()

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=3, framestyle=:box, label=nothing,
        grid=false, fg_legend=false, background_color_legend=false)

ASSETS = joinpath(@__DIR__, "..", "docs", "src", "assets")

p  = compute_cg_parameters()
pg = compute_gravity_parameters()
x  = make_cg_xgrid(p; Nx=2001, xlim=(-15.0, 15.0))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Steady-state: no Rayleigh (Fig 5a Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
sol = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=false))
fig = plot(x, sol.η_s_local .* 1e3; label="local", color="red",
           guidefontsize=16, tickfontsize=14, legendfontsize=14,
           xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
           legend=:outerright, size=(800,400))
plot!(x, sol.η_s_farfield .* 1e3; label="far-field", color="blue")
savefig(fig, joinpath(ASSETS, "steady_no_rayleigh.png"))
println("1/6 steady_no_rayleigh.png")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Steady-state: Rayleigh dissipation (Fig 5b Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
sol_r = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=true))
fig = plot(x, sol_r.η_s_local .* 1e3; label="local", color="red",
           guidefontsize=16, tickfontsize=14, legendfontsize=14,
           xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
           legend=:outerright, size=(800,400))
plot!(x, sol_r.η_s_farfield .* 1e3; label="far-field", color="blue")
savefig(fig, joinpath(ASSETS, "steady_rayleigh.png"))
println("2/6 steady_rayleigh.png")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Pure-gravity IVP at t=183.68 (Fig 6 Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
x_pg = make_gravity_xgrid(pg; Nx=2001)
prof = solve(ForcedGravityProblem(pg, x_pg, 183.68))
fig = plot(x_pg, prof.η .* 1e3; label=L"\eta", color="blue", ls=:dash,
           guidefontsize=16, tickfontsize=14, legendfontsize=14,
           xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-12,12),
           ylims=(-4.8, 8.2), yticks=[-4,0,4,8],
           legend=:outerright, size=(800,400))
plot!(x_pg, prof.η_steady .* 1e3; label=L"\eta_s", color="black")
plot!(x_pg, prof.η_transient .* 1e3; label=L"\eta_{tr}", color="magenta", ls=:dot)
savefig(fig, joinpath(ASSETS, "pure_gravity_fig6.png"))
println("3/6 pure_gravity_fig6.png")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. CG I4 component at t=0.34 (Fig 7 Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
I4 = compute_cg_I4_profile(x, 0.34, p; method=:threaded_vector)
fig = plot(x, -(1/(2π)) .* I4 .* 1e3; color="purple",
           guidefontsize=16, tickfontsize=14,
           xlabel=L"x", ylabel=L"\frac{-\mathbb{I}_{4}}{2\pi} \times 10^{3}",
           xlims=(-10,10), size=(800,400))
savefig(fig, joinpath(ASSETS, "cg_I4_fig7.png"))
println("4/6 cg_I4_fig7.png")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. CG full IVP at t=367.35 (Fig 8 Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
sol_8 = solve(ForcedGCProblem(p, x, 367.35); method=IVP())
fig = plot(x, sol_8.η .* 1e3; label=L"\eta", color="blue", ls=:dash,
           guidefontsize=16, tickfontsize=14, legendfontsize=14,
           xlabel=L"x", ylabel=L"\eta \times 10^{3}",
           xlims=(-10,10), ylims=(-4.8, 8.2), yticks=[-4,0,4,8],
           legend=:outerright, size=(800,400))
plot!(x, sol_8.η_transient .* 1e3; label=L"\eta_{tr}", color="magenta", ls=:dot)
savefig(fig, joinpath(ASSETS, "cg_ivp_fig8.png"))
println("5/6 cg_ivp_fig8.png")

# ═══════════════════════════════════════════════════════════════════════════════
# 6. CG IVP vs simulation at t_dim=25s (Fig 10 Julia-only)
# ═══════════════════════════════════════════════════════════════════════════════
t_sim = 25 / (100 * p.t_c)
sol_sim = solve(ForcedGCProblem(p, x, t_sim); method=IVP())

data_bsk = sortslices(
    readdlm(joinpath(@__DIR__, "..", "notebooks", "if_25.csv"), ',', Float64; skipstart=1),
    dims=1, by=r -> r[6])
x_bsk = data_bsk[:, 6] ./ p.l_c
y_bsk = data_bsk[:, 7] ./ p.l_c

fig = plot(x, sol_sim.η .* 1e3; label=L"\eta", color="blue", ls=:dash,
           guidefontsize=16, tickfontsize=14, legendfontsize=14,
           xlabel=L"x", ylabel=L"\eta \times 10^{3}",
           xlims=(-6,10), ylims=(-6, 8.2), yticks=[-4,0,4,8],
           legend=:outerright, size=(800,400))
plot!(x_bsk, y_bsk .* 1e3; label="Simulation", color="red", ls=:dot)
savefig(fig, joinpath(ASSETS, "cg_sim_fig10.png"))
println("6/6 cg_sim_fig10.png")

println("Done.")
