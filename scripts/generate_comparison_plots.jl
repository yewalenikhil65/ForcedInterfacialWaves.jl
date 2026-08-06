## generate_comparison_plots.jl
#  Generates publication-quality comparison figures.
#  If MATLAB CSV data exists in docs/src/assets/, overlays it.
#  Otherwise plots Julia-only curves.
#
#  Saves SVG and PNG to docs/src/assets/
#  Usage: julia --project=. scripts/generate_comparison_plots.jl

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using WavesDelta
using Plots, LaTeXStrings, Printf, DelimitedFiles

plot_font = "Computer Modern"
default(fontfamily=plot_font, margin=6Plots.mm, linewidth=3,
        framestyle=:box, label=nothing, color="blue", grid=false,
        fg_legend=false, background_color_legend=false)

        scalefontsizes(1.5)
const ASSETS = joinpath(@__DIR__, "..", "docs", "src", "assets")
mkpath(ASSETS)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 10: CG IVP at time_index=300
# ═══════════════════════════════════════════════════════════════════════════════

println("Computing CG profile (MATLAB-exact parameters, threaded vector quadrature)...")
p = compute_cg_parameters()  # k_max=Inf, atol=1e-10, rtol=1e-8
t = 3.0 / p.t_c
x_grid = collect(range(-5.0, 10.0; length=201))

η_ivp_jl = compute_cg_ivp_profile(x_grid, t, p)
η_steady_jl = compute_cg_steady_profile(x_grid, p)

plt = plot(x_grid, 1e3 .* η_steady_jl;
    label="Steady (Julia)", color=:black, ls=:solid, lw=2)
plot!(plt, x_grid, 1e3 .* η_ivp_jl;
    label="IVP (Julia)", color=:blue, ls=:dash, lw=2)

# Overlay MATLAB data if present
matlab_cg = joinpath(ASSETS, "matlab_cg_ivp_t300.csv")
if isfile(matlab_cg)
    data = readdlm(matlab_cg, ',')
    plot!(plt, data[:,1], 1e3 .* data[:,2];
        label="IVP (MATLAB)", color=:red, ls=:dot, lw=2)
end

xlabel!(plt, L"x"); ylabel!(plt, L"10^3\eta")
xlims!(plt, -5, 10); ylims!(plt, -6, 13); yticks!(plt, [-4, 0, 4, 8])


title!(plt, "Figure 10: CG IVP vs Steady (t\\_dim=3s)")

savefig(plt, joinpath(ASSETS, "fig10_comparison.svg"))
savefig(plt, joinpath(ASSETS, "fig10_comparison.png"))
println("  Saved fig10_comparison.{svg,png}")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6: PG analytical + CPV at t_dim=1.0s
# ═══════════════════════════════════════════════════════════════════════════════

println("Computing PG profile...")
pg = compute_gravity_parameters()
t_pg = 1.0 / pg.t_c
x_pg = make_gravity_xgrid(pg; Nx=201)

η_analytical, η_s, η_tr, η_cpv = compute_gravity_profile(x_pg, t_pg, pg)
mask = .!isnan.(η_analytical)

plt2 = plot(x_pg[mask], 1e3 .* η_analytical[mask];
    label=L"\eta_{\mathrm{analytical}}" * " (Julia)", color=:blue, ls=:dashdot, lw=2)
plot!(plt2, x_pg, 1e3 .* η_s;
    label=L"\eta_s" * " (Julia)", color=:black, ls=:solid, lw=1.5)
plot!(plt2, x_pg[mask], 1e3 .* η_tr[mask];
    label=L"\eta_{\mathrm{tr}}" * " (Julia)", color=:magenta, ls=:dot, lw=1.5)
scatter!(plt2, x_pg[mask][1:20:end], 1e3 .* η_cpv[mask][1:20:end];
    label=L"\eta_{\mathrm{CPV}}" * " (Julia)", mc=RGB(0.49, 0.18, 0.56),
    ms=4, msw=0, shape=:circle)

# Overlay MATLAB CPV if present
matlab_pg = joinpath(ASSETS, "matlab_pg_cpv_t1s.csv")
if isfile(matlab_pg)
    data = readdlm(matlab_pg, ',')
    valid = .!isnan.(data[:,2])
    scatter!(plt2, data[valid,1], 1e3 .* data[valid,2];
        label=L"\eta_{\mathrm{CPV}}" * " (MATLAB)", mc=:red,
        ms=2, msw=0.0, shape=:diamond)
end

vline!(plt2, [t_pg]; ls=:dash, lw=1.2, color=:gray, label="")
xlabel!(plt2, L"x"); ylabel!(plt2, L"10^3\eta")
ylims!(plt2, -5, 9); yticks!(plt2, [-4, 0, 4, 8])
title!(plt2, "Figure 6: PG Analytical vs CPV (t\\_dim=1s)")

savefig(plt2, joinpath(ASSETS, "fig6_comparison.svg"))
savefig(plt2, joinpath(ASSETS, "fig6_comparison.png"))
println("  Saved fig6_comparison.{svg,png}")

# ═══════════════════════════════════════════════════════════════════════════════
# Integrand plot
# ═══════════════════════════════════════════════════════════════════════════════

println("Computing integrand plot...")
k_vals = range(0.01, 12.0; length=500)
I_vals = [cg_combined_integrand(k, 3.0, 110.0, p) for k in k_vals]

plt3 = plot(k_vals, I_vals;
    xlabel=L"k", ylabel=L"\mathbb{I}(k;\;x{=}3,\;t{=}110)",
    title="Combined CG integrand", color=:blue, lw=2)
vline!(plt3, [p.k_s, p.k_l]; ls=:dash, lw=1.5, color=:red,
    label=L"k_s,\;k_l")

savefig(plt3, joinpath(ASSETS, "cg_integrand.svg"))
savefig(plt3, joinpath(ASSETS, "cg_integrand.png"))
println("  Saved cg_integrand.{svg,png}")

println("\nAll plots saved in: ", ASSETS)
println("If you have MATLAB CSV files, place them in the same folder and rerun.")
