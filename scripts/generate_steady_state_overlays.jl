"""
    generate_steady_state_overlays.jl

Generates two static overlay PNGs for the capillary-gravity steady state:
  fig5a_overlay.png — without Rayleigh dissipation (eqn. 3.11)
  fig5b_overlay.png — with Rayleigh dissipation (eqn. 3.12)

Julia: continuous lines.
MATLAB: circle markers at every 40th Julia grid point (interpolated from CSV),
        so markers are noticeable but do not crowd the image.

Usage (from repo root):
    julia --project=scripts scripts/generate_steady_state_overlays.jl
"""

using ForcedInterfacialWaves
using Plots, LaTeXStrings
using DelimitedFiles

gr()

# ── Plot defaults (matching the docs page) ────────────────────────────────────
plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=3, framestyle=:box,
        grid=false, fg_legend=false, background_color_legend=false)

ASSETS = joinpath(@__DIR__, "..", "docs", "src", "assets")

# ── Parameters and Julia grid ─────────────────────────────────────────────────
p = compute_cg_parameters()
x = make_cg_xgrid(p; Nx=2001, xlim=(-15.0, 15.0))   # 2000 pts after filtering x=0

# ── Julia solutions ───────────────────────────────────────────────────────────
sol   = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=false))
sol_r = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=true))

# ── Load MATLAB CSVs and interpolate onto Julia grid ─────────────────────────
function load_and_interp(name, xj)
    data = readdlm(joinpath(ASSETS, name), ',', Float64)
    xm = data[:, 1]
    ym = data[:, 2]
    # linear interpolation at each Julia grid point
    [begin
        i = searchsortedfirst(xm, xi)
        if i == 1
            ym[1]
        elseif i > length(xm)
            ym[end]
        else
            t = (xi - xm[i-1]) / (xm[i] - xm[i-1])
            ym[i-1] * (1 - t) + ym[i] * t
        end
    end for xi in xj]
end

mat_nr_ff  = load_and_interp("ssl_no_rayleigh_farfield.csv", x)
mat_nr_loc = load_and_interp("ssl_no_rayleigh_local.csv",    x)
mat_r_ff   = load_and_interp("ssl_rayleigh_farfield.csv",    x)
mat_r_loc  = load_and_interp("ssl_rayleigh_local.csv",       x)

# ── Marker indices: local extrema of far-field, thinned to ≤50 total ─────────
function signal_markers(y)
    return 1:25:length(y)
end

# ── Helper: make one overlay figure ──────────────────────────────────────────
function make_fig(xj, julia_ff, julia_loc, mat_ff, mat_loc)
    mi = signal_markers(julia_ff)   # markers at peaks, troughs and zero-crossings
    fig = plot(; size=(800, 400),
                 guidefontsize=16, tickfontsize=14, legendfontsize=13,
                 xlabel=L"x", ylabel=L"\eta \times 10^{3}",
                 xlims=(-10, 10), legend=:outerright)

    # Julia: continuous dashed lines
    plot!(fig, xj, julia_ff  .* 1e3; color="magenta", ls=:dash,
          label=L"\mathrm{Julia:}\ \eta_{\mathrm{far\text{-}field}}")
    plot!(fig, xj, julia_loc .* 1e3; color="green",   ls=:dash,
          label=L"\mathrm{Julia:}\ \eta_{\mathrm{local}}")

    # MATLAB: markers at extrema
    scatter!(fig, xj[mi], mat_ff[mi]  .* 1e3;
             shape=:circle, ms=5, msw=1.2,
             mc=:blue, markerstrokecolor=:blue,
             label=L"\mathrm{MATLAB:}\ \eta_{\mathrm{far\text{-}field}}")
    # local term: uniform sparse markers (it's nearly flat)
    mi_loc = 1:25:length(xj)
    scatter!(fig, xj[mi_loc], mat_loc[mi_loc] .* 1e3;
             shape=:circle, ms=5, msw=1.2,
             mc=:red, markerstrokecolor=:red,
             label=L"\mathrm{MATLAB:}\ \eta_{\mathrm{local}}")
    return fig
end

# ── Fig. 5a — no Rayleigh dissipation ────────────────────────────────────────
fig5a = make_fig(x, sol.η_s_farfield, sol.η_s_local,
                 mat_nr_ff, mat_nr_loc)
savefig(fig5a, joinpath(ASSETS, "fig5a_overlay.png"))
println("Saved fig5a_overlay.png")

# ── Fig. 5b — Rayleigh dissipation ───────────────────────────────────────────
fig5b = make_fig(x, sol_r.η_s_farfield, sol_r.η_s_local,
                 mat_r_ff, mat_r_loc)
savefig(fig5b, joinpath(ASSETS, "fig5b_overlay.png"))
println("Saved fig5b_overlay.png")
println("Done.")
