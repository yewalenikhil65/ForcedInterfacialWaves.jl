"""
    figure6_pure_gravity.jl

    Pure-gravity (α = 0) IVP: analytical vs numerical CPV evaluation.
    Reproduces Figure 6 of the manuscript at multiple time snapshots.

    Usage:
        julia -t auto --project=scripts scripts/figure6_pure_gravity.jl
"""

# ─── Thread check ──────────────────────────────────────────────────────────────
if Threads.nthreads() == 1
    @warn "Running single-threaded. For speedup: julia -t auto --project=scripts scripts/figure6_pure_gravity.jl"
else
    @info "Using $(Threads.nthreads()) threads"
end

using ForcedInterfacialWaves
using Plots; gr()

# ═══════════════════════════════════════════════════════════════════════════════
# Parameters
# ═══════════════════════════════════════════════════════════════════════════════

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "figure6")
mkpath(OUTPUT_DIR)

const TIME_INDICES   = [1, 7, 15, 40, 60, 100, 300, 500]
const NX_PLOT        = 2001
const MARKER_SPACING = 50

# ═══════════════════════════════════════════════════════════════════════════════
# Setup
# ═══════════════════════════════════════════════════════════════════════════════

pg = compute_gravity_parameters()

println("Pure-gravity parameters:")
println("  β         = ", pg.β)
println("  √β        = ", pg.sqrtβ)
println("  F₀        = ", pg.F₀)
println("  λ_gravity = ", pg.gravity_wavelength)
println("  L         = ", pg.L)

x_grid = make_gravity_xgrid(pg; Nx=NX_PLOT)
half_L = pg.L / 2.0

# ═══════════════════════════════════════════════════════════════════════════════
# Loop over times
# ═══════════════════════════════════════════════════════════════════════════════

for ti in TIME_INDICES
    t_dim = ti / 100.0
    t = t_dim / pg.t_c

    println("\n  time_index = $ti, t_dim = $(t_dim) s, t = $t ...")

    # Solve via CommonSolve API
    sol = solve(ForcedGravityProblem(pg, x_grid, t))
    η_analytical = sol.η
    η_steady     = sol.η_steady
    η_tr         = sol.η_transient

    # Independent CPV verification (pointwise)
    η_cpv = [gravity_numerical_cpv(x, t, pg) for x in x_grid]

    # Error metric (outside front band only)
    valid = .!isnan.(η_analytical) .& .!isnan.(η_cpv)
    if any(valid)
        max_err = maximum(abs.(η_analytical[valid] .- η_cpv[valid]))
        println("    max|η_ana − η_CPV| = ", max_err)
    end

    # ─── Masks for left / right of front ───
    mask_left  = (t .- x_grid) .>  pg.front_band
    mask_right = (t .- x_grid) .< -pg.front_band
    idx_l = findall(mask_left)
    idx_r = findall(mask_right)

    scale = 1.0e3

    # ─── Plot ───
    plt = plot(;
        size       = (500, 340),
        xlims      = (-half_L, half_L),
        ylims      = (-5.0, 9.0),
        yticks     = [-4, 0, 4, 8],
        xlabel     = "\$x\$",
        ylabel     = "\$10^3 \\eta\$",
        legend     = :topleft,
        framestyle = :box,
        grid       = false,
        tickdir    = :out,
    )

    # Analytical total (blue dash-dot)
    if !isempty(idx_l)
        plot!(plt, x_grid[idx_l], scale .* η_analytical[idx_l];
              lw=2.0, ls=:dashdot, lc=:blue, label="\$\\eta\$")
    end
    if !isempty(idx_r)
        plot!(plt, x_grid[idx_r], scale .* η_analytical[idx_r];
              lw=2.0, ls=:dashdot, lc=:blue, label="")
    end

    # Steady (black solid)
    plot!(plt, x_grid, scale .* η_steady;
          lw=1.5, ls=:solid, lc=:black, label="\$\\eta_s\$")

    # Transient (magenta dotted)
    if !isempty(idx_l)
        plot!(plt, x_grid[idx_l], scale .* η_tr[idx_l];
              lw=1.5, ls=:dot, lc=:magenta, label="\$\\eta_{\\mathrm{tr}}\$")
    end
    if !isempty(idx_r)
        plot!(plt, x_grid[idx_r], scale .* η_tr[idx_r];
              lw=1.5, ls=:dot, lc=:magenta, label="")
    end

    # CPV markers (purple circles, sparse)
    cpv_color = RGB(0.49, 0.18, 0.56)
    if !isempty(idx_l)
        mi = idx_l[1:MARKER_SPACING:end]
        scatter!(plt, x_grid[mi], scale .* η_cpv[mi];
                 mc=cpv_color, ms=3, msw=0.8, shape=:circle,
                 markerstrokecolor=cpv_color, label="\$\\eta_{\\mathrm{CPV}}\$")
    end
    if !isempty(idx_r)
        mi = idx_r[1:MARKER_SPACING:end]
        scatter!(plt, x_grid[mi], scale .* η_cpv[mi];
                 mc=cpv_color, ms=3, msw=0.8, shape=:circle,
                 markerstrokecolor=cpv_color, label="")
    end

    # Vertical line at x = t
    vline!(plt, [t]; ls=:dash, lc=:black, lw=1.2, label="")

    # ─── Save ───
    fname = "figure6_t" * lpad(ti, 4, '0')
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".png"))
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".pdf"))
    println("    saved: ", fname)
end

println("\nDone. Figures in: ", OUTPUT_DIR)
