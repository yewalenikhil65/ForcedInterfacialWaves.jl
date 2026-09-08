"""
    figure6_pure_gravity.jl

    Julia equivalent of the MATLAB Fig6 plotting script.
    Pure-gravity (α = 0) IVP: analytical vs numerical CPV evaluation.

    Usage:
        julia -t auto --project=. scripts/figure6_pure_gravity.jl
"""

# Thread check
if Threads.nthreads() == 1
    @warn "Running single-threaded. For speedup, restart with: julia -t auto --project=. scripts/figure6_pure_gravity.jl"
else
    @info "Using $(Threads.nthreads()) threads"
end

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using ForcedInterfacialWaves
using Plots; gr()
using Printf

# ═══════════════════════════════════════════════════════════════════════════════
# Parameters
# ═══════════════════════════════════════════════════════════════════════════════

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "figure6")
mkpath(OUTPUT_DIR)

const TIME_INDICES    = [1, 7, 15, 40, 60, 100, 300, 500]
const NX_PLOT         = 2001
const MARKER_SPACING  = 50

# ═══════════════════════════════════════════════════════════════════════════════
# Setup
# ═══════════════════════════════════════════════════════════════════════════════

pg = compute_gravity_parameters()

@printf("Pure-gravity parameters:\n")
@printf("  β         = %.15e\n", pg.beta)
@printf("  √β        = %.15e\n", pg.sqrt_beta)
@printf("  F₀        = %.6e\n", pg.F0)
@printf("  λ_gravity = %.4f\n", pg.gravity_wavelength)
@printf("  L         = %.4f\n", pg.L)

x_grid = make_gravity_xgrid(pg; Nx=NX_PLOT)
half_L = pg.L / 2.0

# ═══════════════════════════════════════════════════════════════════════════════
# Loop over times
# ═══════════════════════════════════════════════════════════════════════════════

for ti in TIME_INDICES
    t_dim = ti / 100.0
    t = t_dim / pg.t_c

    @printf("\n  time_index = %d, t_dim = %.2f s, t = %.4f ...\n", ti, t_dim, t)

    # Compute profiles
    η_analytical, η_steady, η_tr, η_cpv = compute_gravity_profile(x_grid, t, pg)

    # Error metrics (outside front band only)
    valid = .!isnan.(η_analytical) .& .!isnan.(η_cpv)
    if any(valid)
        abs_err = abs.(η_analytical[valid] .- η_cpv[valid])
        max_err = maximum(abs_err)
        @printf("    max|η_ana − η_CPV| = %.4e\n", max_err)
    end

    # ─── Masks for left / right of front ───
    mask_left  = (t .- x_grid) .>  pg.front_band
    mask_right = (t .- x_grid) .< -pg.front_band

    # Scale for plotting
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
        fontfamily = "Computer Modern",
    )

    # 1. Analytical total (blue dash-dot)
    idx_l = findall(mask_left)
    idx_r = findall(mask_right)

    if !isempty(idx_l)
        plot!(plt, x_grid[idx_l], scale .* η_analytical[idx_l];
              lw=2.0, ls=:dashdot, lc=:blue, label="\$\\eta\$")
    end
    if !isempty(idx_r)
        plot!(plt, x_grid[idx_r], scale .* η_analytical[idx_r];
              lw=2.0, ls=:dashdot, lc=:blue, label="")
    end

    # 2. Steady (black solid)
    plot!(plt, x_grid, scale .* η_steady;
          lw=1.5, ls=:solid, lc=:black, label="\$\\eta_s\$")

    # 3. Transient (magenta dotted)
    if !isempty(idx_l)
        plot!(plt, x_grid[idx_l], scale .* η_tr[idx_l];
              lw=1.5, ls=:dot, lc=:magenta, label="\$\\eta_{\\mathrm{tr}}\$")
    end
    if !isempty(idx_r)
        plot!(plt, x_grid[idx_r], scale .* η_tr[idx_r];
              lw=1.5, ls=:dot, lc=:magenta, label="")
    end

    # 4. Numerical CPV (purple circle markers, sparse)
    cpv_color = RGB(0.49, 0.18, 0.56)
    if !isempty(idx_l)
        marker_idx_l = idx_l[1:MARKER_SPACING:end]
        scatter!(plt, x_grid[marker_idx_l], scale .* η_cpv[marker_idx_l];
                 mc=cpv_color, ms=3, msw=0.8, msc=cpv_color, shape=:circle,
                 markerstrokecolor=cpv_color, label="\$\\eta_{\\mathrm{CPV}}\$")
    end
    if !isempty(idx_r)
        marker_idx_r = idx_r[1:MARKER_SPACING:end]
        scatter!(plt, x_grid[marker_idx_r], scale .* η_cpv[marker_idx_r];
                 mc=cpv_color, ms=3, msw=0.8, msc=cpv_color, shape=:circle,
                 markerstrokecolor=cpv_color, label="")
    end

    # 5. Vertical dashed line at x = t
    vline!(plt, [t]; ls=:dash, lc=:black, lw=1.2, label="")

    # Save
    fname = @sprintf("figure6_t%04d", ti)
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".png"))
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".pdf"))
    @printf("    saved: %s\n", fname)
end

println("\nDone. Figures in: ", OUTPUT_DIR)
