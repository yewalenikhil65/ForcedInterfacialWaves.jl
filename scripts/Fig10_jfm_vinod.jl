"""
    figure10_capillary_gravity.jl

    Julia equivalent of the MATLAB Fig10 plotting script.
    Capillary-gravity (α > 0) IVP solution vs steady state.

    Usage:
        julia -t auto --project=. scripts/figure10_capillary_gravity.jl
"""

# Thread check
if Threads.nthreads() == 1
    @warn "Running single-threaded. For ~5× speedup, restart with: julia -t auto --project=. scripts/figure10_capillary_gravity.jl"
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

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "figure10")
mkpath(OUTPUT_DIR)

const TIME_INDICES = [1, 3, 7, 15, 25, 60, 145, 300]
const NX_PLOT      = 2001

# ═══════════════════════════════════════════════════════════════════════════════
# Setup
# ═══════════════════════════════════════════════════════════════════════════════

p = compute_cg_parameters()  # MATLAB-exact: k_max=Inf, atol=1e-10, rtol=1e-8

@printf("Capillary-gravity parameters:\n")
@printf("  α     = %.6e\n", p.alpha)
@printf("  ρ_r   = %.6e\n", p.rho_r)
@printf("  k_s   = %.12f\n", p.k_s)
@printf("  k_l   = %.12f\n", p.k_l)
@printf("  F₀    = %.6e\n", p.F0)

k_l_dim = p.k_l / p.l_c
k_s_dim = p.k_s / p.l_c
@printf("  λ_capillary = %.4f cm\n", 2π / k_l_dim)
@printf("  λ_gravity   = %.4f cm\n", 2π / k_s_dim)

x_grid = make_cg_xgrid(p; Nx=NX_PLOT)

# ═══════════════════════════════════════════════════════════════════════════════
# Loop over times — threaded in-place vector quadrature over spatial chunks
# ═══════════════════════════════════════════════════════════════════════════════

for ti in TIME_INDICES
    t_dim = ti / 100.0
    t = t_dim / p.t_c

    @printf("\n  time_index = %d, t_dim = %.2f s, t = %.4f ... ", ti, t_dim, t)
    flush(stdout)

    t0 = time()
    η_ivp, η_steady = compute_cg_profile(
        x_grid, t, p; compute_steady=(ti == 300))
    elapsed = time() - t0
    @printf("done in %.1f s\n", elapsed)

    # Scale for plotting
    scale = 1.0e3

    # ─── Plot ───
    plt = plot(;
        size       = (500, 340),
        xlims      = (-5.0, 10.0),
        ylims      = (-6.0, 13.0),
        yticks     = [-4, 0, 4, 8],
        xlabel     = "\$x\$",
        ylabel     = "\$10^3 \\eta\$",
        legend     = :topright,
        framestyle = :box,
        grid       = false,
        tickdir    = :out,
        fontfamily = "Computer Modern",
    )

    if ti == 300
        # Steady (black solid)
        plot!(plt, x_grid, scale .* η_steady;
              lw=1.8, ls=:solid, lc=:black, label="Steady")
        # IVP theory (blue dashed)
        plot!(plt, x_grid, scale .* η_ivp;
              lw=1.8, ls=:dash, lc=:blue, label="Theory")
    else
        # IVP theory (blue dash-dot)
        plot!(plt, x_grid, scale .* η_ivp;
              lw=1.8, ls=:dashdot, lc=:blue, label="Theory")
    end

    # Save
    fname = @sprintf("figure10_t%04d", ti)
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".png"))
    # savefig(plt, joinpath(OUTPUT_DIR, fname * ".pdf"))
    # @printf("    saved: %s\n", fname)
end

println("\nDone. Figures in: ", OUTPUT_DIR)
