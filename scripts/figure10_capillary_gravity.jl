"""
    figure10_capillary_gravity.jl

    Capillary-gravity (α > 0) IVP solution vs steady state.
    Reproduces Figure 10 of the manuscript at multiple time snapshots.

    Usage:
        julia -t auto --project=scripts scripts/figure10_capillary_gravity.jl
"""

# ─── Thread check ──────────────────────────────────────────────────────────────
if Threads.nthreads() == 1
    @warn "Running single-threaded. For ~5× speedup: julia -t auto --project=scripts scripts/figure10_capillary_gravity.jl"
else
    @info "Using $(Threads.nthreads()) threads"
end

using ForcedInterfacialWaves
using Plots; gr()

# ═══════════════════════════════════════════════════════════════════════════════
# Parameters
# ═══════════════════════════════════════════════════════════════════════════════

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "figure10")
const COMPARISON_FILE = joinpath(
    @__DIR__, "..", "docs", "src", "assets", "matlab_cg_ivp_t300.csv"
)
mkpath(OUTPUT_DIR)

const TIME_INDICES = [1, 3, 7, 15, 25, 60, 145, 300]
const NX_PLOT      = 2001

"""Read a headerless two-column `(x, η)` CSV without adding a package dependency."""
function read_comparison_csv(path)
    rows = [split(strip(line), ',') for line in eachline(path) if !isempty(strip(line))]
    all(length(row) >= 2 for row in rows) || error("Expected at least two columns in $path")
    x = [parse(Float64, row[1]) for row in rows]
    η = [parse(Float64, row[2]) for row in rows]
    return x, η
end

comparison = if isfile(COMPARISON_FILE)
    println("CSV comparison: ", COMPARISON_FILE)
    read_comparison_csv(COMPARISON_FILE)
else
    @warn "CSV comparison not found; figures will contain theory only" path=COMPARISON_FILE
    nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Setup
# ═══════════════════════════════════════════════════════════════════════════════

p = compute_cg_parameters()

println("Capillary-gravity parameters:")
println("  α   = ", p.α)
println("  ρᵣ  = ", p.ρᵣ)
println("  kₛ  = ", p.kₛ)
println("  kₗ  = ", p.kₗ)
println("  F₀  = ", p.F₀)
println("  λ_capillary = ", 2π / (p.kₗ / p.l_c), " cm")
println("  λ_gravity   = ", 2π / (p.kₛ / p.l_c), " cm")

x_grid = make_cg_xgrid(p; Nx=NX_PLOT)

# Explicit asymmetric classical steady state for long-time comparison.
# sol.η_steady itself is the symmetric ηₛ from equation (4.5b).
classical_steady = solve(ForcedGCProblem(p, x_grid);
                         method=steady(rayleigh_dissipation=true)).η

# ═══════════════════════════════════════════════════════════════════════════════
# Loop over times
# ═══════════════════════════════════════════════════════════════════════════════

for ti in TIME_INDICES
    t_dim = ti / 100.0
    t = t_dim / p.t_c

    print("\n  time_index = $ti, t_dim = $(t_dim) s, t = $t ... ")
    flush(stdout)

    t0 = time()
    sol = solve(ForcedGCProblem(p, x_grid, t))
    elapsed = time() - t0
    println("done in $(round(elapsed; digits=2)) s")

    η_ivp    = sol.η
    η_s      = sol.η_steady
    η_tr     = sol.η_transient
    scale    = 1.0e3

    # ─── Plot ───
    plt = plot(;
        size       = (500, 340),
        xlims      = (-5.0, 10.0),
        ylims      = (-6.0, 13.0),
        yticks     = [-4, 0, 4, 8, 12],
        xlabel     = "\$x\$",
        ylabel     = "\$10^3 \\eta\$",
        legend     = :topright,
        framestyle = :box,
        grid       = false,
        tickdir    = :out,
    )

    plot!(plt, x_grid, scale .* η_ivp;
          lw=1.8, ls=:dash, lc=:blue, label="\$\\eta\$ (IVP)")
    plot!(plt, x_grid, scale .* η_s;
          lw=1.8, ls=:solid, lc=:black, label="\$\\eta_s\$ (eq. 4.5b)")
    plot!(plt, x_grid, scale .* η_tr;
          lw=1.4, ls=:dot, lc=:green, label="\$\\eta_{\\rm tr}\$")

    if ti == 300
        # Also show the asymmetric long-time radiation solution and archived
        # MATLAB reference; these are not sol.η_steady.
        plot!(plt, x_grid, scale .* classical_steady;
              lw=1.8, ls=:dashdot, lc=:gray, label="classical radiation")
        if comparison !== nothing
            x_csv, η_csv = comparison
            plot!(plt, x_csv, scale .* η_csv;
                  lw=1.6, ls=:dot, lc=:red, label="CSV reference")
        end
    end

    # ─── Save ───
    fname = "figure10_t" * lpad(ti, 4, '0')
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".png"))
    savefig(plt, joinpath(OUTPUT_DIR, fname * ".pdf"))
    println("    saved: ", fname)
end

println("\nDone. Figures in: ", OUTPUT_DIR)
