using QuadGK

"""Typed nondimensional parameters for the two-fluid capillary–gravity IVP."""
struct CapillaryGravityParams
    α::Float64; ρᵣ::Float64; β::Float64; γρ::Float64
    kₗ::Float64; kₛ::Float64; F₀::Float64
    εCPV::Float64; atol::Float64; rtol::Float64
end

function makeCGParams()
    U, g, T = 26.7046, 981.0, 72.0
    ρₗ, ρᵤ = 1.0, 0.001
    ℓc = U^2 / g
    α = T / (ρₗ * U^2 * ℓc)
    ρᵣ = ρᵤ / ρₗ
    β = (1.0 - ρᵣ) / (1.0 + ρᵣ)
    Δ = (1.0 + ρᵣ)^2 - 4.0 * α * (1.0 - ρᵣ)
    kₗ = ((1.0 + ρᵣ) + sqrt(Δ)) / (2.0 * α)
    kₛ = ((1.0 + ρᵣ) - sqrt(Δ)) / (2.0 * α)
    F₀ = 0.01 * T / (ρₗ * U^2 * ℓc)
    return CapillaryGravityParams(α, ρᵣ, β, 1.0/(1.0 + ρᵣ),
                                  kₗ, kₛ, F₀, 1e-6, 1e-10, 1e-8)
end

# Dispersion relation χ(k) = √(βk + γρ αk³).
@inline χ(k::Float64, p::CapillaryGravityParams) =
    sqrt(p.β * k + p.γρ * p.α * k^3)

# Combined integrand in eqn. (4.5a); the kₛ and kₗ pole terms cancel in the sum.
function combinedIntegrand(k::Float64, x::Float64, t::Float64, p::CapillaryGravityParams)
    χk = χ(k, p)
    pole = 1.0 / (p.α * (k - p.kₗ) * (k - p.kₛ))
    dispersion = (1.0 + p.ρᵣ) / (1.0 + p.α * k^2 - p.ρᵣ)
    phase = k * (t - x)
    return pole * (2.0 * cos(k * x) -
           dispersion * (k + χk) * cos(phase - t * χk) -
           dispersion * (k - χk) * cos(phase + t * χk))
end

# CPV quadrature split about kₛ and kₗ.
function cpvParts(x::Float64, t::Float64, p::CapillaryGravityParams)
    f = k -> combinedIntegrand(k, x, t, p)
    I₁, _ = quadgk(f, 0.0, p.kₛ - p.εCPV; atol=p.atol, rtol=p.rtol)
    I₂, _ = quadgk(f, p.kₛ + p.εCPV, p.kₗ - p.εCPV; atol=p.atol, rtol=p.rtol)
    I₃, _ = quadgk(f, p.kₗ + p.εCPV, Inf; atol=p.atol, rtol=p.rtol, order=15)
    return I₁, I₂, I₃
end

# Eqn. (4.7): local steady integral G(x).
function Gₓ(x::Float64, p::CapillaryGravityParams)
    f = k -> (cos(k*x)/(k + p.kₛ) - cos(k*x)/(k + p.kₗ)) / (p.kₗ - p.kₛ)
    return first(quadgk(f, 0.0, Inf; atol=p.atol, rtol=p.rtol))
end

p = makeCGParams()
x, t = 3.0, 110.0
I₁, I₂, I₃ = cpvParts(x, t, p)
η = -p.F₀ * (I₁ + I₂ + I₃) / (2π)
G = Gₓ(x, p)
ηₛ = p.F₀/(p.α*(p.kₗ-p.kₛ)) * (-sin(p.kₛ*abs(x)) + sin(p.kₗ*abs(x))) + p.F₀*G/(π*p.α)
ηₜᵣ = η - ηₛ
ηclassical = p.F₀ * (-2sin(p.kₛ*x)/(p.α*(p.kₗ-p.kₛ)) + G/(π*p.α))

println("𝕀(k=2; x=3, t=110) = ", combinedIntegrand(2.0, x, t, p))
println("I₁ = ", I₁)
println("I₂ = ", I₂)
println("I₃ = ", I₃)
println("η = ", η)
println("ηₛ = ", ηₛ)
println("ηₜᵣ = ", ηₜᵣ)
println("ηclassical = ", ηclassical)
println("Gₓ = ", G)
using DelimitedFiles

using QuadGK, Plots, LaTeXStrings

default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box,
        grid=false, guidefontsize=18, tickfontsize=16, legendfontsize=16)

# Full CG IVP profile: combined CPV integrand (eqn 4.5a–d, 3 pole-split quadratures)
# and G(x) for η_s (eqn 4.5b) are computed per spatial chunk across threads.
function cgIVPProfile(xgrid::Vector{Float64}, t::Float64, p::CapillaryGravityParams)
    α, ρᵣ, F₀ = p.α, p.ρᵣ, p.F₀
    kₗ, kₛ, ε = p.kₗ, p.kₛ, p.εCPV
    atol, rtol = p.atol, p.rtol
    x_grid = xgrid

    N = length(x_grid)
    η = Vector{Float64}(undef, N)
    G = Vector{Float64}(undef, N)  # G(x) for η_s

    nchunks = min(Threads.nthreads(), N); clen = cld(N, nchunks)
    Threads.@threads :static for ci in 1:nchunks
        lo = (ci-1)*clen + 1; hi = min(ci*clen, N); lo > hi && continue
        xc = view(x_grid, lo:hi); M = length(xc)

        # IVP η: combined integrand 𝕀(k;x,t), split at kₛ±ε and kₗ±ε
        I1 = zeros(M); I2 = zeros(M); I3 = zeros(M)
        function ivp!(vals, k)
            c = χ(k, p); invp = 1/(α*(k-kₗ)*(k-kₛ)); invd = (1+ρᵣ)/(1+α*k^2-ρᵣ)
            am = invd*(k+c); ap = invd*(k-c)
            sm, cm = sincos(t*(k-c)); sp, cp = sincos(t*(k+c))
            cc = 2*invp - (am*cm + ap*cp)*invp
            sc =        - (am*sm + ap*sp)*invp
            @inbounds @simd for i in 1:M
                s_kx, c_kx = sincos(k*xc[i]); vals[i] = cc*c_kx + sc*s_kx
            end; vals
        end
        quadgk!(ivp!, I1, 0.0, kₛ-ε; atol=atol, rtol=rtol, norm=v->maximum(abs,v))
        quadgk!(ivp!, I2, kₛ+ε, kₗ-ε; atol=atol, rtol=rtol, norm=v->maximum(abs,v))
        quadgk!(ivp!, I3, kₗ+ε, Inf;   atol=atol, rtol=rtol, order=15, norm=v->maximum(abs,v))
        @inbounds @simd for i in 1:M; η[lo+i-1] = -F₀/(2π)*(I1[i]+I2[i]+I3[i]); end

        # Steady G(x): Lamb integrand (kₗ−kₛ in denominator already cancelled)
        Gc = zeros(M)
        function gint!(vals, k)
            coeff = (1/(k+kₛ) - 1/(k+kₗ)) / (kₗ - kₛ)
            @inbounds @simd for i in 1:M; vals[i] = coeff*cos(k*xc[i]); end; vals
        end
        quadgk!(gint!, Gc, 0.0, Inf; atol=atol, rtol=rtol, norm=v->maximum(abs,v))
        copyto!(view(G, lo:hi), Gc)
    end

    η_s = similar(η)
    @inbounds @simd for i in 1:N
        xv = x_grid[i]
        η_s[i] = F₀/(α*(kₗ-kₛ)) * (-sin(kₛ*abs(xv)) + sin(kₗ*abs(xv))) + F₀*G[i]/(π*α)
    end
    return η, η_s, η .- η_s
end

x_grid = collect(range(-15.0, 15.0; length=2001))
filter!(x -> abs(x) > 1e-12, x_grid)
t_fig8 = 367.35
p = makeCGParams()
η, η_s, η_tr = cgIVPProfile(x_grid, t_fig8, p)

plot(x_grid, η .* 1e3; label=L"\eta", color="blue", ls=:dash,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}",
     xlims=(-10,10), ylims=(-4.8, 8.2), yticks=[-4, 0, 4, 8], legend=:outerright)
plot!(x_grid, η_tr .* 1e3; label=L"\eta_{tr}", color="magenta", ls=:dot)
writedlm(joinpath(@__DIR__, "julia_fig8_profiles.csv"), hcat(x_grid, η, η_tr), ',')
savefig(current(), joinpath(@__DIR__, "..", "cg_ivp_fig8.png"))
