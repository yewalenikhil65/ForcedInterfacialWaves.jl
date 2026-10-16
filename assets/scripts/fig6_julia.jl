using QuadGK, DelimitedFiles
using Plots, LaTeXStrings

struct CapillaryGravityParams
    α::Float64; ρᵣ::Float64; β::Float64; γρ::Float64
    atol::Float64; rtol::Float64
end

function makeCGParams()
    U, g, T = 26.7046, 981.0, 72.0
    ρₗ, ρᵤ = 1.0, 0.001
    ℓc = U^2 / g
    α = T / (ρₗ * U^2 * ℓc)
    ρᵣ = ρᵤ / ρₗ
    β = (1.0 - ρᵣ) / (1.0 + ρᵣ)
    return CapillaryGravityParams(α, ρᵣ, β, 1.0 / (1.0 + ρᵣ), 1e-10, 1e-8)
end

@inline χ(k::Float64, p::CapillaryGravityParams) =
    sqrt(p.β * k + p.γρ * p.α * k^3)

# Eqn. (4.5d): regularized I₄ integrand after cancellation of the pole factor.
function I₄Profile(xgrid::Vector{Float64}, t::Float64, p::CapillaryGravityParams)
    I₄ = Vector{Float64}(undef, length(xgrid))
    nchunks = min(Threads.nthreads(), length(xgrid))
    chunkLength = cld(length(xgrid), nchunks)

    Threads.@threads :static for chunk in 1:nchunks
        first = (chunk - 1) * chunkLength + 1
        last = min(chunk * chunkLength, length(xgrid))
        first > last && continue
        xchunk = @view xgrid[first:last]
        values = zeros(Float64, length(xchunk))
        integrand! = function (out, k)
            k == 0.0 && (fill!(out, 0.0); return out)
            χk = χ(k, p)
            amplitude = k / ((k + χk) * (1.0 + p.α * k^2 - p.ρᵣ))
            sinPhase, cosPhase = sincos(t * (k + χk))
            @inbounds @simd for i in eachindex(xchunk, out)
                sinKx, cosKx = sincos(k * xchunk[i])
                out[i] = amplitude * (cosPhase * cosKx + sinPhase * sinKx)
            end
            out
        end
        quadgk!(integrand!, values, 0.0, Inf; atol=p.atol, rtol=p.rtol,
                order=15, norm=values -> maximum(abs, values))
        copyto!(@view(I₄[first:last]), values)
    end
    return I₄
end

function main(; t::Float64=0.37, Nₓ::Int=2001)
    default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box,
            grid=false, guidefontsize=18, tickfontsize=16, legendfontsize=16)
    p = makeCGParams()
    xgrid = collect(range(-15.0, 15.0; length=Nₓ))
    filter!(x -> abs(x) > 1e-12, xgrid)
    I₄ = I₄Profile(xgrid, t, p)
    assets = @__DIR__
    writedlm(joinpath(assets, "julia_fig7_i4.csv"), hcat(xgrid, -I₄), ',')
    fig = plot(xgrid, -I₄; color="purple", label=L"\mathrm{Julia}",
               xlabel=L"x", ylabel=L"-\mathbb{I}_{4}",
               xlims=(-10, 10), ylims=(-2, 14))
    savefig(fig, joinpath(@__DIR__, "..", "cg_I4_fig7.png"))
    println("Saved Julia Figure 7 profile and cg_I4_fig7.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
