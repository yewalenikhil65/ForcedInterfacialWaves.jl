module WavesDelta

using QuadGK: quadgk, quadgk!, alloc_segbuf
using FresnelIntegrals: fresnelc, fresnels
using Printf

include("parameters.jl")
include("capillary_gravity.jl")
include("pure_gravity.jl")

# Parameters
export CapillaryGravityParams, PureGravityParams
export compute_cg_parameters, compute_gravity_parameters

# Capillary-gravity (Figure 10) — individual integrands for validation
export dispersion_chi, dispersion_omega
export cg_integrand_steady, cg_integrand_I3, cg_integrand_I4
export cg_combined_integrand, cg_partial_integrals
export cg_Gx_integral
export ivp_surface_elevation, steady_surface_elevation
export cg_profile_integrand!, cg_profile_partial_integrals
export compute_cg_ivp_profile, compute_cg_steady_profile
export compute_cg_profile, make_cg_xgrid

# Pure gravity (Figure 6) — individual T-terms for validation
export fresnel_C, fresnel_S
export gravity_T0
export gravity_T1_left, gravity_T2_left, gravity_T3_left, gravity_T4_left
export gravity_T1_right, gravity_T2_right, gravity_T3_right, gravity_T4_right
export gravity_analytical_left, gravity_analytical_right
export gravity_combined_integrand, gravity_numerical_cpv
export compute_gravity_profile, make_gravity_xgrid

end # module
