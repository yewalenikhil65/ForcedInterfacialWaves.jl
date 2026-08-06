"""
    Parameters and nondimensionalization for interfacial wave IVP.

    Dimensional parameters are in CGS units:
    - U: base flow speed [cm/s]
    - g: gravitational acceleration [cm/s²]
    - T: surface tension [dyn/cm]
    - rho_l: lower fluid density [g/cm³]
    - rho_u: upper fluid density [g/cm³]
"""

"""
    CapillaryGravityParams

All nondimensional parameters needed for the capillary-gravity (α > 0) IVP.
"""
struct CapillaryGravityParams
    # Dimensional
    U::Float64
    g::Float64
    T::Float64
    rho_l::Float64
    rho_u::Float64
    L::Float64          # Domain length [cm]

    # Characteristic scales
    l_c::Float64
    t_c::Float64
    F_c::Float64

    # Nondimensional
    alpha::Float64
    rho_r::Float64
    beta::Float64       # Atwood number
    gamma_rho::Float64

    # Wave roots
    k_l::Float64        # Capillary root (large k)
    k_s::Float64        # Gravity root (small k)

    # Forcing
    F0::Float64

    # Quadrature parameters
    epsilon_pv::Float64
    atol::Float64
    rtol::Float64
    k_max::Float64
    k_max_steady::Float64
end

"""
    PureGravityParams

All nondimensional parameters needed for the pure-gravity (α = 0) IVP.
"""
struct PureGravityParams
    # Dimensional
    U::Float64
    g::Float64
    rho_l::Float64
    rho_u::Float64

    # Characteristic scales
    l_c::Float64
    t_c::Float64
    F_c::Float64

    # Nondimensional
    rho_r::Float64
    beta::Float64
    sqrt_beta::Float64

    # Forcing
    F0::Float64

    # Domain
    gravity_wavelength::Float64
    L::Float64

    # Quadrature
    epsilon_cpv::Float64
    atol_cpv::Float64
    rtol_cpv::Float64
    atol_transformed::Float64
    rtol_transformed::Float64
    k_max_analytical::Float64
    k_max_cpv::Float64

    # Front exclusion band
    front_band::Float64
end

"""
    compute_cg_parameters(; U, g, T, rho_l, rho_u, L, F0_factor, epsilon_pv, atol, rtol, k_max, k_max_steady)

Compute all capillary-gravity parameters from dimensional inputs.
"""
function compute_cg_parameters(;
    U::Float64 = 26.7046,
    g::Float64 = 981.0,
    T::Float64 = 72.0,
    rho_l::Float64 = 1.0,
    rho_u::Float64 = 0.001,
    L::Float64 = 75.5996,
    F0_factor::Float64 = 0.01,
    epsilon_pv::Float64 = 1.0e-6,
    atol::Float64 = 1.0e-10,
    rtol::Float64 = 1.0e-8,
    k_max::Float64 = Inf,
    k_max_steady::Float64 = Inf
)
    # Characteristic scales
    l_c = U^2 / g
    t_c = U / g
    F_c = rho_l * U^2 * l_c

    # Nondimensional parameters
    alpha = T / (rho_l * U^2 * l_c)
    rho_r = rho_u / rho_l
    beta = (1.0 - rho_r) / (1.0 + rho_r)
    gamma_rho = 1.0 / (1.0 + rho_r)

    # Gravity and capillary wave roots
    discriminant = (1.0 + rho_r)^2 - 4.0 * alpha * (1.0 - rho_r)
    if discriminant <= 0.0
        error("The steady capillary-gravity roots are not distinct positive real numbers. discriminant = $discriminant")
    end

    k_l = ((1.0 + rho_r) + sqrt(discriminant)) / (2.0 * alpha)
    k_s = ((1.0 + rho_r) - sqrt(discriminant)) / (2.0 * alpha)

    # Forcing amplitude
    F0_dim = F0_factor * T
    F0 = F0_dim / F_c

    return CapillaryGravityParams(
        U, g, T, rho_l, rho_u, L,
        l_c, t_c, F_c,
        alpha, rho_r, beta, gamma_rho,
        k_l, k_s,
        F0,
        epsilon_pv, atol, rtol, k_max, k_max_steady
    )
end

"""
    compute_gravity_parameters(; U, g, rho_l, rho_u, sigma_reference, F0_factor, ...)

Compute all pure-gravity parameters from dimensional inputs.
"""
function compute_gravity_parameters(;
    U::Float64 = 26.7046,
    g::Float64 = 981.0,
    rho_l::Float64 = 1.0,
    rho_u::Float64 = 0.001,
    sigma_reference::Float64 = 72.0,
    F0_factor::Float64 = 0.01,
    epsilon_cpv::Float64 = 1.0e-6,
    atol_cpv::Float64 = 1.0e-10,
    rtol_cpv::Float64 = 1.0e-8,
    atol_transformed::Float64 = 1.0e-10,
    rtol_transformed::Float64 = 1.0e-8,
    k_max_analytical::Float64 = 100.0,
    k_max_cpv::Float64 = 100.0,
    front_band::Float64 = 1.0
)
    # Characteristic scales
    l_c = U^2 / g
    t_c = U / g
    F_c = rho_l * U^2 * l_c

    # Nondimensional parameters
    rho_r = rho_u / rho_l
    beta = (1.0 - rho_r) / (1.0 + rho_r)
    sqrt_beta = sqrt(beta)

    # Forcing amplitude
    F_nd = sigma_reference / F_c
    F0 = F0_factor * F_nd

    # Domain
    gravity_wavelength = 2.0 * pi / beta
    L = 4.0 * gravity_wavelength

    return PureGravityParams(
        U, g, rho_l, rho_u,
        l_c, t_c, F_c,
        rho_r, beta, sqrt_beta,
        F0,
        gravity_wavelength, L,
        epsilon_cpv, atol_cpv, rtol_cpv,
        atol_transformed, rtol_transformed,
        k_max_analytical, k_max_cpv,
        front_band
    )
end
