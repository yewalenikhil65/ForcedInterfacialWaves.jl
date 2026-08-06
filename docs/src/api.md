# API Reference

## Parameters

```@docs
CapillaryGravityParams
PureGravityParams
compute_cg_parameters
compute_gravity_parameters
```

## Pure Gravity (Figure 6)

### Fresnel integrals

```@docs
fresnel_C
fresnel_S
```

### Analytical terms

```@docs
gravity_T0
gravity_T1_left
gravity_T2_left
gravity_T3_left
gravity_T4_left
gravity_T1_right
gravity_T2_right
gravity_T3_right
gravity_T4_right
gravity_analytical_left
gravity_analytical_right
```

### Numerical CPV

```@docs
gravity_combined_integrand
gravity_numerical_cpv
```

### Profiles

```@docs
compute_gravity_profile
make_gravity_xgrid
```

## Capillary–Gravity (Figure 10)

### Dispersion and integrands

```@docs
dispersion_chi
cg_integrand_steady
cg_integrand_I3
cg_integrand_I4
cg_combined_integrand
```

### Quadrature and profiles

```@docs
cg_partial_integrals
ivp_surface_elevation
cg_Gx_integral
steady_surface_elevation
compute_cg_ivp_profile
compute_cg_steady_profile
compute_cg_profile
make_cg_xgrid
```
