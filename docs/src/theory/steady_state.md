# Steady-State Decomposition (Manuscript §3.4)

This page derives the steady-state interfacial response from the general Fourier integral solution, establishing the far-field sinusoidal waves and localised deformation that form the basis for both the pure-gravity and capillary-gravity treatments that follow.

## Formal solution

A localised pressure $\tilde{p}_e = \tilde{F}_0\delta(\tilde{x})$ force is applied at the interface between two inviscid, incompressible, fluids of infinite depth, both streams moving at uniform speed $U$ rightwards. In the absence of this forcing, the interface is flat and remains at $\tilde{z}=0$. After non-dimensionalisation, the interfacial displacement resulting from the forcing i.e. $\eta(x,t)$ is obtained by evaluating the following time-dependent Fourier integrals (eqns. $3.7$ and $3.8$ in the manuscript)

```math
\dfrac{\sqrt{2\pi}\;\bar{\eta}(k,t)}{F_0} = \left(\dfrac{1}{-\alpha|k|^2 + \left(1+\rho_r\right)|k| - \left(1-\rho_r\right)}\right) - \dfrac{1}{\alpha |k|^2 + \left(1 - \rho_r\right)} \left(\dfrac{k}{2}\right) \Bigg( \dfrac{\exp\left[-it\lambda_2(k)\right]}{\lambda_2(k)} + \dfrac{\exp\left[-it\lambda_1(k)\right]}{\lambda_1(k)} \Bigg) \tag{3.7}
```
The inverse Fourier integral, leads to the (non-dimensional) interface displacement as a function of time:
```math
\eta(x,t) = \dfrac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}dk\;\exp\left(ikx\right)\bar{\eta}(k,t) \tag{3.8}
```
where $\lambda_{1,2}(k) \equiv k\mp\sqrt{\dfrac{\alpha|k|^3}{1+\rho_r}\;+\;\beta|k|}$ and

```math
\alpha = \frac{gT}{\rho_l U^4}, \quad
\rho_r = \frac{\rho_u}{\rho_l}, \quad
\beta = \frac{1-\rho_r}{1+\rho_r}.
```

## Steady-state decomposition

In this section, the manuscript shows that neglecting the time-dependent terms in eqns. (3.7) and (3.8) and $\textit{without}$ using any Rayleigh dissipation, the steady-state response turns out to be (we exclude all the Dirac delta function terms in eqn. $3.9$ in the manuscript) the following. For proof of this, see [Steady-state proof](../steady_proof.md).
```math
\dfrac{\eta_{s}(x)}{F_0} =\dfrac{ \eta^{\text{far-field}}_{s}(x)}{F_0} + \dfrac{\eta^{\text{local}}_{s}(x)}{F_0} \tag{3.11}
```
where ,
```math
\dfrac{\eta^{\text{far-field}}_{s}(x)}{F_0} \equiv \dfrac{1}{\alpha(k_l-k_s)}\bigg\{-\sin(k_s|x|)\quad + \quad \sin(k_l|x|)\bigg\}
```
```math
\dfrac{\eta^{\text{local}}_{s}(x)}{F_0} \equiv  \dfrac{\left(k_l+k_s\right)}{\pi\alpha}\int_{0}^{\infty}dy \dfrac{y\exp\left(-|x|y\right)}{\left(y^2 + k_l^2\right)\left(y^2+k_s^2\right)},\quad x\neq 0,
```

```math
k_{l,s} = \frac{1+\rho_r}{2\alpha}\left[1\pm\sqrt{1-\frac{4\alpha\beta}{1+\rho_r}}\right].
```

The integral expression for $\frac{\eta^{\text{local}}_{s}(x)}{F_0}$ above is solved numerically (in Julia and MATLAB) using the following codes. The Julia side uses [`compute_cg_parameters`](@ref) for the nondimensional parameters and the unicode `∫` alias (exported by the package) for the local-term quadrature.

```@example steady_state
using ForcedInterfacialWaves: ∫
using ForcedInterfacialWaves

p = compute_cg_parameters()

# Spatial grid (nondimensional), excluding x = 0
x = filter(xi -> abs(xi) > 1e-12, collect(-10:0.01:10))

# Far-field steady term
η_far = @. p.F₀ / (p.α * (p.kₗ - p.kₛ)) *
           (-sin(p.kₛ * abs(x)) + sin(p.kₗ * abs(x)))

# Local steady term via direct quadrature
η_local = [p.F₀ * (p.kₗ + p.kₛ) / (π * p.α) *
           ∫(y -> y * exp(-abs(xi) * y) / ((y^2 + p.kₗ^2) * (y^2 + p.kₛ^2)),
             0.0, Inf; atol=p.atol, rtol=p.rtol)
           for xi in x]

extrema(η_far), extrema(η_local)
```

The same components can be plotted directly. Because Plots.jl is available in the documentation
environment, this block is executed and rendered by Documenter:

```@example steady_state
using Plots
plot(x, η_far .* 1e3; label="η_s far-field", lw=1.5)
plot!(x, η_local .* 1e3; label="η_s local", lw=1.5,
      xlabel="x", ylabel="η × 10³")
```

```matlab
%% Capillary-gravity steady-state decomposition
U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0;  rho_u = 0.001;

% Characteristic scale and nondimensional roots
l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;

discriminant = (1.0 + rho_r)^2 - 4.0 * alpha * (1.0 - rho_r);
k_l = ((1.0 + rho_r) + sqrt(discriminant)) / (2.0 * alpha);
k_s = ((1.0 + rho_r) - sqrt(discriminant)) / (2.0 * alpha);

F0 = 0.01 * T / (rho_l * U^2 * l_c);

% Spatial grid (nondimensional), excluding x = 0
x = -10:0.01:10;
x(x == 0) = [];

% Far-field steady term
eta_far = F0 / (alpha * (k_l - k_s)) .* ...
    (-sin(k_s * abs(x)) + sin(k_l * abs(x)));

% Local steady term via direct quadrature
eta_local = zeros(size(x));

for i = 1:length(x)
    integrand = @(y) y .* exp(-abs(x(i)) .* y) ./ ...
        ((y.^2 + k_l^2) .* (y.^2 + k_s^2));

    I = integral(integrand, 0, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-8);

    eta_local(i) = F0 * (k_l + k_s) / (pi * alpha) * I;
end

% Plot both components (x1e3 for readability)
figure;
hold on;
plot(x, eta_far * 1e3, 'b-', 'LineWidth', 1.5, ...
    'DisplayName', '\eta_s far-field');
plot(x, eta_local * 1e3, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', '\eta_s local');
xlabel('x');
ylabel('\eta \times 10^3');
legend('Location', 'best');
grid on;
```

The resulting figure from both Julia and MATLAB are superimposed  in the following figure, which is a reproduction of fig. $5a$ in the manuscript.

```@raw html
<figure style="text-align:center;">
  <img src="../assets/Fig5a.png" alt="Figure 5a comparison" style="max-width:60%; height:auto;">
  <figcaption>Figure 5a: Comparison of steady wave profiles computed in MATLAB vs Julia.</figcaption>
</figure>
```
The steady-state response $\eta(x)$ using the Rayleigh dissipation approach is obtained as:
```math
\begin{aligned}
\dfrac{\eta(x)}{F_0} &= -\dfrac{2}{\alpha\left(k_{l} - k_{s}\right)}
\begin{cases}
    \sin(k_{l}x), & x<0\\
    \sin(k_{s}x), & x>0
\end{cases} + \dfrac{G(x)}{\pi\alpha}, \\
\text{with}\quad G(x) &\equiv \dfrac{1}{k_{l}-k_{s}}\int_{0}^{\infty}\;dk\;\left(\dfrac{\cos(kx)}{k+k_{s}}-\dfrac{\cos(kx)}{k+k_{l}}\right).
\end{aligned} \tag{3.12}
```

It is shown [here](../capillary_gravity_rayleigh_dissipation.md) that $\dfrac{G(x)}{\pi\alpha}$ in eqn. 3.12  is identical to $\dfrac{\eta_s^{\text{local}}(x)}{F_0}$ in eqn. 3.11; the latter expression being preferable compared to $G(x)$ in 3.12 due to the apparence of its local nature via the explicit exponential decay term. Fig. 5b confirms that the steady-state response from the Rayleigh dissipation approach is qualitatively different from that of fig. 5a. Notably, the response employing Rayleigh dissipation is asymmetric about $x=0$ (fig. 5b), consistent with observations. In the next sections, we show that similar results will be arrived at through the IVP approach.
