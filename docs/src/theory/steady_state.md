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

The integral expression for $\frac{\eta^{\text{local}}_{s}(x)}{F_0}$ above is solved numerically (in Julia and MATLAB) using the following codes. The Julia code writes out equation (3.11) term by term — exactly as the MATLAB reference does — so each line maps directly onto the mathematics: the far-field is the closed-form sine combination, and the local part is the exponentially-decaying integral evaluated with adaptive quadrature (`quadgk`).

```julia
using QuadGK          # adaptive Gauss–Kronrod quadrature (∫₀^∞ …)
using Plots, LaTeXStrings

# ─── Nondimensional parameters (see the Overview page) ───
U, g, T = 26.7046, 981.0, 72.0
ρₗ, ρᵤ  = 1.0, 0.001
l_c = U^2 / g
α   = T / (ρₗ * U^2 * l_c)
ρᵣ  = ρᵤ / ρₗ
F₀  = 0.01 * T / (ρₗ * U^2 * l_c)

# Steady capillary/gravity roots  k_{l,s}  (roots of  α k² − (1+ρᵣ)k + (1−ρᵣ) = 0)
Δ  = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
kₗ = ((1 + ρᵣ) + √Δ) / (2α)      # long-wavelength (capillary) root
kₛ = ((1 + ρᵣ) - √Δ) / (2α)      # short-wavelength (gravity) root

# ─── Equation (3.11):  ηₛ(x)/F₀ = far-field + local ───

# Far-field: closed-form sine combination
η_farfield(x) = F₀ / (α * (kₗ - kₛ)) * (-sin(kₛ * abs(x)) + sin(kₗ * abs(x)))

# Local: exponentially-decaying integral  (kₗ+kₛ)/(πα) ∫₀^∞ y e^{-|x|y} / [(y²+kₗ²)(y²+kₛ²)] dy
# quadgk returns (value, error); `first` keeps the value.
η_local(x) = F₀ * (kₗ + kₛ) / (π * α) *
             first(quadgk(y -> y * exp(-abs(x) * y) / ((y^2 + kₗ^2) * (y^2 + kₛ^2)),  0, Inf; atol = 1e-10, rtol = 1e-8))

# ─── Evaluate on a grid (skip x = 0, where the local integrand is singular) ───
x = range(-15.0, 15.0; length = 2001)
x = filter(xi -> abs(xi) > 1e-12, collect(x))

η_ff  = η_farfield.(x)
η_loc = η_local.(x)

plot(x, η_loc .* 1e3; label="local", color="red",
     guidefontsize=16, tickfontsize=14, legendfontsize=14,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
     legend=:outerright, size=(800,400))
plot!(x, η_ff .* 1e3; label="far-field", color="blue")
```

```matlab
%% Steady-state decomposition (no Rayleigh dissipation) — eqn. (3.11)
U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0; rho_u = 0.001;
l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
disc  = (1 + rho_r)^2 - 4*alpha*(1 - rho_r);
k_l   = ((1 + rho_r) + sqrt(disc)) / (2*alpha);   % long-wavelength (capillary) root
k_s   = ((1 + rho_r) - sqrt(disc)) / (2*alpha);   % short-wavelength (gravity) root
F0    = 0.01*T / (rho_l*U^2*l_c);

x = linspace(-15, 15, 2001);
x(abs(x) < 1e-12) = [];

% Far-field: closed-form sine combination
eta_far = F0/(alpha*(k_l - k_s)) .* (-sin(k_s*abs(x)) + sin(k_l*abs(x)));

% Local: exponentially-decaying integral
eta_local = zeros(size(x));
for i = 1:length(x)
    eta_local(i) = F0*(k_l + k_s)/(pi*alpha) * ...
        integral(@(y) y.*exp(-abs(x(i)).*y)./((y.^2+k_l^2).*(y.^2+k_s^2)), 0, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-8);
end

figure; hold on;
plot(x, eta_local*1e3, 'r-', 'LineWidth', 2);
plot(x, eta_far*1e3, 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('\eta \times 10^3');
legend('local','far-field'); xlim([-10 10]);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/steady_no_rayleigh.png" alt="Fig 5a(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 5a(i).</strong> Steady-state far-field and local components without Rayleigh dissipation (Julia).</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig5a_overlay.png" alt="Fig 5a(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 5a(ii).</strong> Julia (lines) and MATLAB (markers at peaks/troughs) overlay — without Rayleigh dissipation.</figcaption>
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

The only change from the code above is the far-field term, which now follows the asymmetric
form of equation (3.12): short waves are selected upstream ($x<0$) and long waves downstream
($x>0$). The local term is unchanged (its equivalence to $G(x)/\pi\alpha$ is proved in the
linked page), so we reuse `η_local` from the previous block.

```julia
# ─── Equation (3.12):  Rayleigh-dissipation far-field (asymmetric) ───
# Upstream (x<0): long-wavelength kₗ;  Downstream (x>0): short-wavelength kₛ.
# (At x = 0 both branches give sin(0) = 0, so no separate case is needed.)
η_farfield_rayleigh(x) = -2F₀ / (α * (kₗ - kₛ)) * sin((x < 0 ? kₗ : kₛ) * x)

η_ff_r  = η_farfield_rayleigh.(x)
η_loc_r = η_local.(x)          # identical local term as in eqn. (3.11)

plot(x, η_loc_r .* 1e3; label="local", color="red",
     guidefontsize=16, tickfontsize=14, legendfontsize=14,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlims=(-10,10),
     legend=:outerright, size=(800,400))
plot!(x, η_ff_r .* 1e3; label="far-field", color="blue")
```

```matlab
%% Steady-state with Rayleigh dissipation (eqn. 3.12)
% Far-field: asymmetric sinusoidal terms
eta_far_r = zeros(size(x));
eta_far_r(x < 0) = -2*F0/(alpha*(k_l - k_s)) .* sin(k_l*x(x < 0));
eta_far_r(x > 0) = -2*F0/(alpha*(k_l - k_s)) .* sin(k_s*x(x > 0));

% Local: identical local term as the non-Rayleigh case (eqn. 3.11)
eta_local_r = eta_local;

figure; hold on;
plot(x, eta_local_r*1e3, 'r-', 'LineWidth', 2);
plot(x, eta_far_r*1e3, 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('\eta \times 10^3');
legend('local','far-field'); xlim([-10 10]);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/steady_rayleigh.png" alt="Fig 5b(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 5b(i).</strong> Steady-state far-field and local components with Rayleigh dissipation (Julia). The asymmetric far-field response selects short waves upstream (<em>x</em>&lt;0) and long waves downstream (<em>x</em>&gt;0).</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig5b_overlay.png" alt="Fig 5b(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 5b(ii).</strong> Julia (lines) and MATLAB (markers at peaks/troughs) overlay — Rayleigh dissipation.</figcaption>
</figure>
```
