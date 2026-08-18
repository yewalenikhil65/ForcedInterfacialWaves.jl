# ForcedInterfacialWaves.jl

*Julia implementation of the initial value problem (IVP) for pressure-forced interfacial waves in a two-fluid system. Please see the manuscript file below for details*

Based on:

> Kadari, V.K., Yewale, N., Farsoiya, P.K., Mayya, Y.S. & Dasgupta, R. (2026).
> Interfacial waves from pressure forcing: revisiting classical theories from an IVP perspective.
> *arXiv preprint* [arXiv:2605.12254](https://arxiv.org/abs/2605.12254).
## Citation

```bibtex
@article{kadari2026interfacial,
  title={Interfacial waves from pressure forcing: revisiting classical
         theories from an {IVP} perspective},
  author={Kadari, Vinod Kumar and Yewale, Nikhil and Farsoiya, Palas Kumar
          and Mayya, Y. S. and Dasgupta, Ratul},
  journal={arXiv preprint arXiv:2605.12254},
  year={2026}
}
```

## Problem


```@raw html
<figure style="text-align:center;">
  <img src="assets/Fig3.png" alt="Figure 10 comparison" style="max-width:100%;">
  <!--<figcaption>Figure 3 (manuscript): localised pressure forcing at the interface.</figcaption>-->
</figure>
```

As shown above (Fig. $3$ in the manuscript), a localised pressure $\tilde{p}_e = \tilde{F}_0\delta(\tilde{x})$ force is applied at the interface between two inviscid, incompressible, fluids of infinite depth, both streams moving at uniform speed $U$ rightwards. In the absence of this forcing, the interface is flat and remains at $\tilde{z}=0$. After non-dimensionalisation, the interfacial displacement resulting from the forcing i.e. $\eta(x,t)$ is obtained by evaluating the following time-dependent Fourier integrals (eqns. $3.7$ and $3.8$ in the manuscript)
```math
\begin{aligned}
\dfrac{\sqrt{2\pi}\;\bar{\eta}(k,t)}{F_0} &= \left(\dfrac{1}{-\alpha|k|^2 + \left(1+\rho_r\right)|k| - \left(1-\rho_r\right)}\right) \\
&\quad - \dfrac{1}{\alpha |k|^2 + \left(1 - \rho_r\right)} \left(\dfrac{k}{2}\right) \Bigg( \dfrac{\exp\left[-it\lambda_2(k)\right]}{\lambda_2(k)} + \dfrac{\exp\left[-it\lambda_1(k)\right]}{\lambda_1(k)} \Bigg) \tag{3.7}
\end{aligned}
```
The inverse Fourier integral, leads to the (non-dimensional) interface displacement as a function of time:
```math
\begin{aligned}
	\eta(x,t) = \dfrac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}dk\;\exp\left(ikx\right)\bar{\eta}(k,t), 
  \end{aligned}  \tag{3.8}
```
where $\lambda_{1,2}(k) \equiv k\mp\sqrt{\dfrac{\alpha|k|^3}{1+\rho_r}\;+\;\beta|k|}$ and 

```math
\alpha = \frac{gT}{\rho_l U^4}, \quad
\rho_r = \frac{\rho_u}{\rho_l}, \quad
\beta = \frac{1-\rho_r}{1+\rho_r}.  
```

### Steady-state decomposition (manuscript §3.4)
In this section, the manuscript shows that neglecting the time-dependent terms in eqns. (3.7) and (3.8) and $\textit{without}$ using any Rayleigh dissipation, the steady-state response turns out to be (we exclude all the Dirac delta function terms in eqn. $3.9$ in the manuscript) the following. For proof of this, see [Steady-state proof](steady_proof.md).
```math
\begin{aligned}
\dfrac{\eta_{s}(x)}{F_0} =\dfrac{ \eta^{\text{far-field}}_{s}(x)}{F_0} + \dfrac{\eta^{\text{local}}_{s}(x)}{F_0},\tag{3.11}
\end{aligned}
```
where , 
```math
\begin{aligned}
\dfrac{\eta^{\text{far-field}}_{s}(x)}{F_0} \equiv \dfrac{1}{\alpha(k_l-k_s)}\bigg\{-\sin(k_s|x|)\quad + \quad \sin(k_l|x|)\bigg\}, \\
\dfrac{\eta^{\text{local}}_{s}(x)}{F_0} \equiv  \dfrac{\left(k_l+k_s\right)}{\pi\alpha}\int_{0}^{\infty}dy \dfrac{y\exp\left(-|x|y\right)}{\left(y^2 + k_l^2\right)\left(y^2+k_s^2\right)},\quad x\neq 0,
\end{aligned}
```

```math
k_{l,s}
=
\frac{1+\rho_r}{2\alpha}
\left[
1\pm
\sqrt{
1-\frac{4\alpha\beta}{1+\rho_r}
}
\right].
```

The integral expression for $\frac{\eta^{\text{local}}_{s}(x)}{F_0}$ above is solved numerically (in Julia and MATLAB) using the following codes:

```@raw html
<table style="width:100%"><tr>
<td style="vertical-align:top; width:50%">
<strong>Julia</strong>
<pre><code class="language-julia">using QuadGK, Plots
using ForcedInterfacialWaves

# Local manuscript-style alias for QuadGK.quadgk
const ∫ = quadgk

p = compute_cg_parameters()

# Spatial grid (nondimensional)
x = collect(-10:0.01:10)
filter!(xi -> abs(xi) > 1e-12, x)

# Far-field steady
η_far = @. p.F0 / (p.alpha * (p.k_l - p.k_s)) *
             (-sin(p.k_s * abs(x)) + sin(p.k_l * abs(x)))

# Local steady
η_local = similar(η_far)
for i in eachindex(x)
    integrand = y -> y * exp(-abs(x[i]) * y) /
        ((y^2 + p.k_l^2) * (y^2 + p.k_s^2))
    I, _ = ∫(integrand, 0.0, Inf;
                  atol=p.atol, rtol=p.rtol)
    η_local[i] = p.F0 * (p.k_l + p.k_s) /
                   (π * p.alpha) * I
end

# Plot (×10³)
plot(x, η_far .* 1e3,
     label="η_s^{far-field}", lw=1.5)
plot!(x, η_local .* 1e3,
      label="η_s^{local}", lw=1.5)
xlabel!("x"); ylabel!("η × 10³")
savefig("eta_steady_components.png")
</code></pre>
</td>
<td style="vertical-align:top; width:50%">
<strong>MATLAB</strong>
<pre><code class="language-matlab">% Parameters (CGS)
U = 26.7046; g = 981; T = 72;
rho_l = 1; rho_u = 0.001;
l_c = U^2/g;
alpha = T/(rho_l*U^2*l_c);
rho_r = rho_u/rho_l;
disc = (1+rho_r)^2 - 4*alpha*(1-rho_r);
k_l = ((1+rho_r) + sqrt(disc))/(2*alpha);
k_s = ((1+rho_r) - sqrt(disc))/(2*alpha);
F0 = 0.01*T / (rho_l*U^2*l_c);

% Spatial grid (nondimensional)
x = -10:0.01:10;
x(x == 0) = [];

% Far-field steady
eta_far = F0 / (alpha*(k_l - k_s)) * ...
    (-sin(k_s*abs(x)) + sin(k_l*abs(x)));

% Local steady
eta_local = zeros(size(x));
for i = 1:length(x)
    integrand = @(y) y .* exp(-abs(x(i)).*y) ...
        ./ ((y.^2+k_l^2) .* (y.^2+k_s^2));
    I = integral(integrand, 0, Inf, ...
        'AbsTol', 1e-10, 'RelTol', 1e-8);
    eta_local(i) = F0*(k_l+k_s)/(pi*alpha)*I;
end

% Plot (×10³)
figure; hold on;
plot(x, eta_far*1e3, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', '\eta_s^{far-field}');
plot(x, eta_local*1e3, 'r-', 'LineWidth', 1.5, ...
     'DisplayName', '\eta_s^{local}');
xlabel('x'); ylabel('\eta \times 10^3');
legend('Location','best'); grid on;
</code></pre>
</td>
</tr></table>
```

The resulting figure from both Julia and MATLAB are superimposed  in the following figure, which is a reproduction of fig. $5a$ in the manuscript.

![Figure 5a comparison](assets/Fig5a.png)
*Figure 5a: Comparison of steady wave profiles computed in MATLAB vs Julia.*

In eqn 3.11 , $\eta^{\text{local}}_{s}(x)$ ,  can be shown to be exactly equivalent to Lamb's(1932) classical local function $G(x)$(for the proof, see [here](lamb_gx_equivalence.md)). 

$G(x) \equiv \dfrac{1}{k_{l}-k_{s}}\int_{0}^{\infty}\;dk\;\left(\dfrac{\cos(kx)}{k+k_{s}}-\dfrac{\cos(kx)}{k+k_{l}}\right)$

## Manuscript §4 :
### Zero capillarity limit ($\alpha=  0$) - section  §4.1 in the manuscript

Below eqns. (4.1) and (4.2a,b,c) are obtained by substituting $\alpha=0$ in eqn. (3.8) in the manuscript,

```math
		\begin{align}
				\eta(x,t) &=& \eta_{s}(x) + \eta_{tr}^{(1)}(x,t) + \eta_{tr}^{(2)}(x,t), 
		\end{align}\tag{4.1}
```

```math
      \begin{align}	
        \dfrac{\eta_{s}(x)}{F_0} &\equiv \dfrac{1}{2\pi\left(1+\rho_r\right)}\int_{-\infty}^{\infty}dk\; \dfrac{\exp\left(ikx\right)}{|k| - \beta},\quad 0 < \beta \leq 1  \nonumber\\
        \dfrac{\eta_{tr}^{(1)}(x,t)}{F_0} &\equiv - \dfrac{1}{4\pi(1-\rho_r)}\int_{-\infty}^{\infty}dk\;\dfrac{k\;\exp\left[-i\left(k(t-x) + t\sqrt{\beta|k|}\right)\right]}{k+ \sqrt{\beta|k|}},  \nonumber\\
        \dfrac{\eta_{tr}^{(2)}(x,t)}{F_0} &\equiv - \dfrac{1}{4\pi(1-\rho_r)}\int_{-\infty}^{\infty}dk\;\dfrac{k\;\exp\left[-i\left(k(t-x) - t\sqrt{\beta|k|}\right)\right]}{k- \sqrt{\beta|k|}}. 
      \end{align} \tag{4.2a,b,c}
```
It may be further shown using principal value techniques that (see proof [here](pure_gravity_steady_proof.md))

```math
\begin{aligned}
\frac{\eta_s(x)}{F_0}
={}&
\frac{1}{\pi(1+\rho_r)}
\left[
-\pi\sin\!\left(\beta|x|\right)
+
\int_0^\infty
\frac{y\,e^{-|x|y}}
{\beta^2+y^2}\,dy
\right],
\\[0.4em]
&\qquad
0<\beta\leq1,
\qquad
-\infty<x<\infty ,
\end{aligned}
\tag{4.3}
```

In expression (4.3), $\eta_{s}(x)$ is a symmetric function of $x$, implying a symmetric response upstream and downstream of the forcing at $x=0$. However, a contribution to the steady-state *also comes* from the time-dependent term in eqns. (4.2b),(4.2c). As shown in the  proof([here](pure_gravity_transient_proof.md)) these transient terms may be further simplified to obtain the following analytical representation valid for all $x,t$ i.e.

```math
\begin{aligned}
\eta_{\mathrm{tr}}(x,t)
&\equiv
\eta_{\mathrm{tr}}^{(1)}(x,t)
+
\eta_{\mathrm{tr}}^{(2)}(x,t)
\\
&=
\left\{
\mathbb{T}_1(x)
+\mathbb{T}_2(x,t)
+\mathbb{T}_3(x,t)
+\mathbb{T}_4(x,t)
\right\}F_0,
\end{aligned}
\tag{4.4a}
```

where,

```math
\mathbb{T}_1(x)
\equiv
\mp\frac{1}{1+\rho_r}\sin(\beta x),
\tag{4.4b}
```

```math
\begin{aligned}
\mathbb{T}_2(x,t)
\equiv{}&
-\frac{4\beta^{-1}}{\pi(1+\rho_r)}
\int_0^\infty dv\;
v^2
\frac{
\exp\!\left(\mp2av^2\pm2bv\right)
}{
\beta+\left(2v-\beta^{1/2}\right)^2
}
\\
&\times
\Bigg[
\beta^{1/2}
\cos\!\left(\beta^{1/2}tv\right)
\pm
\left(2v-\beta^{1/2}\right)
\sin\!\left(\beta^{1/2}tv\right)
\Bigg],
\end{aligned}
\tag{4.4c}
```

```math
\begin{aligned}
\mathbb{T}_3(x,t)
\equiv{}&
\frac{\beta^{-1/2}}{\pi(1+\rho_r)}
\left(
1+\frac{t}{2(t-x)}
\right)
\sqrt{\frac{\pi}{2|a|}}
\\
&\times
\Bigg[
\cos\!\left(\frac{b^2}{|a|}\right)
\left\{
\frac{1}{2}
\mp
\mathrm{C}\!\left(
b\sqrt{\frac{2}{\pi|a|}}
\right)
\right\}
\\
&\qquad+
\sin\!\left(\frac{b^2}{|a|}\right)
\left\{
\frac{1}{2}
\mp
\mathrm{S}\!\left(
b\sqrt{\frac{2}{\pi|a|}}
\right)
\right\}
\Bigg],
\end{aligned}
\tag{4.4d}
```

```math
\mathbb{T}_4(x,t)
\equiv
-\frac{1}{\pi(1+\rho_r)}
\int_0^\infty dv\;
\frac{
\cos\!\left(av^2+2bv\right)
}{
v+\beta^{1/2}
},
\tag{4.4e}
```

where, $a \equiv t-x, \; b \equiv \dfrac{t\sqrt{\beta}}{2}$, the upper signs in $\mathbb{T}_1(x),\mathbb{T}_2(x,t)$ are used for $x<t$, while lower signs are for $x > t$. The Fresnel integrals $\mathrm{C}(\cdot)$ and $\mathrm{S}(\cdot)$ in eqn. (4.4d) are defined as

```math
\mathrm{C}\left(b\sqrt{\frac{2}{\pi |a|}}\right) \equiv  \int_{0}^{b\sqrt{\frac{2}{\pi |a|}}}dt\;  \cos\left(\frac{\pi t^2}{2}\right),\quad 
	\mathrm{S}\left(b\sqrt{\frac{2}{\pi |a|}}\right) \equiv  \int_{0}^{b\sqrt{\frac{2}{\pi |a|}}}dt\;  \sin\left(\frac{\pi t^2}{2}\right).\tag{4.4f}
```

In expressions (4.4), the *time-independent term*, $\mathbb{T}_1(x)$ (eqn. 4.4b), is asymmetric with the same amplitude as the first term on the right hand side of eqn. (4.3). As a result, these two terms reinforce each other for $x>0$ but cancel for $x<0$. Further, we note that as $\rho_r\rightarrow 1$ $\left(\beta = \dfrac{1-\rho_r}{1 + \rho_r}\rightarrow 0\right)$, the terms diverge. This is physically reasonable because in this limit ($\rho_r\rightarrow 1$), gravity vanishes and in the absence of capillary forces as well (i.e. $\alpha=0$ that we are currently assuming), there remains no restoring force to resist deformation due to the external pressure. The analytical strategy is clear now: provided one can show that $\mathbb{T}_2(x,t\rightarrow\infty)\rightarrow0$, $\mathbb{T}_3(x,t\rightarrow\infty)\rightarrow0$ and $\mathbb{T}_4(x,t\rightarrow\infty)\rightarrow0$, one obtains the expected steady-state lacking waves upstream (except for small localised deformation of the interface due to the localised integral in eqn. (4.3)) and sinusoidal waves downstream ($x>0$) with wavenumber $\beta$.

The integral expressions $\eta(x,t)$ (eqn. (4.1) of the manuscript), $\eta_s(x)$ (eqn. (4.3) of the manuscript) and $\eta_{tr}(x,t)$ (eqn. (4.4) of the manuscript) are solved numerically (in Julia and MATLAB) using the following codes:

```@raw html
<table style="width:100%"><tr>
<td style="vertical-align:top; width:50%">
<strong>Julia</strong>
<pre><code class="language-julia">using QuadGK, Plots
using ForcedInterfacialWaves

# Local manuscript-style alias for QuadGK.quadgk
const ∫ = quadgk

p = compute_cg_parameters()

# Spatial regions
front_band = 1.0

# Times and spatial resolution
# The entries in time_indices are hundredths of a second.
# time_indices = [1, 7, 15, 40, 60, 100, 300, 500]
time_indices = [1]
Nx_plot = 2001

# Dimensional physical parameters in CGS units
g     = 981.0
rho_l = 1.0
rho_u = 0.001
U     = 26.7046

# Characteristic scales
l_c = U^2 / g
t_c = U / g
F_c = rho_l * U^2 * l_c

# Nondimensional forcing amplitude
T_reference = 72.0
F_nd = T_reference / F_c
F0   = 0.01 * F_nd

# Density parameters
rho_r     = rho_u / rho_l
beta      = (1.0 - rho_r) / (1.0 + rho_r)
sqrt_beta = sqrt(beta)
beta_sq   = beta^2

# Plotting domain
gravity_wavelength = 2.0 * pi / beta
L = 4.0 * gravity_wavelength

# Integral parameters
KMAX_analytical = 100.0
KMAX_cpv        = 100.0

AbsTol_ana = 1.0e-10
RelTol_ana = 1.0e-8

epsilon_cpv = 1.0e-6
AbsTol_cpv  = 1.0e-10
RelTol_cpv  = 1.0e-8

# QuadGK integration
function qintegral(f, a, b; atol=1e-10, rtol=1e-8)
    value, _ = quadgk(f, a, b; atol=atol, rtol=rtol)
    return value
end

for time_index in time_indices

    # Dimensional and nondimensional time
    t_dim = time_index / 100.0
    t     = t_dim / t_c

    # Nondimensional spatial grid
    x = collect(range(-L / 2.0, L / 2.0; length=Nx_plot))

    epsilon_x = 1.0e-6
    x = x[abs.(x) .>= epsilon_x]

    N = length(x)

    # Regions relative to x = t
    a = t .- x

    mask_left  = a .>  front_band
    mask_right = a .< -front_band

    x_abs = abs.(x)

    # Preallocate arrays
    eta_analytical = fill(NaN, N)
    eta_transient  = fill(NaN, N)

    T1 = fill(NaN, N)
    T2 = fill(NaN, N)
    T3 = fill(NaN, N)
    T4 = fill(NaN, N)

    T1_xgt = fill(NaN, N)
    T2_xgt = fill(NaN, N)
    T3_xgt = fill(NaN, N)
    T4_xgt = fill(NaN, N)

    # Analytical steady contribution
    T0a = -pi .* sin.(beta .* x_abs)

    T0b = [
        qintegral(
            y -> exp(-y * xabs_val) * y / (beta_sq + y^2),
            0.0, Inf;
            atol=AbsTol_ana,
            rtol=RelTol_ana
        )
        for xabs_val in x_abs
    ]

    T0 = (T0a .+ T0b) ./ (pi * (1.0 + rho_r))
    eta_steady = F0 .* T0

    # Analytical solution for x < t - front_band
    if any(mask_left)

        idx_left = findall(mask_left)
        a_left   = a[idx_left]

        T1[idx_left] .=
            -sin.(beta .* x[idx_left]) ./ (1.0 + rho_r)

        T2_integrand(v, aa) =
            exp(-2.0 * v^2 * aa + v * t * sqrt_beta) *
            v^2 *
            (
                sqrt_beta * cos(v * t * sqrt_beta) +
                (2.0 * v - sqrt_beta) * sin(v * t * sqrt_beta)
            ) /
            (beta + (2.0 * v - sqrt_beta)^2)

        T2[idx_left] .=
            -4.0 / (pi * (1.0 + rho_r) * beta) .* [
                qintegral(
                    v -> T2_integrand(v, aa),
                    0.0, KMAX_analytical;
                    atol=AbsTol_ana,
                    rtol=RelTol_ana
                )
                for aa in a_left
            ]

        b_left = 0.5 * t * sqrt_beta

        X_left =
            b_left .* sqrt.(2.0 ./ (pi .* a_left))

        prefactor_left =
            (1.0 / (pi * (1.0 + rho_r) * sqrt_beta)) .*
            (1.0 .+ t ./ (2.0 .* a_left))

        T3[idx_left] .=
            prefactor_left .*
            sqrt.(pi ./ (2.0 .* a_left)) .*
            (
                cos.(b_left^2 ./ a_left) .*
                    (0.5 .- fresnelc.(X_left)) .+
                sin.(b_left^2 ./ a_left) .*
                    (0.5 .- fresnels.(X_left))
            )

        T4_integrand(v, aa) =
            cos(v^2 * aa + v * t * sqrt_beta) /
            (v + sqrt_beta)

        T4[idx_left] .=
            -1.0 / (pi * (1.0 + rho_r)) .* [
                qintegral(
                    v -> T4_integrand(v, aa),
                    0.0, KMAX_analytical;
                    atol=AbsTol_ana,
                    rtol=RelTol_ana
                )
                for aa in a_left
            ]

        eta_transient[idx_left] .=
            F0 .* (
                T1[idx_left] .+
                T2[idx_left] .+
                T3[idx_left] .+
                T4[idx_left]
            )

        eta_analytical[idx_left] .=
            eta_steady[idx_left] .+
            eta_transient[idx_left]
    end

    # Analytical solution for x > t + front_band
    if any(mask_right)

        idx_right = findall(mask_right)
        a_right   = a[idx_right]

        T1_xgt[idx_right] .=
            sin.(beta .* x[idx_right]) ./ (1.0 + rho_r)

        T2_xgt_integrand(v, aa) =
            exp(2.0 * v^2 * aa - v * t * sqrt_beta) *
            v^2 *
            (
                sqrt_beta * cos(v * t * sqrt_beta) -
                (2.0 * v - sqrt_beta) * sin(v * t * sqrt_beta)
            ) /
            (beta + (2.0 * v - sqrt_beta)^2)

        T2_xgt[idx_right] .=
            -4.0 / (pi * (1.0 + rho_r) * beta) .* [
                qintegral(
                    v -> T2_xgt_integrand(v, aa),
                    0.0, KMAX_analytical;
                    atol=AbsTol_ana,
                    rtol=RelTol_ana
                )
                for aa in a_right
            ]

        b_right = 0.5 * t * sqrt_beta

        X_right =
            b_right .* sqrt.(2.0 ./ (pi .* abs.(a_right)))

        prefactor_right =
            (1.0 / (pi * (1.0 + rho_r) * sqrt_beta)) .*
            (1.0 .+ t ./ (2.0 .* a_right))

        T3_xgt[idx_right] .=
            prefactor_right .*
            sqrt.(pi ./ (2.0 .* abs.(a_right))) .*
            (
                cos.(b_right^2 ./ abs.(a_right)) .*
                    (0.5 .+ fresnelc.(X_right)) .+
                sin.(b_right^2 ./ abs.(a_right)) .*
                    (0.5 .+ fresnels.(X_right))
            )

        T4_xgt_integrand(v, aa) =
            cos(v^2 * aa + v * t * sqrt_beta) /
            (v + sqrt_beta)

        T4_xgt[idx_right] .=
            -1.0 / (pi * (1.0 + rho_r)) .* [
                qintegral(
                    v -> T4_xgt_integrand(v, aa),
                    0.0, KMAX_analytical;
                    atol=AbsTol_ana,
                    rtol=RelTol_ana
                )
                for aa in a_right
            ]

        eta_transient[idx_right] .=
            F0 .* (
                T1_xgt[idx_right] .+
                T2_xgt[idx_right] .+
                T3_xgt[idx_right] .+
                T4_xgt[idx_right]
            )

        eta_analytical[idx_right] .=
            eta_steady[idx_right] .+
            eta_transient[idx_right]
    end

    # Direct numerical CPV solution
    eta_cpv_numerical = fill(NaN, N)

    for ix in eachindex(x)

        xx = x[ix]

        combined_integrand(k) =
            cos(k * xx) /
            (pi * (1.0 + rho_r) * (k - beta)) -
            k * cos(k * (t - xx) - t * sqrt(beta * k)) /
            (
                2.0 * pi * (1.0 - rho_r) *
                (k - sqrt(beta * k))
            ) -
            k * cos(k * (t - xx) + t * sqrt(beta * k)) /
            (
                2.0 * pi * (1.0 - rho_r) *
                (k + sqrt(beta * k))
            )

        integral_below_pole = qintegral(
            combined_integrand,
            0.0,
            beta - epsilon_cpv;
            atol=AbsTol_cpv,
            rtol=RelTol_cpv
        )

        integral_above_pole = qintegral(
            combined_integrand,
            beta + epsilon_cpv,
            KMAX_cpv;
            atol=AbsTol_cpv,
            rtol=RelTol_cpv
        )

        eta_cpv_numerical[ix] =
            F0 * (integral_below_pole + integral_above_pole)
    end

    # Plot
    idx_left  = findall(mask_left)
    idx_right = findall(mask_right)

    p = plot(
        xlabel=L"x",
        ylabel=L"10^3 y",
        xlims=(-L / 2.0, L / 2.0),
        ylims=(-5.0, 9.0),
        yticks=[-4.0, 0.0, 4.0, 8.0],
        legend=:topleft,
        framestyle=:box,
        size=(500, 335)
    )

    plot!(
        p,
        x[idx_left],
        1.0e3 .* eta_analytical[idx_left];
        linestyle=:dashdot,
        linewidth=2.0,
        label=L"\eta"
    )

    if !isempty(idx_right)
        plot!(
            p,
            x[idx_right],
            1.0e3 .* eta_analytical[idx_right];
            linestyle=:dashdot,
            linewidth=2.0,
            label=""
        )
    end

    plot!(
        p,
        x,
        1.0e3 .* eta_steady;
        linestyle=:solid,
        linewidth=1.5,
        label=L"\eta_s"
    )

    plot!(
        p,
        x[idx_left],
        1.0e3 .* eta_transient[idx_left];
        linestyle=:dot,
        linewidth=1.5,
        label=L"\eta_{\mathrm{tr}}"
    )

    if !isempty(idx_right)
        plot!(
            p,
            x[idx_right],
            1.0e3 .* eta_transient[idx_right];
            linestyle=:dot,
            linewidth=1.5,
            label=""
        )
    end

    marker_spacing = 50

    cpv_left_idx = idx_left[1:marker_spacing:end]

    scatter!(
        p,
        x[cpv_left_idx],
        1.0e3 .* eta_cpv_numerical[cpv_left_idx];
        markersize=3,
        label=L"\eta_{\mathrm{CPV}}"
    )

    if !isempty(idx_right)

        cpv_right_idx = idx_right[1:marker_spacing:end]

        scatter!(
            p,
            x[cpv_right_idx],
            1.0e3 .* eta_cpv_numerical[cpv_right_idx];
            markersize=3,
            label=""
        )
    end

    vline!(
        p,
        [t];
        linestyle=:dash,
        linewidth=1.2,
        label=""
    )

    annotate!(
        p,
        t,
        8.3,
        text(L"x=t", 10, :center)
    )

    display(p)

    # savefig(p, "pure_gravity_ivp_t$(time_index).png")
end
</code></pre>
</td>
<td style="vertical-align:top; width:50%">
<strong>MATLAB</strong>
<pre><code class="language-matlab">% Parameters (CGS)
%% Two-fluid pure-gravity IVP: analytical and numerical CPV
clear;
clc;
close all;

%% Spatial regions
% Analytical expressions are plotted outside a band around x = t.
front_band = 1.0;

%% Times and spatial resolution
% The entries in time_indices are hundredths of a second.
% time_indices = [1, 7, 15, 40, 60, 100, 300, 500]; % Uncomment if you want to reproduce fig.6 of the manuscript.
time_indices = 1; % To plot at a single time
Nx_plot      = 2001;

%% Dimensional physical parameters in CGS units
g     = 981.0;       % Gravitational acceleration [cm/s^2]
rho_l = 1.0;         % Lower-fluid density [g/cm^3]
rho_u = 0.001;       % Upper-fluid density [g/cm^3]
U     = 26.7046;     % Uniform base flow speed [cm/s]

%% Characteristic scales
l_c = U^2/g;
t_c = U/g;
F_c = rho_l*U^2*l_c;

%% Nondimensional forcing amplitude
T_reference = 72.0;
F_nd             = T_reference/F_c;
F0               = 0.01*F_nd;

%% Density parameters
rho_r     = rho_u/rho_l;
beta      = (1.0 - rho_r)/(1.0 + rho_r);
sqrt_beta = sqrt(beta);
beta_sq   = beta^2;

%% Fourier normalisation factors
sqrt_2pi     = sqrt(2.0*pi);
inv_sqrt_2pi = 1.0/sqrt_2pi;

%% Plotting domain
gravity_wavelength = 2.0*pi/beta;
L                  = 4.0*gravity_wavelength;

%% Integral parameters
KMAX_analytical = 100.0;
KMAX_cpv        = 100.0;

AbsTol_ana = 1.0e-10;
RelTol_ana = 1.0e-8;

epsilon_cpv = 1.0e-6;
AbsTol_cpv = 1.0e-10;
RelTol_cpv = 1.0e-8;

%% Loop over requested times
for jt = 1:numel(time_indices)

    time_index = time_indices(jt);

    % Dimensional time [s] and nondimensional time
    t_dim = time_index/100.0;
    t     = t_dim/t_c;

    %% Nondimensional spatial grid
    x = linspace(-L/2.0, L/2.0, Nx_plot);

    % Remove x = 0 because the steady transformed expression
    % becomes numerically sensitive there.
    epsilon_x = 1.0e-6;
    x(abs(x) < epsilon_x) = [];

    N = numel(x);

    %% Regions relative to x = t point
    a = t - x;

    mask_left  = a >  front_band;
    mask_right = a < -front_band;
    mask_outer = mask_left | mask_right;

    x_abs = abs(x);

    %% Preallocate analytical arrays
    eta_analytical = nan(1, N);
    eta_transient  = nan(1, N);

    T1 = nan(1, N);
    T2 = nan(1, N);
    T3 = nan(1, N);
    T4 = nan(1, N);

    T1_xgt = nan(1, N);
    T2_xgt = nan(1, N);
    T3_xgt = nan(1, N);
    T4_xgt = nan(1, N);

    %% Analytical steady contribution T0
    T0a = -pi.*sin(beta.*x_abs);

    T0b = arrayfun(@(x_abs_value) ...
        integral( ...
        @(y) exp(-y.*x_abs_value) ...
        .*y./(beta_sq + y.^2), ...
        0.0, ...
        Inf, ...
        'AbsTol', AbsTol_ana, ...
        'RelTol', RelTol_ana), ...
        x_abs);

    T0 = (T0a + T0b)./(pi.*(1.0 + rho_r));

    eta_steady = F0.*T0;

    %% Analytical solution for x < t - front_band
    if any(mask_left)

        a_left = a(mask_left);

        % T1
        T1(mask_left) = ...
    		-sin(beta.*x(mask_left))./(1.0 + rho_r);

        % T2
        T2_integrand = @(v, aa) ...
    		exp(-2.0.*v.^2.*aa + v.*t.*sqrt_beta) ...
    		.*v.^2 ...
    		.*( ...
    		sqrt_beta.*cos(v.*t.*sqrt_beta) ...
    		+(2.0.*v - sqrt_beta) ...
    		.*sin(v.*t.*sqrt_beta) ...
    		) ...
    		./(beta + (2.0.*v - sqrt_beta).^2);

        T2(mask_left) = ...
    		-4.0./(pi.*(1.0 + rho_r).*beta) ...
    		.*arrayfun(@(aa) ...
    		integral( ...
    		@(v) T2_integrand(v, aa), ...
    		0.0, ...
    		KMAX_analytical, ...
    		'AbsTol', AbsTol_ana, ...
    		'RelTol', RelTol_ana), ...
    		a_left);

        % T3
        b_left = 0.5.*t.*sqrt_beta;

        X_left = b_left ...
    		.*sqrt(2.0./(pi.*a_left));

        prefactor_left = ...
    		(1.0./(pi*(1+rho_r)*sqrt_beta)) ...
    		.*(1.0 + t./(2.0.*a_left));

        T3(mask_left) = ...
    	    prefactor_left ...
    		.*sqrt(pi./(2.0.*a_left)) ...
    		.*( ...
    		cos(b_left.^2./a_left) ...
    		.*(0.5 - fresnelc(X_left)) ...
    		+ ...
    		sin(b_left.^2./a_left) ...
    		.*(0.5 - fresnels(X_left)) ...
    		);

        % T4
        T4_integrand = @(v, aa) ...
    		cos(v.^2.*aa + v.*t.*sqrt_beta) ...
    		./(v + sqrt_beta);

        T4(mask_left) = ...
    		-1.0./(pi.*(1.0 + rho_r)) ...
    		.*arrayfun(@(aa) ...
    		integral( ...
    		@(v) T4_integrand(v, aa), ...
    		0.0, ...
    		KMAX_analytical, ...
    		'AbsTol', AbsTol_ana, ...
    		'RelTol', RelTol_ana), ...
    		a_left);

        eta_transient(mask_left) = F0.*( ...
    		T1(mask_left) ...
    		+ T2(mask_left) ...
    		+ T3(mask_left) ...
    		+ T4(mask_left));

        eta_analytical(mask_left) = ...
    		eta_steady(mask_left) ...
    		+ eta_transient(mask_left);
    end

    %% Analytical solution for x > t + front_band
    if any(mask_right)

        % In the source formulation, a_right = t - x < 0.
        a_right = a(mask_right);

        % T1_xgt
        T1_xgt(mask_right) = ...
    		sin(beta.*x(mask_right))./(1.0 + rho_r);

        % T2_xgt
        T2_xgt_integrand = @(v, aa) ...
    		exp(2.0.*v.^2.*aa - v.*t.*sqrt_beta) ...
    		.*v.^2 ...
    		.*( ...
    		sqrt_beta.*cos(v.*t.*sqrt_beta) ...
    		-(2.0.*v - sqrt_beta) ...
    		.*sin(v.*t.*sqrt_beta) ...
    		) ...
    		./(beta + (2.0.*v - sqrt_beta).^2);

        T2_xgt(mask_right) = ...
    		-4.0./(pi.*(1.0 + rho_r).*beta) ...
    		.*arrayfun(@(aa) ...
    		integral( ...
    		@(v) T2_xgt_integrand(v, aa), ...
    		0.0, ...
    		KMAX_analytical, ...
    		'AbsTol', AbsTol_ana, ...
    		'RelTol', RelTol_ana), ...
    		a_right);

        % T3_xgt
        b_right = 0.5.*t.*sqrt_beta;

        X_right = b_right ...
    		.*sqrt(2.0./(pi.*abs(a_right)));

        prefactor_right = ...
    		(1.0./(pi*(1+rho_r)*sqrt_beta)) ...
    		.*(1.0 + t./(2.0.*a_right));

        T3_xgt(mask_right) = ...
    		prefactor_right ...
    		.*sqrt(pi./(2.0.*abs(a_right))) ...
    		.*( ...
    		cos(b_right.^2./abs(a_right)) ...
    		.*(0.5 + fresnelc(X_right)) ...
    		+ ...
    		sin(b_right.^2./abs(a_right)) ...
    		.*(0.5 + fresnels(X_right)) ...
    		);

        % T4_xgt
        T4_xgt_integrand = @(v, aa) ...
    		cos(v.^2.*aa + v.*t.*sqrt_beta) ...
    		./(v + sqrt_beta);

        T4_xgt(mask_right) = ...
    		-1.0./(pi.*(1.0 + rho_r)) ...
    		.*arrayfun(@(aa) ...
    		integral( ...
    		@(v) T4_xgt_integrand(v, aa), ...
    		0.0, ...
    		KMAX_analytical, ...
    		'AbsTol', AbsTol_ana, ...
    		'RelTol', RelTol_ana), ...
    		a_right);

        eta_transient(mask_right) = F0.*( ...
    		T1_xgt(mask_right) ...
    		+ T2_xgt(mask_right) ...
    		+ T3_xgt(mask_right) ...
    		+ T4_xgt(mask_right));

        eta_analytical(mask_right) = ...
    		eta_steady(mask_right) ...
    		+ eta_transient(mask_right);
    end

    %% Direct numerical CPV solution
    eta_cpv_numerical = nan(1, N);

    for ix = 1:N

        xx = x(ix);

        % Form the complete integrand before integration so that
        % the singular contributions are treated together.
        combined_integrand = @(k) ...
    		cos(k.*xx) ...
    		./(pi.*(1.0 + rho_r).*(k - beta)) ...
    		...
    		-k.*cos( ...
    		k.*(t - xx) - t.*sqrt(beta.*k)) ...
    		./( ...
    		2.0.*pi.*(1.0 - rho_r) ...
    		.*(k - sqrt(beta.*k)) ...
    		) ...
    		...
    		-k.*cos( ...
    		k.*(t - xx) + t.*sqrt(beta.*k)) ...
    		./( ...
    		2.0.*pi.*(1.0 - rho_r) ...
    		.*(k + sqrt(beta.*k)) ...
    		);

        integral_below_pole = integral( ...
    		combined_integrand, ...
    		0.0, ...
    		beta - epsilon_cpv, ...
    		'RelTol', RelTol_cpv, ...
    		'AbsTol', AbsTol_cpv, ...
    		'ArrayValued', true);

        integral_above_pole = integral( ...
    		combined_integrand, ...
    		beta + epsilon_cpv, ...
    		KMAX_cpv, ...
    		'RelTol', RelTol_cpv, ...
    		'AbsTol', AbsTol_cpv, ...
    		'ArrayValued', true);

        eta_cpv_numerical(ix) = F0.*( ...
    		integral_below_pole ...
    		+ integral_above_pole);
    end

    %% Plot analytical and numerical CPV solutions
    figure_size_in = [3.45, 2.30];

    axis_line_width = 0.9;
    curve_line_width = 1.8;

    axis_font_size  = 10;
    label_font_size = 12;
    legend_font_size = 10;

    marker_size    = 3;
    marker_spacing = 50;

    fig = figure( ...
        'Visible', 'on', ...
        'Color', 'w', ...
        'Renderer', 'painters');

    set(fig, ...
        'Units', 'inches', ...
        'Position', [1, 1, figure_size_in]);

    layout = tiledlayout( ...
        fig, ...
        1, ...
        1, ...
        'Padding', 'compact', ...
        'TileSpacing', 'compact');

    ax = nexttile(layout);

    hold(ax, 'on');

    set(ax, ...
        'TickDir', 'out', ...
        'LineWidth', axis_line_width, ...
        'FontSize', axis_font_size, ...
        'XAxisLocation', 'bottom', ...
        'YAxisLocation', 'left');

    %% Plot analytical total solution first
    h_analytical_left = plot( ...
        ax, ...
        x(mask_left), ...
        1.0e3.*eta_analytical(mask_left), ...
        'b-.', ...
        'LineWidth', 2.0);

    if any(mask_right)
        h_analytical_right = plot( ...
    		ax, ...
    		x(mask_right), ...
    		1.0e3.*eta_analytical(mask_right), ...
    		'b-.', ...
    		'LineWidth', 2.0);

        h_analytical_right.HandleVisibility = 'off';
    end
    %% Plot steady solution
    h_steady = plot( ...
        ax, ...
        x, ...
        1.0e3.*eta_steady, ...
        'k-', ...
        'LineWidth', 1.5);

    %% Plot analytical transient contribution
    h_transient_left = plot( ...
        ax, ...
        x(mask_left), ...
        1.0e3.*eta_transient(mask_left), ...
        'm:', ...
        'LineWidth', 1.5);

    if any(mask_right)
        h_transient_right = plot( ...
    		ax, ...
    		x(mask_right), ...
    		1.0e3.*eta_transient(mask_right), ...
    		'm:', ...
    		'LineWidth', 1.5);

        h_transient_right.HandleVisibility = 'off';
    end

    %% Plot direct numerical CPV solution using markers
    index_left = find(mask_left);

    h_numerical_cpv = plot( ...
        ax, ...
        x(index_left), ...
        1.0e3.*eta_cpv_numerical(index_left), ...
        'LineStyle', 'none', ...
        'Color', [0.49, 0.18, 0.56], ...
        'Marker', 'o', ...
        'MarkerIndices', ...
        1:marker_spacing:numel(index_left), ...
        'MarkerSize', marker_size, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'auto', ...
        'LineWidth', curve_line_width - 1.0);

    if any(mask_right)

        index_right = find(mask_right);

        h_numerical_cpv_right = plot( ...
    		ax, ...
    		x(index_right), ...
    		1.0e3.*eta_cpv_numerical(index_right), ...
    		'LineStyle', 'none', ...
    		'Color', [0.49, 0.18, 0.56], ...
    		'Marker', 'o', ...
    		'MarkerIndices', ...
    		1:marker_spacing:numel(index_right), ...
    		'MarkerSize', marker_size, ...
    		'MarkerFaceColor', 'w', ...
    		'MarkerEdgeColor', 'auto', ...
    		'LineWidth', curve_line_width - 1.0);

        h_numerical_cpv_right.HandleVisibility = 'off';
    end

    %% Axis limits
    xlim(ax, [-L/2.0, L/2.0]);
    ylim(ax, [-5.0, 9.0]);
    yticks(ax, [-4.0, 0.0, 4.0, 8.0]);

    %% Mark the front x = t
    y_limits = ylim(ax);

    line( ...
        ax, ...
        [t, t], ...
        y_limits, ...
        'LineStyle', '--', ...
        'Color', [0, 0, 0], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    y_text = y_limits(2) ...
        -0.05.*(y_limits(2) - y_limits(1));

    text( ...
        ax, ...
        t, ...
        y_text, ...
        '$x=t$', ...
        'Interpreter', 'latex', ...
        'FontSize', axis_font_size, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);

    %% Draw the upper and right frame lines
    set(ax, 'Box', 'off');

    delete(findall(ax, 'Tag', 'TopRightFrame'));

    x_limits = xlim(ax);
    y_limits = ylim(ax);

    line( ...
        ax, ...
        [x_limits(1), x_limits(2)], ...
        [y_limits(2), y_limits(2)], ...
        'Color', 'k', ...
        'LineWidth', axis_line_width, ...
        'Clipping', 'off', ...
        'HandleVisibility', 'off', ...
        'Tag', 'TopRightFrame');

    line( ...
        ax, ...
        [x_limits(2), x_limits(2)], ...
        [y_limits(1), y_limits(2)], ...
        'Color', 'k', ...
        'LineWidth', axis_line_width, ...
        'Clipping', 'off', ...
        'HandleVisibility', 'off', ...
        'Tag', 'TopRightFrame');

    %% Legend
    legend_handle = legend( ...
        ax, ...
        [ ...
        h_analytical_left, ...
        h_steady, ...
        h_transient_left, ...
        h_numerical_cpv ...
        ], ...
        { ...
        '$\eta$', ...
        '$\eta_s$', ...
        '$\eta_{\mathrm{tr}}$', ...
        '$\eta_{\mathrm{CPV}}$' ...
        }, ...
        'Location', 'northwest', ...
        'Interpreter', 'latex', ...
        'FontSize', legend_font_size, ...
        'Box', 'on');

    legend_handle.EdgeColor = 0.4.*[1, 1, 1];
    legend_handle.LineWidth = 0.6;
    legend_handle.Color = 'w';
    legend_handle.ItemTokenSize = [12, 9];

    %% Axis labels
    xlabel( ...
        ax, ...
        '$x$', ...
        'Interpreter', 'latex', ...
        'FontSize', label_font_size);

    ylabel( ...
        ax, ...
        '$10^3 y$', ...
        'Interpreter', 'latex', ...
        'FontSize', label_font_size);

    drawnow;
end

</code></pre>
</td>
</tr></table>
```

### Finite-capillarity : ($\alpha > 0$) section §4.2 in the manuscript

We now turn to the case of $\alpha>0$. As $\rho_r<1$, we have both capillary and gravitational forces. For this case, eqn. (3.8) of the manuscript can be rewritten as (Note that the integrals in eqn.(3.8) are folded onto the positive k-axis to get rid of the $|k|$ terms),

```math
\eta(x,t)
=
\eta_s(x)+\eta_{\mathrm{tr}}(x,t),
\tag{4.5a}
```

where,

```math
\frac{\eta_s(x)}{F_0}
\equiv
-\frac{1}{\pi}
\int_0^\infty dk\;
\frac{\cos(kx)}
{\alpha (k-k_l)(k-k_s)},
\tag{4.5b}
```

```math
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0}
\equiv
-\frac{1}{2\pi}
\left[
\mathbb{I}_3(x,t)+\mathbb{I}_4(x,t)
\right],
\tag{4.5c}
```

```math
\begin{aligned}
\mathbb{I}_{3,4}(x,t)
\equiv{}&
-\frac{1+\rho_r}{\alpha}
\int_0^\infty
dk\,
\frac{
\left(k\pm\chi(k)\right)
\cos\left[t\left(k\mp\chi(k)\right)-kx\right]
}{
\left(1+\alpha k^2-\rho_r\right)
\left(k-k_l\right)
\left(k-k_s\right)
},
\end{aligned}
\tag{4.5d}
```

```math
\chi(k)
\equiv
\sqrt{
\beta k+\frac{\alpha}{1+\rho_r}k^3
},
```

```math
k_{l,s}
=
\frac{1+\rho_r}{2\alpha}
\left[
1\pm
\sqrt{
1-\frac{4\alpha\beta}{1+\rho_r}
}
\right].
```

Similar to the previous section ($\alpha=0$ case), the interface shape due to the time-independent response $\eta_s(x)$, given by eqn. (4.5b), is also symmetric about $x=0$, (see the [Steady-state proof](steady_proof.md)). Note that eqn. (4.5b) is identical to eqn. (3.9) of the manuscript after excluding all terms arising from the Dirac delta function. This symmetric response is shown by the black solid curve in panel (a) of Figure 8 of the manuscript.

We now turn to a formal demonstration of the asymmetric cancellations about $x=0$ as $t\rightarrow\infty$. The expression for $\eta_s(x)$ in eqn. (4.5b), after application of principal-value techniques, may be written as (see the [Steady-state proof](steady_proof.md)):

```math
\begin{aligned}
\frac{\eta_s(x)}{F_0}
={}&
\frac{1}{\alpha(k_l-k_s)}
\left[
-\sin(k_s|x|)
+\sin(k_l|x|)
\right]
\\
&+
\left(
\frac{k_l+k_s}{\alpha\pi}
\right)
\int_0^\infty
dy\,
\frac{
y\exp\left(-|x|y\right)
}{
\left(y^2+k_l^2\right)
\left(y^2+k_s^2\right)
}.
\end{aligned}
\tag{4.7}
```

After lengthy calculations involving contour integration and stationary-phase approximation (see the [Capillary-gravity asymmetric cancellation proof](capillary_gravity_asymmetric_cancellation.md)), we may show that 

```math
\begin{aligned}
\frac{\eta_{\mathrm{tr}}(x,t\rightarrow\infty)}{F_0}
={}&
\frac{1}{\alpha(k_l-k_s)}
\left[
-\sin(k_sx)-\sin(k_lx)
\right],
\\
&\qquad
x\in(-\infty,\infty).
\end{aligned}
\tag{4.8}
```

This contribution to the steady state essentially stems from the term $\mathbb{I}_3(x,t)$ in eqn. (4.5d), whereas the term $\mathbb{I}_4(x,t)$ in the same equation tends to zero as $t\rightarrow\infty$. Figure 7 of the manuscript confirms this decay for large time, $t\gg1$.

The sum of eqns. (4.7) and (4.8) yields the final form of the steady-state interface at all $x$.

The asymmetric cancellation in $\eta(x,t\rightarrow\infty)=\eta_s(x)+\eta_{\mathrm{tr}}(x,t\rightarrow\infty)$, upstream ($x<0$) and downstream ($x>0$) of the forcing, may readily be observed by comparing these expressions. We reiterate that the short waves for $x<0$ and the long waves for $x>0$ seen at steady state result from this cancellation.

Unlike the $\alpha=0$ case, it was not possible to obtain closed-form expressions in terms of real integrals for the $\eta_{\mathrm{tr}}(x,t)$ terms in eqn. (4.5d). Hence, we have evaluated these integrals directly numerically in the principal-value sense around the pole(s). 

The integral expressions for $\eta(x,t)$ [eqn. (4.5a) of the manuscript], $\eta_s(x)$ [eqn. (4.5b)], and $\eta_{\mathrm{tr}}(x,t)$ [eqn. (4.5c)] are evaluated numerically using both Julia and MATLAB with the codes provided below. The integrals are computed using a numerical Cauchy principal value (CPV) procedure, in which a small neighborhood of width $\epsilon=10^{-6}$ around each pole, $k=k_s$ and $k=k_l$, is excluded from the numerical integration to avoid direct evaluation at the singularities.

```@raw html
<table style="width:100%"><tr>
<td style="vertical-align:top; width:50%">
<strong>Julia</strong>
<pre><code class="language-julia">using QuadGK, Plots
using ForcedInterfacialWaves

# Local manuscript-style alias for QuadGK.quadgk
const ∫ = quadgk

p = compute_cg_parameters()

# Parameters (CGS)
# Two-fluid capillary-gravity IVP solution

using QuadGK
using Plots
using LaTeXStrings
using Printf

# Number of spatial points
Nx_plot = 2001

# Dimensional physical parameters in CGS units
U     = 26.7046       # Uniform base flow speed [cm/s]
g     = 981.0         # Gravitational acceleration [cm/s^2]
T     = 72.0          # Surface tension [dyn/cm]
rho_l = 1.0           # Lower-fluid density [g/cm^3]
rho_u = 0.001         # Upper-fluid density [g/cm^3]
L     = 75.5996       # Dimensional domain length [cm]

# Characteristic scales
l_c = U^2/g
t_c = U/g
F_c = rho_l*U^2*l_c

# Nondimensional parameters
alpha = T/(rho_l*U^2*l_c)
rho_r = rho_u/rho_l

beta  = (1.0 - rho_r)/(1.0 + rho_r)
gamma_rho = 1.0/(1.0 + rho_r)

# Gravity and capillary wave roots
discriminant = (1.0 + rho_r)^2 -
    4.0*alpha*(1.0 - rho_r)

k_l = ((1.0 + rho_r) + sqrt(discriminant))/(2.0*alpha)
k_s = ((1.0 + rho_r) - sqrt(discriminant))/(2.0*alpha)

# Dimensional wavenumbers and wavelengths
k_l_dim = k_l/l_c
k_s_dim = k_s/l_c

lambda_c = 2.0*pi/k_l_dim    # Capillary wavelength
lambda_g = 2.0*pi/k_s_dim    # Gravity wavelength

@printf("alpha    = %.2e\n", alpha)
@printf("rho_r    = %.2e\n", rho_r)
@printf("k_s      = %.2e\n", k_s)
@printf("k_l      = %.2e\n", k_l)
@printf("lambda_c = %.2e cm\n", lambda_c)
@printf("lambda_g = %.2e cm\n", lambda_g)

# Nondimensional forcing amplitude
F0_dim = 0.01*T
F0    = F0_dim/F_c

# Nondimensional spatial grid
x_grid = collect(range(-L/2.0, L/2.0, length=Nx_plot))/l_c

# Remove x = 0 from the grid
x_grid = x_grid[abs.(x_grid) .>= eps(Float64)]

Nx = length(x_grid)

# Quadrature parameters
epsilon_pv = 1.0e-6
AbsTol     = 1.0e-10
RelTol     = 1.0e-8

k_max        = Inf
k_max_steady = Inf
# k_max_steady = 100.0

# Simulation times stored as hundredths of a second.
# time_indices = [1, 3, 7, 15, 25, 60, 145, 300]
time_indices = [1] # To plot at a single time

for it = 1:length(time_indices)

    time_index = time_indices[it]

    # Dimensional and nondimensional times
    t_dim = time_index/100.0
    t     = t_dim/t_c

    # Allocate theoretical solutions
    eta_ivp    = zeros(size(x_grid))
    eta_s      = zeros(size(x_grid))
    eta_tr     = zeros(size(x_grid))

    # Evaluate the theoretical solution
    for ix = 1:Nx

        x = x_grid[ix]

        # Two-fluid dispersion function
        chi(k) = sqrt(
            beta*k + gamma_rho*alpha*k^3
        )

        # Nominally time-independent contribution
        integrand_eta_s(k) =
            2.0*cos(k*x) /
            (alpha*(k - k_l)*(k - k_s))

        # First time-dependent contribution
        integrand_I3(k) =
            -((1.0 + rho_r)*(k + chi(k))) /
            (1.0 - rho_r + alpha*k^2) *
            cos(k*(t - x) - t*chi(k)) /
            (alpha*(k - k_l)*(k - k_s))

        # Second time-dependent contribution
        integrand_I4(k) =
            -((1.0 + rho_r)*(k - chi(k))) /
            (1.0 - rho_r + alpha*k^2) *
            cos(k*(t - x) + t*chi(k)) /
            (alpha*(k - k_l)*(k - k_s))

        # eta: Complete solution
        # Combine the terms before integration
        total_integrand(k) =
            integrand_eta_s(k) +
            integrand_I3(k) +
            integrand_I4(k)

        # Integration below k_s
        integral_1, _ = quadgk(
            total_integrand,
            0.0,
            k_s - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration between k_s and k_l
        integral_2, _ = quadgk(
            total_integrand,
            k_s + epsilon_pv,
            k_l - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration above k_l
        integral_3, _ = quadgk(
            total_integrand,
            k_l + epsilon_pv,
            k_max,
            atol=AbsTol,
            rtol=RelTol
        )

        eta_ivp[ix] = -F0/(2.0*pi) *
            (integral_1 + integral_2 + integral_3)

        # eta_s: Time-independent solution
        # Integration below k_s
        integral_1_s, _ = quadgk(
            integrand_eta_s,
            0.0,
            k_s - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration between k_s and k_l
        integral_2_s, _ = quadgk(
            integrand_eta_s,
            k_s + epsilon_pv,
            k_l - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration above k_l
        integral_3_s, _ = quadgk(
            integrand_eta_s,
            k_l + epsilon_pv,
            k_max,
            atol=AbsTol,
            rtol=RelTol
        )

        eta_s[ix] = -F0/(2.0*pi) *
            (integral_1_s + integral_2_s + integral_3_s)

        # eta_tr: Time-dependent solution
        # Combine the I3 and I4 terms before integration
        total_integrand_tr(k) =
            integrand_I3(k) +
            integrand_I4(k)

        # Integration below k_s
        integral_1_tr, _ = quadgk(
            total_integrand_tr,
            0.0,
            k_s - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration between k_s and k_l
        integral_2_tr, _ = quadgk(
            total_integrand_tr,
            k_s + epsilon_pv,
            k_l - epsilon_pv,
            atol=AbsTol,
            rtol=RelTol
        )

        # Integration above k_l
        integral_3_tr, _ = quadgk(
            total_integrand_tr,
            k_l + epsilon_pv,
            k_max,
            atol=AbsTol,
            rtol=RelTol
        )

        eta_tr[ix] = -F0/(2.0*pi) *
            (integral_1_tr + integral_2_tr + integral_3_tr)

    end

    # Plotting eta, eta_s and eta_tr
    FIG_SIZE_IN = [3.45, 2.30]   # width x height in inches (good for 2-per-row)
    AX_LW   = 0.9 # axis box
    LW      = 1.8 # curves

    AX_FS   = 10
    FS      = 10    # tick numbers
    LABLE_FS = 12   # x/y labels
    LEG_FS  = 10    # legend

    p = plot(
        x_grid,
        1.0e3*eta_ivp,
        label=L"\eta",
        linestyle=:dashdot,
        linewidth=LW,
        xlabel=L"x",
        ylabel=L"10^3 y",
        xlims=(-10, 10),
        ylims=(-6, 15),
        yticks=[-5, 0, 5, 10, 15],
        legend=:topright,
        framestyle=:box,
        grid=false,
        linewidth_subplot=AX_LW,
        tickfontsize=AX_FS,
        guidefontsize=LABLE_FS,
        legendfontsize=LEG_FS,
        size=(round(Int, FIG_SIZE_IN[1]*144),
              round(Int, FIG_SIZE_IN[2]*144))
    )

    # Theoretical IVP solution (eta_s: Time-independent solution)
    plot!(
        p,
        x_grid,
        1.0e3*eta_s,
        label=L"\eta_s",
        linestyle=:solid,
        linewidth=LW
    )

    # Theoretical IVP solution (eta_tr: Time-dependent solution)
    plot!(
        p,
        x_grid,
        1.0e3*eta_tr,
        label=L"\eta_{tr}",
        linestyle=:dot,
        linewidth=LW
    )

    display(p)

end

</code></pre>
</td>
<td style="vertical-align:top; width:50%">
<strong>MATLAB</strong>
<pre><code class="language-matlab">% Parameters (CGS)
%% Two-fluid capillary-gravity IVP solution
clear;
clc;
close all;

%% Number of spatial points
Nx_plot = 2001;

%% Dimensional physical parameters in CGS units
U     = 26.7046;       % Uniform base flow speed [cm/s]
g     = 981.0;         % Gravitational acceleration [cm/s^2]
T = 72.0;          	   % Surface tension [dyn/cm]
rho_l = 1.0;           % Lower-fluid density [g/cm^3]
rho_u = 0.001;         % Upper-fluid density [g/cm^3]
L     = 75.5996;       % Dimensional domain length [cm]

%% Characteristic scales
l_c = U^2/g;
t_c = U/g;
F_c = rho_l*U^2*l_c;

%% Nondimensional parameters
alpha = T/(rho_l*U^2*l_c);
rho_r = rho_u/rho_l;

beta  = (1.0 - rho_r)/(1.0 + rho_r);
gamma_rho = 1.0/(1.0 + rho_r);

%% Gravity and capillary wave roots
discriminant = (1.0 + rho_r)^2 ...
    - 4.0*alpha*(1.0 - rho_r);

if discriminant <= 0.0
    error(['The steady capillary-gravity roots are not ', ...
        'distinct positive real numbers.']);
end

k_l = ((1.0 + rho_r) + sqrt(discriminant))/(2.0*alpha);
k_s = ((1.0 + rho_r) - sqrt(discriminant))/(2.0*alpha);

%% Dimensional wavenumbers and wavelengths
k_l_dim = k_l/l_c;
k_s_dim = k_s/l_c;

lambda_c = 2.0*pi/k_l_dim;	% Capillary wavelength
lambda_g = 2.0*pi/k_s_dim;	% Gravity wavelength

fprintf('alpha    = %.2e\n', alpha);
fprintf('rho_r    = %.2e\n', rho_r);
fprintf('k_s      = %.2e\n', k_s);
fprintf('k_l      = %.2e\n', k_l);
fprintf('lambda_c = %.2e cm\n', lambda_c);
fprintf('lambda_g = %.2e cm\n', lambda_g);

%% Nondimensional forcing amplitude
F0_dim = 0.01*T;
F0    = F0_dim/F_c;

%% Nondimensional spatial grid
x_grid = linspace(-L/2.0, L/2.0, Nx_plot)/l_c;

% Remove x = 0 from the grid
x_grid(abs(x_grid) < eps) = [];

Nx = numel(x_grid);

%% Quadrature parameters
epsilon_pv = 1.0e-6;
AbsTol     = 1.0e-10;
RelTol     = 1.0e-8;

k_max        = Inf;
k_max_steady = Inf;
% k_max_steady = 100.0;

%% Simulation times stored as hundredths of a second.
% time_indices = [1, 3, 7, 15, 25, 60, 145, 300];
time_indices = 1; % To plot at a single time

for it = 1:numel(time_indices)

    time_index = time_indices(it);

    % Dimensional and nondimensional times
    t_dim = time_index/100.0;
    t     = t_dim/t_c;

    %% Allocate theoretical solutions
    eta_ivp    = zeros(size(x_grid));
    eta_s    = zeros(size(x_grid));
    eta_tr    = zeros(size(x_grid));

    %% Evaluate the theoretical solution
    for ix = 1:Nx

        x = x_grid(ix);

        % Two-fluid dispersion function
        chi = @(k) sqrt( ...
            beta.*k + gamma_rho.*alpha.*k.^3);

        % Nominally time-independent contribution
        integrand_eta_s = @(k) ...
            2.0 .* cos(k.*x) ./ ...
            (alpha.*(k - k_l).*(k - k_s));

        % First time-dependent contribution
        integrand_I3 = @(k) ...
            -((1.0 + rho_r).*(k + chi(k))) ./ ...
            (1.0 - rho_r + alpha.*k.^2) .* ...
            cos(k.*(t - x) - t.*chi(k)) ./ ...
            (alpha.*(k - k_l).*(k - k_s));

        % Second time-dependent contribution
        integrand_I4 = @(k) ...
            -((1.0 + rho_r).*(k - chi(k))) ./ ...
            (1.0 - rho_r + alpha.*k.^2) .* ...
            cos(k.*(t - x) + t.*chi(k)) ./ ...
            (alpha.*(k - k_l).*(k - k_s));

        %% eta: Complete solution
        % Combine the terms before integration
        total_integrand = @(k) ...
            integrand_eta_s(k) ...
            + integrand_I3(k) ...
            + integrand_I4(k);

        % Integration below k_s
        integral_1 = integral( ...
            total_integrand, ...
            0.0, ...
            k_s - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration between k_s and k_l
        integral_2 = integral( ...
            total_integrand, ...
            k_s + epsilon_pv, ...
            k_l - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration above k_l
        integral_3 = integral( ...
            total_integrand, ...
            k_l + epsilon_pv, ...
            k_max, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        eta_ivp(ix) = -F0/(2.0*pi) ...
            * (integral_1 + integral_2 + integral_3);

        %% eta_s: Time-independent solution
        % Integration below k_s
        integral_1_s = integral( ...
            integrand_eta_s, ...
            0.0, ...
            k_s - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration between k_s and k_l
        integral_2_s = integral( ...
            integrand_eta_s, ...
            k_s + epsilon_pv, ...
            k_l - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration above k_l
        integral_3_s = integral( ...
            integrand_eta_s, ...
            k_l + epsilon_pv, ...
            k_max, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        eta_s(ix) = -F0/(2.0*pi) ...
            * (integral_1_s + integral_2_s + integral_3_s);

        %% eta_tr: Time-dependent solution
        % Combine the I3 and I4 terms before integration
        total_integrand_tr = @(k) ...
            integrand_I3(k) ...
            + integrand_I4(k);

        % Integration below k_s
        integral_1_tr = integral( ...
            total_integrand_tr, ...
            0.0, ...
            k_s - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration between k_s and k_l
        integral_2_tr = integral( ...
            total_integrand_tr, ...
            k_s + epsilon_pv, ...
            k_l - epsilon_pv, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        % Integration above k_l
        integral_3_tr = integral( ...
            total_integrand_tr, ...
            k_l + epsilon_pv, ...
            k_max, ...
            'AbsTol', AbsTol, ...
            'RelTol', RelTol);

        eta_tr(ix) = -F0/(2.0*pi) ...
            * (integral_1_tr + integral_2_tr + integral_3_tr);

    end

    %% Plotting eta, eta_s and eta_tr
    FIG_SIZE_IN = [3.45 2.30];   % width x height in inches (good for 2-per-row)
    AX_LW   = 0.9; % axis box
    LW      = 1.8; % curves

    AX_FS   = 10;
    FS      = 10;    % tick numbers
    LABLE_FS = 12;   % x/y labels
    LEG_FS  = 10;     % legend

    fig = figure('Visible','on','Color','w','Renderer','painters');
    set(fig,'Units','inches','Position',[1 1 FIG_SIZE_IN]);
    tl = tiledlayout(fig,1,1,'Padding','compact','TileSpacing','compact');
    ax = nexttile(tl); hold(ax,'on');
    set(ax,'TickDir','out','LineWidth',AX_LW,'FontSize',AX_FS);
    set(ax,'XAxisLocation','bottom','YAxisLocation','left');

    hold on;

    % Theoretical IVP solution (eta: Complete solution)
    h1 = plot( ...
        x_grid, ...
        1.0e3*eta_ivp, ...
        'b-.', ...
        'LineWidth', LW);

    % Theoretical IVP solution (eta_s: Time-independent solution)
    h2 = plot( ...
        x_grid, ...
        1.0e3*eta_s, ...
        'k-', ...
        'LineWidth', LW);

    % Theoretical IVP solution (eta_tr: Time-dependent solution)
    h3 = plot( ...
        x_grid, ...
        1.0e3*eta_tr, ...
        'm:', ...
        'LineWidth', LW);

    lgd = legend( ...
        [h1 h2 h3], ...
        {'$\eta$', '$\eta_s$', '$\eta_{tr}$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', LEG_FS, ...
        'Location', 'northeast');

    % Shorter line samples inside legend
    lgd.ItemTokenSize = [12, 9];

    xlabel( ...
        '$x$', ...
        'Interpreter', 'latex', ...
        'FontSize', LABLE_FS);

    ylabel( ...
        '$10^3 y$', ...
        'Interpreter', 'latex', ...
        'FontSize', LABLE_FS);

    xlim([-10, 10]);
    ylim([-6, 15]);
    yticks([-5, 0, 5, 10, 15]);

    set(ax,'Box','off');  % important

    % draw top/right frame lines
    delete(findall(ax,'Tag','TopRightFrame'));

    xl = xlim(ax);
    yl = ylim(ax);

    line(ax, ...
        [xl(1) xl(2)], ...
        [yl(2) yl(2)], ...
        'Color','k', ...
        'LineWidth',AX_LW, ...
        'Clipping','off', ...
        'HandleVisibility','off', ...
        'Tag','TopRightFrame');

    line(ax, ...
        [xl(2) xl(2)], ...
        [yl(1) yl(2)], ...
        'Color','k', ...
        'LineWidth',AX_LW, ...
        'Clipping','off', ...
        'HandleVisibility','off', ...
        'Tag','TopRightFrame');

end

</code></pre>
</td>
</tr></table>
```



This package evaluates these integrals numerically for two physical regimes:

| Case | Surface tension | Poles | Figure |
|:-----|:--------------:|:------|:------:|
| Pure gravity ($\alpha = 0$) | absent | CPV at $k = \beta$ | 6 |
| Capillary–gravity ($\alpha > 0$) | active | Removable at $k_s$, $k_l$ (combined integrand cancels) | 10 |

## What the code computes

**Pure gravity (α = 0):** The solution is decomposed analytically into $T_0$ (steady), $T_1^\pm$–$T_4^\pm$ (transient), where $T_3$ uses the Fresnel cosine and sine integrals (evaluated in closed form via `FresnelIntegrals.jl`, no quadrature). An independent direct numerical Cauchy principal value (CPV) evaluation verifies the analytical decomposition.

**Capillary–gravity (α > 0):** The combined integrand $\mathbb{I}(k; x, t)$ — which sums the steady part $\eta_s$ and the two transient integrands $\mathbb{I}_3$, $\mathbb{I}_4$ so that their singularities at the gravity root $k_s$ and capillary root $k_l$ cancel — is integrated over $[0,\infty)$ split around both poles. The classical steady solution via residues and a $G(x)$ integral is also computed for large-time comparison.




## Quick start
### Install julia from terminal 
```julia
curl -fsSL https://install.julialang.org | sh
```
```julia
using ForcedInterfacialWaves

# ─── Pure gravity ───
pg = compute_gravity_parameters()       # returns PureGravityParams struct
t_pg = 1.0 / pg.t_c

# Analytical decomposition (left of wavefront):
η_total, η_steady, η_transient = gravity_analytical_left(-2.0, t_pg, pg)

# Independent numerical CPV verification:
η_cpv = gravity_numerical_cpv(-2.0, t_pg, pg)

# Individual T-terms for validation:
T0 = gravity_T0(-2.0, pg)
T3 = gravity_T3_left(-2.0, t_pg, t_pg + 2.0, pg)  # uses fresnelc/fresnels

# ─── Capillary–gravity parameters (all CGS, matching the paper) ───
p = compute_cg_parameters()   # returns CapillaryGravityParams struct

# Evaluate the IVP surface displacement at a single (x, t) point:
η = ivp_surface_elevation(3.0, 110.0, p)

# Compute a full spatial profile using threaded vector-valued quadrature:
x_grid = make_cg_xgrid(p; Nx=2001)
t = 3.0 / p.t_c                        # t_dim = 3 s → nondimensional
η_ivp    = compute_cg_ivp_profile(x_grid, t, p)
η_steady = compute_cg_steady_profile(x_grid, p)
```

## Generating the paper figures

```bash
julia -t auto --project=. scripts/figure6_pure_gravity.jl
julia -t auto --project=. scripts/figure10_capillary_gravity.jl
```

## Contents

```@contents
Pages = ["theory.md", "julia_matlab.md", "validation.md", "api.md"]
Depth = 2
```

