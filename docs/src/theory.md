# Theory

Following the notation of the supplementary material (*Interfacial waves from pressure forcing*, Kadari et al.).

## Pure-Gravity IVP ($\alpha = 0$)

Setting $\alpha = 0$ in the general formulation. The dispersion simplifies to $\chi(k) = \sqrt{\beta k}$ and the single pole is at $k = \beta$.

### Analytical decomposition

The solution is decomposed piecewise around the disturbance front $x = t$:

```math
\eta^-(x,t) = F_0[T_0(x) + T_1^- + T_2^- + T_3^- + T_4^-], \quad x < t - \delta_x
```

```math
\eta^+(x,t) = F_0[T_0(x) + T_1^+ + T_2^+ + T_3^+ + T_4^+], \quad x > t + \delta_x
```

where $T_3^\pm$ involve the Fresnel cosine and sine integrals $C(X)$, $S(X)$ (evaluated in closed form, no quadrature), and $\delta_x = 1$.


### Direct numerical CPV

```math
\eta_{\mathrm{CPV}}(x,t) = F_0\left[
  \int_0^{\beta-\varepsilon} G(k;x,t)\,dk +
  \int_{\beta+\varepsilon}^{K_{\max}} G(k;x,t)\,dk
\right]
```

with

```math
G(k;x,t) = \frac{\cos(kx)}{\pi(1+\rho_r)(k-\beta)}
  - \frac{k\cos[k(t-x)-t\sqrt{\beta k}]}{2\pi(1-\rho_r)(k-\sqrt{\beta k})}
  - \frac{k\cos[k(t-x)+t\sqrt{\beta k}]}{2\pi(1-\rho_r)(k+\sqrt{\beta k})}
```

The analytical and CPV evaluations agree to $O(10^{-6})$ in the left region and $O(10^{-5})$ in the right region, providing independent verification of the derivation.

## Capillary–Gravity IVP ($\alpha > 0$)

### Problem setup (§1.1 of supplementary)

A localized pressure $p_e(x,z=0^+) = F_0\,\delta(x)$ acts at the interface of two inviscid, incompressible, irrotational fluids of infinite depth. The nondimensionalisation uses

```math
l_c = \frac{U^2}{g}, \quad t_c = \frac{U}{g}, \quad p_c = \rho_l U^2
```

with nondimensional groups (eq. 2.3, 2.25 of supplementary):

```math
\alpha = \frac{g\widetilde T}{\rho_l U^4}, \quad
\rho_r = \frac{\rho_u}{\rho_l}, \quad
\beta = \frac{1-\rho_r}{1+\rho_r}, \quad
\Delta = \sqrt{1-\frac{4\alpha\beta}{1+\rho_r}}
```

### Dispersion function (eq. 2.13 of supplementary)

```math
\chi(k) = \sqrt{\frac{\alpha|k|^3}{1+\rho_r} + \beta|k|}
```

### Formal solution (eq. 2.17 of supplementary)

```math
\eta(x,t) = \eta_s(x) + \eta_{\mathrm{tr}}(x,t)
```

### Time-independent part (eq. 2.22, 2.26 of supplementary)

```math
\frac{\eta_s(x)}{F_0} = -\frac{1}{2\pi}\left[\mathbb{I}_1(x) + \mathbb{I}_2(x)\right]
```

where

```math
\mathbb{I}_1(x) = \int_0^\infty dk\;\frac{e^{ikx}}{\alpha(k-k_l)(k-k_s)}, \quad
\mathbb{I}_2(x) = \int_0^\infty dk\;\frac{e^{-ikx}}{\alpha(k-k_l)(k-k_s)}
```

Equivalently (the numerical implementation uses this real form):

```math
\frac{\eta_s(x)}{F_0} = -\frac{1}{\pi}\int_0^\infty dk\;\frac{\cos(kx)}{\alpha(k-k_l)(k-k_s)}
```

### Time-dependent part (eq. 2.39–2.41 of supplementary)

```math
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0} = -\frac{1}{2\pi}\left[\mathbb{I}_3(x,t) + \mathbb{I}_4(x,t)\right]
```

```math
\mathbb{I}_3(x,t) = -\frac{(1+\rho_r)}{\alpha}\int_0^\infty dk\;
\frac{(k+\chi(k))\cos[t(k-\chi(k))-kx]}{(\alpha k^2+1-\rho_r)(k-k_l)(k-k_s)}
```

```math
\mathbb{I}_4(x,t) = \int_0^\infty dk\;\frac{k}{\alpha k^2+1-\rho_r}\cdot
\frac{\cos[k(t-x)+t\chi(k)]}{k+\chi(k)}
```

### Poles (eq. 2.23–2.24 of supplementary)

The roots of $\alpha k^2 - (1+\rho_r)k + (1-\rho_r) = 0$ are

```math
k_l = \frac{(1+\rho_r)}{2\alpha}\left[1 + \Delta\right], \quad
k_s = \frac{(1+\rho_r)}{2\alpha}\left[1 - \Delta\right]
```

Both are positive and real for $\alpha \le \alpha_{\max} = \frac{(1+\rho_r)^2}{4(1-\rho_r)}$.

### Numerical evaluation (combined integrand)

The implementation combines all three contributions before integration so that the poles cancel (see `jfm_matlab_codes.tex`, eq. cg\_total\_integrand):

```math
\mathbb{I}(k;x,t) = \frac{2\cos(kx)}{\alpha(k-k_l)(k-k_s)}
  - \frac{(1+\rho_r)(k+\chi)}{1+\alpha k^2-\rho_r}\,
    \frac{\cos[k(t-x)-t\chi]}{\alpha(k-k_l)(k-k_s)}
  - \frac{(1+\rho_r)(k-\chi)}{1+\alpha k^2-\rho_r}\,
    \frac{\cos[k(t-x)+t\chi]}{\alpha(k-k_l)(k-k_s)}
```

```math
\eta(x,t) = -\frac{F_0}{2\pi}\left[
  \int_0^{k_s-\varepsilon}\!\mathbb{I}\,dk +
  \int_{k_s+\varepsilon}^{k_l-\varepsilon}\!\mathbb{I}\,dk +
  \int_{k_l+\varepsilon}^{\infty}\!\mathbb{I}\,dk
\right], \quad \varepsilon = 10^{-6}
```

### Steady solution (evaluated analytically via residues)

```math
G(x) = \frac{1}{k_l-k_s}\int_0^\infty\left[\frac{\cos(kx)}{k+k_s} - \frac{\cos(kx)}{k+k_l}\right]dk
```

```math
\frac{\eta_{\text{steady}}(x)}{F_0} = \begin{cases}
  -\dfrac{2}{\alpha(k_l-k_s)}\sin(k_s x) + \dfrac{G(x)}{\pi\alpha}, & x > 0 \\[6pt]
  -\dfrac{2}{\alpha(k_l-k_s)}\sin(k_l x) + \dfrac{G(x)}{\pi\alpha}, & x < 0
\end{cases}
```

### Steady-state decomposition (manuscript §3.4)
In this section, the manuscript shows that neglecting the time-dependent terms in eqns. (3.7) and (3.8) and $\textit{without}$ using any Rayleigh dissipation, the steady-state response turns out to be (we exclude all the Dirac delta function terms in eqn. $3.9$ in the manuscript) the following. For proof, see [Steady-state proof](steady_proof.md).
```math
\begin{aligned}
\dfrac{\eta_{s}(x)}{F_0} =\dfrac{ \eta^{\text{far-field}}_{s}(x)}{F_0} + \dfrac{\eta^{\text{local}}_{s}(x)}{F_0},\tag{3.11}
\end{aligned}
```
where , 
```math
\begin{aligned}
\dfrac{\eta^{\text{far-field}}_{s}(x)}{F_0} \equiv \dfrac{1}{\alpha(k_l-k_s)}\bigg\{-\sin(k_s|x|)\quad + \quad \sin(k_l|x|)\bigg\} \\
\dfrac{\eta^{\text{local}}_{s}(x)}{F_0} \equiv  \dfrac{\left(k_l+k_s\right)}{\pi\alpha}\int_{0}^{\infty}dy \dfrac{y\exp\left(-|x|y\right)}{\left(y^2 + k_l^2\right)\left(y^2+k_s^2\right)},\quad x\neq 0
\end{aligned}
```
The integral in $\frac{\eta^{\text{local}}_{s}(x)}{F_0}$ is solved numerically (in Julia and MATLAB) as follows:

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
