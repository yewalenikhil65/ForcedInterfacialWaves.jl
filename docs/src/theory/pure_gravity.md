# Pure Gravity — Zero Capillarity Limit (Manuscript §4.1)

Setting $\alpha = 0$ removes capillary forces entirely. The dispersion simplifies to $\chi(k) = \sqrt{\beta k}$ and the steady-state integral acquires a single Cauchy principal-value pole at $k = \beta$. This page presents the full analytical decomposition into $T_0$–$T_4^\pm$ terms and an independent numerical CPV verification.

Below eqns. (4.1) and (4.2a,b,c) are obtained by substituting $\alpha=0$ in eqn. (3.8) in the manuscript,

```math
\eta(x,t) = \eta_{s}(x) + \eta_{tr}^{(1)}(x,t) + \eta_{tr}^{(2)}(x,t) \tag{4.1}
```

```math
      \begin{align}
        \dfrac{\eta_{s}(x)}{F_0} &\equiv \dfrac{1}{2\pi\left(1+\rho_r\right)}\int_{-\infty}^{\infty}dk\; \dfrac{\exp\left(ikx\right)}{|k| - \beta},\quad 0 < \beta \leq 1  \nonumber\\
        \dfrac{\eta_{tr}^{(1)}(x,t)}{F_0} &\equiv - \dfrac{1}{4\pi(1-\rho_r)}\int_{-\infty}^{\infty}dk\;\dfrac{k\;\exp\left[-i\left(k(t-x) + t\sqrt{\beta|k|}\right)\right]}{k+ \sqrt{\beta|k|}},  \nonumber\\
        \dfrac{\eta_{tr}^{(2)}(x,t)}{F_0} &\equiv - \dfrac{1}{4\pi(1-\rho_r)}\int_{-\infty}^{\infty}dk\;\dfrac{k\;\exp\left[-i\left(k(t-x) - t\sqrt{\beta|k|}\right)\right]}{k- \sqrt{\beta|k|}}.
      \end{align} \tag{4.2a,b,c}
```
It may be further shown using principal value techniques that (see proof [here](../pure_gravity_steady_proof.md))

```math
\frac{\eta_s(x)}{F_0} = \frac{1}{\pi(1+\rho_r)}\left[-\pi\sin\!\left(\beta|x|\right) + \int_0^\infty \frac{y\,e^{-|x|y}}{\beta^2+y^2}\,dy\right], \quad 0<\beta\leq1,\; -\infty<x<\infty \tag{4.3}
```

In expression (4.3), $\eta_{s}(x)$ is a symmetric function of $x$, implying a symmetric response upstream and downstream of the forcing at $x=0$. However, a contribution to the steady-state *also comes* from the time-dependent term in eqns. (4.2b),(4.2c). As shown in the  proof([here](../pure_gravity_transient_proof.md)) these transient terms may be further simplified to obtain the following analytical representation valid for all $x,t$ i.e.

```math
\eta_{\mathrm{tr}}(x,t) \equiv \eta_{\mathrm{tr}}^{(1)}(x,t) + \eta_{\mathrm{tr}}^{(2)}(x,t) = \left\{\mathbb{T}_1(x) +\mathbb{T}_2(x,t) +\mathbb{T}_3(x,t) +\mathbb{T}_4(x,t)\right\}F_0, \tag{4.4a}
```

where,

```math
\mathbb{T}_1(x) \equiv \mp\frac{1}{1+\rho_r}\sin(\beta x) \tag{4.4b}
```

```math
\mathbb{T}_2(x,t) \equiv -\frac{4\beta^{-1}}{\pi(1+\rho_r)} \int_0^\infty dv\; v^2 \frac{\exp\!\left(\mp2av^2\pm2bv\right)}{\beta+\left(2v-\beta^{1/2}\right)^2} \left[\beta^{1/2}\cos\!\left(\beta^{1/2}tv\right) \pm \left(2v-\beta^{1/2}\right)\sin\!\left(\beta^{1/2}tv\right)\right] \tag{4.4c}
```

```math
\mathbb{T}_3(x,t) \equiv \frac{\beta^{-1/2}}{\pi(1+\rho_r)}\left(1+\frac{t}{2(t-x)}\right)\sqrt{\frac{\pi}{2|a|}} \left[\cos\!\left(\frac{b^2}{|a|}\right)\left\{\frac{1}{2}\mp\mathrm{C}\!\left(b\sqrt{\frac{2}{\pi|a|}}\right)\right\} + \sin\!\left(\frac{b^2}{|a|}\right)\left\{\frac{1}{2}\mp\mathrm{S}\!\left(b\sqrt{\frac{2}{\pi|a|}}\right)\right\}\right] \tag{4.4d}
```

```math
\mathbb{T}_4(x,t) \equiv -\frac{1}{\pi(1+\rho_r)}\int_0^\infty dv\;\frac{\cos\!\left(av^2+2bv\right)}{v+\beta^{1/2}} \tag{4.4e}
```

where, $a \equiv t-x, \; b \equiv \dfrac{t\sqrt{\beta}}{2}$, the upper signs in $\mathbb{T}_1(x),\mathbb{T}_2(x,t)$ are used for $x<t$, while lower signs are for $x > t$. The Fresnel integrals $\mathrm{C}(\cdot)$ and $\mathrm{S}(\cdot)$ in eqn. (4.4d) are defined as

```math
\mathrm{C}\left(b\sqrt{\frac{2}{\pi |a|}}\right) \equiv  \int_{0}^{b\sqrt{\frac{2}{\pi |a|}}}dt\;  \cos\left(\frac{\pi t^2}{2}\right),\quad
	\mathrm{S}\left(b\sqrt{\frac{2}{\pi |a|}}\right) \equiv  \int_{0}^{b\sqrt{\frac{2}{\pi |a|}}}dt\;  \sin\left(\frac{\pi t^2}{2}\right).\tag{4.4f}
```

In expressions (4.4), the *time-independent term*, $\mathbb{T}_1(x)$ (eqn. 4.4b), is asymmetric with the same amplitude as the first term on the right hand side of eqn. (4.3). As a result, these two terms reinforce each other for $x>0$ but cancel for $x<0$. Further, we note that as $\rho_r\rightarrow 1$ $\left(\beta = \dfrac{1-\rho_r}{1 + \rho_r}\rightarrow 0\right)$, the terms diverge. This is physically reasonable because in this limit ($\rho_r\rightarrow 1$), gravity vanishes and in the absence of capillary forces as well (i.e. $\alpha=0$ that we are currently assuming), there remains no restoring force to resist deformation due to the external pressure. The analytical strategy is clear now: provided one can show that $\mathbb{T}_2(x,t\rightarrow\infty)\rightarrow0$, $\mathbb{T}_3(x,t\rightarrow\infty)\rightarrow0$ and $\mathbb{T}_4(x,t\rightarrow\infty)\rightarrow0$, one obtains the expected steady-state lacking waves upstream (except for small localised deformation of the interface due to the localised integral in eqn. (4.3)) and sinusoidal waves downstream ($x>0$) with wavenumber $\beta$.

Putting the $T_0$–$T_4$ terms together gives the full transient-plus-steady solution $\eta(x,t)$ [eqn. (4.1)], and, as an independent cross-check, the original CPV integral can be evaluated directly, bypassing the $T_0$–$T_4$ decomposition entirely. The small difference between $\eta$ and $\eta_{\mathrm{cpv}}$ at $x=-2$ reflects that this point lies just inside the front-exclusion band, where the analytical branch and the direct CPV evaluation are not expected to agree to full precision (see the comparison table below).

## Numerical evaluation

The integral expressions $\eta(x,t)$ (eqn. (4.1) of the manuscript), $\eta_s(x)$ (eqn. (4.3) of the manuscript) and $\eta_{tr}(x,t)$ (eqn. (4.4) of the manuscript) are evaluated numerically, at $x=-2$ and $t_{\dim}=1\,\mathrm{s}$, using both Julia and MATLAB with the codes below. Each analytical term $T_0$–$T_4$ (eqns. 4.3, 4.4b–e) is written out directly; $\eta_s = F_0\,T_0$, $\eta_{tr} = F_0\,(T_1+T_2+T_3+T_4)$, and $\eta = \eta_s + \eta_{tr}$. As an independent check, $\eta$ is also compared against a direct Cauchy-principal-value evaluation $\eta_{\mathrm{CPV}}$ of the combined integrand.

```julia
using QuadGK, FresnelIntegrals

# ─── Nondimensional parameters (α = 0) ───
U, g   = 26.7046, 981.0
ρₗ, ρᵤ = 1.0, 0.001
ρᵣ = ρᵤ / ρₗ
β  = (1 - ρᵣ) / (1 + ρᵣ)
√β = sqrt(β)
F₀ = 0.01 * 72.0 / (ρₗ * U^2 * (U^2 / g))
atol, rtol, K = 1e-10, 1e-8, 100.0

t = 1.0 / (U / g)      # t_dim = 1 s in nondimensional time
x = -2.0

# ─── Analytical terms T₀–T₄ (eqns 4.3, 4.4b–e); sign s: +1 for x<t, −1 for x>t ───
T₀(x) = (-π*sin(β*abs(x)) +
         first(quadgk(y -> exp(-y*abs(x))*y/(β^2 + y^2), 0, Inf; atol=atol, rtol=rtol))) /
        (π*(1 + ρᵣ))
T₁(x, t) = (x < t ? -1 : 1) * sin(β*x) / (1 + ρᵣ)
function T₂(x, t)
    a = t - x; s = a > 0 ? 1.0 : -1.0
    -4/(π*(1 + ρᵣ)*β) * first(quadgk(v ->
        v^2 * exp(-s*2v^2*a + s*v*t*√β) *
        (√β*cos(v*t*√β) + s*(2v - √β)*sin(v*t*√β)) / (β + (2v - √β)^2),
        0, K; atol=atol, rtol=rtol))
end
function T₃(x, t)
    a = t - x; absa = abs(a); s = a > 0 ? 1.0 : -1.0
    b = 0.5*t*√β; X = b*sqrt(2/(π*absa))
    1/(π*(1 + ρᵣ)*√β) * (1 + t/(2a)) * sqrt(π/(2absa)) *
        (cos(b^2/absa)*(0.5 - s*fresnelc(X)) + sin(b^2/absa)*(0.5 - s*fresnels(X)))
end
T₄(x, t) = -1/(π*(1 + ρᵣ)) *
    first(quadgk(v -> cos(v^2*(t - x) + v*t*√β)/(v + √β), 0, K; atol=atol, rtol=rtol, order=15))

println("T₀ = ", T₀(x))
println("T₁ = ", T₁(x, t))
println("T₂ = ", T₂(x, t))
println("T₃ = ", T₃(x, t))
println("T₄ = ", T₄(x, t))

# Assembled solution
η_steady    = F₀ * T₀(x)
η_transient = F₀ * (T₁(x, t) + T₂(x, t) + T₃(x, t) + T₄(x, t))
η           = η_steady + η_transient
println("η         = ", η)
println("η_steady  = ", η_steady)
println("η_transient = ", η_transient)

# ─── Independent CPV check: direct combined integrand G(k;x,t) ───
function η_cpv(x, t)
    ε = 1e-6
    G = k -> begin
        χ = sqrt(β*k); ph = k*(t - x); tχ = t*χ
        cos(k*x)/(π*(1 + ρᵣ)*(k - β)) -
        k*cos(ph - tχ)/(2π*(1 - ρᵣ)*(k - χ)) -
        k*cos(ph + tχ)/(2π*(1 - ρᵣ)*(k + χ))
    end
    lo = first(quadgk(G, 0, β - ε; atol=atol, rtol=rtol))
    hi = first(quadgk(G, β + ε, K; atol=atol, rtol=rtol, order=15))
    F₀ * (lo + hi)
end
ηc = η_cpv(x, t)
println("η_CPV     = ", ηc)
println("|diff|    = ", abs(η - ηc))
```

```text
T₀ = -0.8639502272578199
T₁ = 0.910043043935118
T₂ = 0.005830863275361377
T₃ = 0.0007044858552279276
T₄ = -0.0007003237060223437
η         = 7.212029103433244e-5
η_steady  = -0.0011999023896811143
η_transient = 0.0012720226807154467
η_CPV     = 6.979529999850479e-5
|diff|    = 2.324991035827648e-6
```

```matlab
%% Two-fluid pure-gravity IVP: analytical T0-T4 decomposition
U     = 26.7046;   % Base flow speed [cm/s]
g     = 981.0;     % Gravitational acceleration [cm/s^2]
rho_l = 1.0;       % Lower-fluid density [g/cm^3]
rho_u = 0.001;     % Upper-fluid density [g/cm^3]

% Characteristic scales and nondimensional density ratio
t_c   = U / g;
rho_r = rho_u / rho_l;
beta  = (1.0 - rho_r) / (1.0 + rho_r);
sqrt_beta = sqrt(beta);
F0    = 0.01 * 72.0 / (rho_l * U^2 * (U^2 / g));

% Quadrature tolerances
AbsTol = 1.0e-10;
RelTol = 1.0e-8;
K_MAX  = 100.0;

% Evaluation point and time
t = 1.0 / t_c;
x = -2.0;
a = t - x;          % a > 0 => point lies behind the wavefront

%% T0: steady (time-independent) contribution
T0a = -pi * sin(beta * abs(x));
T0b = integral(@(y) exp(-y * abs(x)) .* y ./ (beta^2 + y.^2), ...
    0, Inf, 'AbsTol', AbsTol, 'RelTol', RelTol);
T0 = (T0a + T0b) / (pi * (1.0 + rho_r));

%% T1: closed-form transient term (left region)
T1 = -sin(beta * x) / (1.0 + rho_r);

%% T2: exponentially-damped oscillatory integral (left region)
T2_integrand = @(v) exp(-2.0 * v.^2 * a + v * t * sqrt_beta) .* v.^2 .* ...
    (sqrt_beta * cos(v * t * sqrt_beta) + ...
     (2.0 * v - sqrt_beta) .* sin(v * t * sqrt_beta)) ./ ...
    (beta + (2.0 * v - sqrt_beta).^2);
T2 = -4.0 / (pi * (1.0 + rho_r) * beta) * ...
    integral(T2_integrand, 0, K_MAX, 'AbsTol', AbsTol, 'RelTol', RelTol);

%% T3: Fresnel-integral contribution (left region)
b = 0.5 * t * sqrt_beta;
X = b * sqrt(2.0 / (pi * a));
prefactor = (1.0 / (pi * (1.0 + rho_r) * sqrt_beta)) * (1.0 + t / (2.0 * a));

T3 = prefactor * sqrt(pi / (2.0 * a)) * ...
    (cos(b^2 / a) * (0.5 - fresnelc(X)) + ...
     sin(b^2 / a) * (0.5 - fresnels(X)));

%% T4: oscillatory integral (left region)
T4_integrand = @(v) cos(v.^2 * a + v * t * sqrt_beta) ./ (v + sqrt_beta);
T4 = -1.0 / (pi * (1.0 + rho_r)) * ...
    integral(T4_integrand, 0, K_MAX, 'AbsTol', AbsTol, 'RelTol', RelTol);

%% Assembled solution
eta_s  = F0 * T0;
eta_tr = F0 * (T1 + T2 + T3 + T4);
eta    = eta_s + eta_tr;

%% Independent CPV check: direct combined integrand G(k;x,t)
eps_pv = 1e-6;
G = @(k) cos(k.*x) ./ (pi*(1+rho_r)*(k - beta)) ...
       - k .* cos(k*(t - x) - t*sqrt(beta*k)) ./ (2*pi*(1-rho_r)*(k - sqrt(beta*k))) ...
       - k .* cos(k*(t - x) + t*sqrt(beta*k)) ./ (2*pi*(1-rho_r)*(k + sqrt(beta*k)));
lo = integral(G, 0, beta - eps_pv, 'AbsTol', AbsTol, 'RelTol', RelTol);
hi = integral(G, beta + eps_pv, K_MAX, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_cpv = F0 * (lo + hi);

fprintf('T0 = %.16g\n', T0);
fprintf('T1 = %.16g\n', T1);
fprintf('T2 = %.16g\n', T2);
fprintf('T3 = %.16g\n', T3);
fprintf('T4 = %.16g\n', T4);
fprintf('eta         = %.16g\n', eta);
fprintf('eta_steady  = %.16g\n', eta_s);
fprintf('eta_transient = %.16g\n', eta_tr);
fprintf('eta_CPV     = %.16g\n', eta_cpv);
fprintf('|diff|      = %.16g\n', abs(eta - eta_cpv));
```

```text
T0 = -0.8639502272578199
T1 = 0.910043043935118
T2 = 0.005830863275361377
T3 = 0.0007044858552279276
T4 = -0.0007003237060223437
eta         = 7.212029103433244e-05
eta_steady  = -0.001199902389681114
eta_transient = 0.001272022680715447
eta_CPV     = 6.979529999850479e-05
|diff|      = 2.324991035827648e-06
```

The full spatial profile at $t = 183.68$ (last panel of manuscript Fig. 6) is computed below by
integrating equations (4.1)–(4.2) directly: the combined-integrand Fourier/CPV inversion gives
the total $\eta(x)$, equation (4.3) gives the steady part $\eta_s(x)$, and $\eta_{tr} = \eta - \eta_s$.
The scripts below plot all three ($\eta$, $\eta_s$, $\eta_{tr}$).

```julia
using QuadGK, FresnelIntegrals
using Plots, LaTeXStrings

# ─── Nondimensional parameters (α = 0) ───
U, g   = 26.7046, 981.0
ρₗ, ρᵤ = 1.0, 0.001
ρᵣ = ρᵤ / ρₗ
β  = (1 - ρᵣ) / (1 + ρᵣ)
F₀ = 0.01 * 72.0 / (ρₗ * U^2 * (U^2 / g))
ε, atol, rtol, K = 1e-6, 1e-10, 1e-8, 100.0

t_fig6 = 183.68
x_grid = collect(range(-12.0, 12.0; length = 2001))
filter!(x -> abs(x) > 1e-6, x_grid)      # exclude x = 0 (integrand singular there)
N = length(x_grid)

# ─── η(x,t): invert the Fourier integral (eqns 4.1–4.2 with α=0) directly ───
# The combined CPV integrand G(k;x,t) below cancels the k=β pole between its three
# terms. Writing each cos(k(t−x) ∓ t√(βk)) via angle-addition exposes cos(kx) and
# sin(kx), so one adaptive quadrature over k serves the whole grid: for each node k
# we compute χ, phases and amplitudes once, then only a cheap sincos(kx) per point.
function η_integrand!(vals, k)
    χ   = sqrt(β * k)                 # dispersion  χ(k) = √(β|k|)
    φ   = t_fig6 * χ
    kt  = k * t_fig6
    a1  = 1.0 / (π * (1 + ρᵣ) * (k - β))          # steady-pole term  → cos(kx)
    a2  = -k / (2π * (1 - ρᵣ) * (k - χ))          # transient branch  (−χ)
    a3  = -k / (2π * (1 - ρᵣ) * (k + χ))          # transient branch  (+χ)
    s_m, c_m = sincos(kt - φ)
    s_p, c_p = sincos(kt + φ)
    cos_coeff = a1 + a2 * c_m + a3 * c_p
    sin_coeff =      a2 * s_m + a3 * s_p
    @inbounds @simd for i in 1:N
        s_kx, c_kx = sincos(k * x_grid[i])
        vals[i] = cos_coeff * c_kx + sin_coeff * s_kx
    end
    return vals
end

acc = zeros(N); buf = zeros(N)
quadgk!(η_integrand!, buf, 0.0, β - ε; atol=atol, rtol=rtol, norm=v->maximum(abs, v)); acc .+= buf
quadgk!(η_integrand!, buf, β + ε, K;   atol=atol, rtol=rtol, order=15, norm=v->maximum(abs, v)); acc .+= buf
η = F₀ .* acc

# ─── η_s(x): steady part (eqn 4.3), one vector quadrature for the local integral ───
loc = zeros(N)
local_integrand!(vals, y) = (@inbounds @simd for i in 1:N
        vals[i] = exp(-y * abs(x_grid[i])) * y / (β^2 + y^2)
    end; vals)
quadgk!(local_integrand!, loc, 0.0, Inf; atol=atol, rtol=rtol, norm=v->maximum(abs, v))
η_steady    = @. F₀ * (-π * sin(β * abs(x_grid)) + loc) / (π * (1 + ρᵣ))
η_transient = η .- η_steady

plot(x_grid, η .* 1e3, label=L"\eta", color="blue", linestyle=:dash,
     guidefontsize=32, tickfontsize=20, xlabel=L"x", ylabel=L"\eta \times 10^{3}",
     legend=:outerright, size=(800,400), legendfontsize=15,
     ylims=(-4.8, 8.2), yticks=[-4, 0, 4, 8])
plot!(x_grid, η_steady .* 1e3, label=L"\eta_s", color="black",
      guidefontsize=32, tickfontsize=20, xlabel=L"x", ylabel=L"\eta \times 10^{3}", xlim=(-12, 12))
plot!(x_grid, η_transient .* 1e3, label=L"\eta_{tr}", color="magenta", linestyle=:dot,
      guidefontsize=32, tickfontsize=20, xlabel=L"x", ylabel=L"\eta \times 10^{3}")
```

```matlab
%% Pure-gravity spatial profile at t = 183.68 (direct Fourier/CPV inversion)
U = 26.7046; g = 981.0;
rho_l = 1.0; rho_u = 0.001;
rho_r = rho_u / rho_l;
beta  = (1 - rho_r) / (1 + rho_r);
F0    = 0.01 * 72.0 / (rho_l * U^2 * (U^2 / g));
eps_pv = 1e-6; AbsTol = 1e-10; RelTol = 1e-8; K_MAX = 100.0;

t = 183.68;
Nx = 2001;
x_grid = linspace(-12, 12, Nx);
x_grid(abs(x_grid) < 1e-6) = [];      % exclude x = 0 (integrand singular there)

% Combined CPV integrand G(k;x,t): the k=beta pole cancels between the three terms.
% 'ArrayValued' integrates the whole x_grid in one adaptive quadrature over k.
G = @(k) cos(k*x_grid) ./ (pi*(1+rho_r)*(k - beta)) ...
       - k .* cos(k*(t - x_grid) - t*sqrt(beta*k)) ./ (2*pi*(1-rho_r)*(k - sqrt(beta*k))) ...
       - k .* cos(k*(t - x_grid) + t*sqrt(beta*k)) ./ (2*pi*(1-rho_r)*(k + sqrt(beta*k)));

I_lo = integral(G, 0, beta - eps_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
I_hi = integral(G, beta + eps_pv, K_MAX, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta  = F0 * (I_lo + I_hi);

% Steady part eta_s (eqn 4.3): closed-form far term + local integral
loc  = integral(@(y) exp(-y.*abs(x_grid)).*y./(beta^2 + y.^2), 0, Inf, ...
                'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_s  = F0 * (-pi*sin(beta*abs(x_grid)) + loc) / (pi*(1 + rho_r));
eta_tr = eta - eta_s;

figure('Position', [100 100 800 400]); hold on;
plot(x_grid, eta*1e3,    'b--', 'LineWidth', 2);
plot(x_grid, eta_s*1e3,  'k-',  'LineWidth', 2);
plot(x_grid, eta_tr*1e3, 'm:',  'LineWidth', 2);
xlabel('x', 'FontSize', 32); ylabel('\eta \times 10^3', 'FontSize', 32);
set(gca, 'FontSize', 20);
xlim([-12 12]); ylim([-4.8 8.2]); yticks([-4 0 4 8]);
legend({'\eta','\eta_s','\eta_{tr}'}, 'Location', 'eastoutside', 'FontSize', 15);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/pure_gravity_fig6.png" alt="Fig 6(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 6(i).</strong> Pure-gravity IVP at <em>t</em> = 183.68 (Julia): total displacement <em>&eta;</em>, steady part <em>&eta;<sub>s</sub></em>, and transient part <em>&eta;<sub>tr</sub></em>.</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig6_comparison.svg" alt="Fig 6(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 6(ii).</strong> Overlay of the Julia and MATLAB profiles, demonstrating that the two computations are equivalent — they agree across the full <em>x</em>-range (shown here at <em>t</em><sub>dim</sub> = 1 s).</figcaption>
</figure>
```

## Asymptotic form of the transient integral(s) and discussion  - section 3.1.1 in the manuscript.

In order to see the transient behaviour of the time-dependent integrals $\mathbb{T}_2(x,t),\mathbb{T}_3(x,t)$ and $\mathbb{T}_4(x,t)$, it is useful to revert to eqns. 3.2(b), (c) in the manuscript which maybe re-written as $\int_0^\infty \frac{k}{k+\sqrt{\beta k}}	\cos\left\{t\left(k+\sqrt{\beta k}\right)-kx \right\}\,dk$ and $\int_0^\infty\frac{k}{k-\sqrt{\beta k}}	\cos\left\{	t\left(k-\sqrt{\beta k}\right)-kx \right\}\,dk$ respectively. It may then be shown using the method of stationary-phase (see [large-time transient-asymptotics proof](../pure_gravity_transient_asymptotic_proof.md)) that for $t\to\infty$:

```math
\begin{align}
        \int_0^\infty \frac{k}{k+\sqrt{\beta k}} \cos\left\{t\left(k+\sqrt{\beta k}\right)-kx\right\}\,dk &\sim -\left(\dfrac{60}{\beta^3}\right)t^{-4} \tag{3.5a}\\[8pt]
        \int_0^\infty \frac{k}{k-\sqrt{\beta k}} \cos\left\{t\left(k-\sqrt{\beta k}\right)-kx\right\}\,dk &\sim 2\pi\beta\sin(\beta x) - \sqrt{\pi\beta}\;t^{-1/2}\cos\left[\frac{\beta(t+x)}{4}-\frac{\pi}{4}\right] \tag{3.5b}
\end{align}
```



Among the three integrals in (4.4c–e), the $t\to\infty$ limit is particularly interesting for $\mathbb{T}_2(x,t)$. For $a\equiv t-x>0$ and in the limit $\beta=1$ (zero density ratio), it decays algebraically as $t^{-1/2}$ at every finite $x$; the other, more general cases may be treated similarly.

At sufficiently large $t\gg1$, the asymptotic form of the cosine integral in $\mathbb{T}_2(x,t\to\infty)$ is obtained from the real part of an integral of the standard form

```math
\mathbb{I}(t)=\int_0^\infty g(\nu)\exp\!\left[t f(\nu)\right]d\nu,
\qquad
f(\nu)=-2\nu^2+\nu+i\nu,
\qquad
g(\nu)=\frac{\nu^2}{1+(2\nu-1)^2},
```

where $\nu$ is continued into the complex plane. The convergence at large time is not immediately apparent because

```math
\exp\!\left[t f(\nu)\right]
=\exp(t/8)\,
 \exp\!\left[-t\left\{\left(\sqrt{2}\nu-\frac{1}{2\sqrt{2}}\right)^2-i\nu\right\}\right].
```

The saddle point is $\nu_0=(1+i)/4$, with $f(\nu_0)=i/4$, $g(\nu_0)=(-1+2i)/20$, and $f''(\nu_0)=-4$. Standard saddle-point integration gives

```math
\mathbb{I}(t\gg1)
\sim
\left(\frac{\pi}{2}\right)^{1/2}
\left(\frac{-1+2i}{20}\right)t^{-1/2}
\exp\!\left(\frac{it}{4}\right).
```

The asymptotic form for the sine term in $\mathbb{T}_2(x,t\to\infty)$ can similarly be obtained by considering the imaginary part of $\mathbb{I}(t)$ with a modified $g(z)$; the algebraic decay of $t^{-1/2}$ is apparent in these results.

For the complete fixed-$x$ asymptotic proof of the two transient integrals, including the pole, stationary-phase, endpoint, and contour-rotation analyses, see the [large-time transient-asymptotics proof](../pure_gravity_transient_asymptotic_proof.md).