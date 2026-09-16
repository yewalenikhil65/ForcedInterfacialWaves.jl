# Capillary–Gravity (Manuscript §4.2)

With surface tension present ($\alpha > 0$), both gravity and capillary restoring forces act on the interface. The dispersion relation gains a cubic term, producing two distinct real poles $k_s$ (gravity root) and $k_l$ (capillary root). The key numerical technique is combining all integrand contributions so that the singularities cancel *before* quadrature.

## Formulation

We now turn to the case of $\alpha>0$. As $\rho_r<1$, we have both capillary and gravitational forces. For this case, eqn. (3.8) of the manuscript can be rewritten as (Note that the integrals in eqn.(3.8) are folded onto the positive k-axis to get rid of the $|k|$ terms),

```math
\eta(x,t) = \eta_s(x)+\eta_{\mathrm{tr}}(x,t) \tag{4.5a}
```

where,

```math
\frac{\eta_s(x)}{F_0} \equiv -\frac{1}{\pi}\int_0^\infty dk\;\frac{\cos(kx)}{\alpha (k-k_l)(k-k_s)} \tag{4.5b}
```

```math
\frac{\eta_{\mathrm{tr}}(x,t)}{F_0} \equiv -\frac{1}{2\pi}\left[\mathbb{I}_3(x,t)+\mathbb{I}_4(x,t)\right] \tag{4.5c}
```

```math
\mathbb{I}_{3,4}(x,t) \equiv -\frac{1+\rho_r}{\alpha}\int_0^\infty dk\,\frac{\left(k\pm\chi(k)\right)\cos\left[t\left(k\mp\chi(k)\right)-kx\right]}{\left(1+\alpha k^2-\rho_r\right)\left(k-k_l\right)\left(k-k_s\right)} \tag{4.5d}
```

Summing the steady term (4.5b) and both transient terms (4.5d) into a single combined integrand before quadrature cancels the poles at $k_s,k_l$ analytically, leaving a smooth function to integrate numerically. The combined integrand is then split into three pieces around the (now removable) singularities at $k_s$ and $k_l$ and integrated separately. The CG IVP `WaveSolution` uses this symmetric $\eta_s$ from (4.5b) as `sol.η_steady`; its `sol.η_transient` is the remainder $\eta-\eta_s$.

The asymmetric classical radiation solution is a separate steady reference, obtained with either `IVP(asym_cancel=true)` for the full time-dependent decomposition or `steady(rayleigh_dissipation=true)` for the time-independent profile. It is the long-time result after the transient contribution supplies the asymmetric cancellation described below.

```math
\chi(k) \equiv \sqrt{\beta k+\frac{\alpha}{1+\rho_r}k^3},
```

```math
k_{l,s} = \frac{1+\rho_r}{2\alpha}\left[1\pm\sqrt{1-\frac{4\alpha\beta}{1+\rho_r}}\right].
```

Similar to the previous section ($\alpha=0$ case), the interface shape due to the time-independent response $\eta_s(x)$, given by eqn. (4.5b), is also symmetric about $x=0$, (see the [Steady-state proof](../steady_proof.md)). Note that eqn. (4.5b) is identical to eqn. (3.9) of the manuscript after excluding all terms arising from the Dirac delta function. This symmetric response is shown by the black solid curve in panel (a) of Figure 8 of the manuscript.

We now turn to a formal demonstration of the asymmetric cancellations about $x=0$ as $t\rightarrow\infty$. The expression for $\eta_s(x)$ in eqn. (4.5b), after application of principal-value techniques, may be written as (see the [Steady-state proof](../steady_proof.md)):

```math
\frac{\eta_s(x)}{F_0} = \frac{1}{\alpha(k_l-k_s)}\left[-\sin(k_s|x|)+\sin(k_l|x|)\right] + \frac{k_l+k_s}{\alpha\pi}\int_0^\infty dy\,\frac{y\exp\left(-|x|y\right)}{\left(y^2+k_l^2\right)\left(y^2+k_s^2\right)} \tag{4.7}
```

The second, integral term above is Lamb's $G(x)$ function (up to the prefactor); it is evaluated numerically in the code section below via [`cg_Gx_integral`](@ref).

After lengthy calculations involving contour integration and stationary-phase approximation (see the [Capillary-gravity asymmetric cancellation proof](../capillary_gravity_asymmetric_cancellation.md)), we may show that

```math
\frac{\eta_{\mathrm{tr}}(x,t\rightarrow\infty)}{F_0} = \frac{1}{\alpha(k_l-k_s)}\left[-\sin(k_sx)-\sin(k_lx)\right], \quad x\in(-\infty,\infty). \tag{4.8}
```

This contribution to the steady state essentially stems from the term $\mathbb{I}_3(x,t)$ in eqn. (4.5d), whereas the term $\mathbb{I}_4(x,t)$ in the same equation tends to zero as $t\rightarrow\infty$. Figure 7 of the manuscript confirms this decay for large time, $t\gg1$.

The sum of eqns. (4.7) and (4.8) yields the final form of the steady-state interface at all $x$. The asymmetric cancellation in $\eta(x,t\rightarrow\infty)=\eta_s(x)+\eta_{\mathrm{tr}}(x,t\rightarrow\infty)$, upstream ($x<0$) and downstream ($x>0$) of the forcing, may readily be observed by comparing these expressions. We reiterate that the short waves for $x<0$ and the long waves for $x>0$ seen at steady state result from this cancellation.

Unlike the $\alpha=0$ case, it was not possible to obtain closed-form expressions in terms of real integrals for the $\eta_{\mathrm{tr}}(x,t)$ terms in eqn. (4.5d). Hence, we have evaluated these integrals directly numerically in the principal-value sense around the pole(s).

## Numerical evaluation

The integral expressions for $\eta(x,t)$ [eqn. (4.5a) of the manuscript], $\eta_s(x)$ [eqn. (4.5b)], and $\eta_{\mathrm{tr}}(x,t)$ [eqn. (4.5c)] are evaluated numerically using both Julia and MATLAB with the codes provided below, at $x=3$ and $t=110$. The integrals are computed using a numerical Cauchy principal value (CPV) procedure, in which a small neighborhood of width $\epsilon=10^{-6}$ around each pole, $k=k_s$ and $k=k_l$, is excluded from the numerical integration to avoid direct evaluation at the singularities.

```julia
using QuadGK

# ─── Nondimensional parameters (α > 0) ───
U, g, T = 26.7046, 981.0, 72.0
ρₗ, ρᵤ  = 1.0, 0.001
l_c = U^2 / g
α   = T / (ρₗ * U^2 * l_c)
ρᵣ  = ρᵤ / ρₗ
β   = (1 - ρᵣ) / (1 + ρᵣ)
γᵨ  = 1 / (1 + ρᵣ)
F₀  = 0.01 * T / (ρₗ * U^2 * l_c)

# Gravity (kₛ) and capillary (kₗ) roots
Δ  = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
kₗ = ((1 + ρᵣ) + √Δ) / (2α)
kₛ = ((1 + ρᵣ) - √Δ) / (2α)

ε, atol, rtol = 1e-6, 1e-10, 1e-8
x, t = 3.0, 110.0

# Two-fluid dispersion χ(k)
χ(k) = sqrt(β*k + γᵨ*α*k^3)

# Combined integrand 𝕀(k;x,t): steady (4.5b) + 𝕀₃ + 𝕀₄ (4.5d).
# Summing before quadrature cancels the poles at kₛ, kₗ.
function 𝕀(k, x, t)
    invₚ = 1 / (α * (k - kₗ) * (k - kₛ))          # pole factor
    inv_d = (1 + ρᵣ) / (1 + α*k^2 - ρᵣ)           # dispersion factor
    c = χ(k); ph = k*(t - x); tc = t*c
    2cos(k*x)*invₚ - inv_d*(k + c)*cos(ph - tc)*invₚ - inv_d*(k - c)*cos(ph + tc)*invₚ
end

println("𝕀(k=2; x=3, t=110) = ", 𝕀(2.0, x, t))

# CPV split around the two removable poles kₛ, kₗ
I₁ = first(quadgk(k -> 𝕀(k, x, t), 0, kₛ - ε; atol=atol, rtol=rtol))
I₂ = first(quadgk(k -> 𝕀(k, x, t), kₛ + ε, kₗ - ε; atol=atol, rtol=rtol))
I₃ = first(quadgk(k -> 𝕀(k, x, t), kₗ + ε, Inf; atol=atol, rtol=rtol, order=15))
println("I₁ = ", I₁)
println("I₂ = ", I₂)
println("I₃ = ", I₃)

η = -F₀ / (2π) * (I₁ + I₂ + I₃)          # full IVP, eqn (4.5a)

# Steady η_s (4.5b) via Lamb's G(x): symmetric far-field + local integral (eqn 4.7)
Gₓ = first(quadgk(k -> cos(k*x)/(k + kₛ) - cos(k*x)/(k + kₗ), 0, Inf; atol=atol, rtol=rtol)) / (kₗ - kₛ)
η_steady    = F₀/(α*(kₗ - kₛ)) * (-sin(kₛ*abs(x)) + sin(kₗ*abs(x))) + F₀*Gₓ/(π*α)
η_transient = η - η_steady

println("η_ivp          = ", η)
println("η_s (4.5b)     = ", η_steady)
println("η_transient    = ", η_transient)

# Asymmetric classical radiation steady (long-time reference), x > 0
η_classical = F₀ * (-2/(α*(kₗ - kₛ)) * sin(kₛ*x) + Gₓ/(π*α))
println("η_classical    = ", η_classical)
println("G(3) = ", Gₓ)
```

```text
𝕀(k=2; x=3, t=110) = -2.4925901110201236
I₁ = -10.469830630677315
I₂ = -3.3709095739317627
I₃ = 5.152909149338946
η_ivp          = 0.001920386718354745
η_s (4.5b)     = -0.0005772670229863901
η_transient    = 0.002497653741341135
η_classical    = 0.0018388025677065956
G(3) = 0.011715104659267
```

```matlab
%% Two-fluid capillary-gravity IVP via CPV integration
U     = 26.7046;   % Base flow speed [cm/s]
g     = 981.0;     % Gravitational acceleration [cm/s^2]
T     = 72.0;      % Surface tension [dyn/cm]
rho_l = 1.0;       % Lower-fluid density [g/cm^3]
rho_u = 0.001;     % Upper-fluid density [g/cm^3]

% Characteristic scales and nondimensional parameters
l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
beta      = (1.0 - rho_r) / (1.0 + rho_r);
gamma_rho = 1.0 / (1.0 + rho_r);

% Gravity and capillary wave roots
discriminant = (1.0 + rho_r)^2 - 4.0 * alpha * (1.0 - rho_r);
k_l = ((1.0 + rho_r) + sqrt(discriminant)) / (2.0 * alpha);
k_s = ((1.0 + rho_r) - sqrt(discriminant)) / (2.0 * alpha);

F0 = 0.01 * T / (rho_l * U^2 * l_c);

% Quadrature parameters
epsilon_pv = 1.0e-6;
AbsTol     = 1.0e-10;
RelTol     = 1.0e-8;
k_max      = Inf;

% Evaluation point and time
x = 3.0;
t = 110.0;

% Two-fluid dispersion function
chi = @(k) sqrt(beta * k + gamma_rho * alpha * k.^3);

% Combined integrand: steady term + two transient terms.
% Summing before integration cancels the poles at k_s, k_l.
total_integrand = @(k) ...
    2.0 * cos(k * x) ./ (alpha * (k - k_l) .* (k - k_s)) ...
    - (1.0 + rho_r) * (k + chi(k)) ./ (1.0 - rho_r + alpha * k.^2) .* ...
      cos(k * (t - x) - t * chi(k)) ./ (alpha * (k - k_l) .* (k - k_s)) ...
    - (1.0 + rho_r) * (k - chi(k)) ./ (1.0 - rho_r + alpha * k.^2) .* ...
      cos(k * (t - x) + t * chi(k)) ./ (alpha * (k - k_l) .* (k - k_s));

% Split the CPV integral around the two removable poles
I1 = integral(total_integrand, 0, k_s - epsilon_pv, ...
    'AbsTol', AbsTol, 'RelTol', RelTol);
I2 = integral(total_integrand, k_s + epsilon_pv, k_l - epsilon_pv, ...
    'AbsTol', AbsTol, 'RelTol', RelTol);
I3 = integral(total_integrand, k_l + epsilon_pv, k_max, ...
    'AbsTol', AbsTol, 'RelTol', RelTol);

eta_ivp = -F0 / (2.0 * pi) * (I1 + I2 + I3);

% Steady solution via Lamb's G(x)
G_integrand = @(k) cos(k * x) ./ (k + k_s) - cos(k * x) ./ (k + k_l);
G_x = integral(G_integrand, 0, Inf, 'AbsTol', AbsTol, 'RelTol', RelTol) / ...
    (k_l - k_s);

eta_s = F0 / (alpha * (k_l - k_s)) * ...
    (-sin(k_s * abs(x)) + sin(k_l * abs(x))) + ...
    F0 * G_x / (pi * alpha);
eta_transient = eta_ivp - eta_s;

% Asymmetric classical radiation solution for long-time comparison
eta_classical = F0 * (-2.0 / (alpha * (k_l - k_s)) * sin(k_s * x) + ...
    G_x / (pi * alpha));

fprintf('I(k=2; x=3, t=110) = %.16g\n', total_integrand(2));
fprintf('I1 = %.16g\n', I1);
fprintf('I2 = %.16g\n', I2);
fprintf('I3 = %.16g\n', I3);
fprintf('eta_ivp          = %.16g\n', eta_ivp);
fprintf('eta_s (4.5b)     = %.16g\n', eta_s);
fprintf('eta_transient    = %.16g\n', eta_transient);
fprintf('eta_classical    = %.16g\n', eta_classical);
fprintf('G(3) = %.16g\n', G_x);
```

```text
I(k=2; x=3, t=110) = -2.492590111020124
I1 = -10.46983063067731
I2 = -3.370909573932243
I3 = 5.152903874765406
eta_ivp          = 0.001920387884263914
eta_s (4.5b)     = -0.0005772670409069884
eta_transient    = 0.002497654925170902
eta_classical    = 0.001838802549785997
G(3) = 0.011715099029345
```

The full spatial profile — IVP solution, its steady part, and the transient remainder — reproduces Figure 10 of the manuscript.

## $\mathbb{I}_4$ transient (Fig. 7)

Figure 7 of the manuscript shows the transient component $-\mathbb{I}_4(x,t)$ at an early time $t = 0.37$, confirming the decay of $\mathbb{I}_4$ as $t\to\infty$.

```julia
using QuadGK, Plots, LaTeXStrings

# Standalone regularized 𝕀₄ profile (eqn 4.5d). The pole at k=β-type roots is
# cancelled analytically, so 𝕀₄ is a single integral over [0,∞). The integrand
# factors as amp(k)·[cos(t(k+χ))cos(kx) + sin(t(k+χ))sin(kx)], so one vector-valued
# quadrature over k serves a whole spatial chunk (χ and the time phase are computed
# once per node); the grid is split across threads for speed.
function cg_I4_profile(x_grid, t)
    U, g, T = 26.7046, 981.0, 72.0
    ρₗ, ρᵤ  = 1.0, 0.001
    l_c = U^2 / g
    α   = T / (ρₗ * U^2 * l_c)
    ρᵣ  = ρᵤ / ρₗ
    β   = (1 - ρᵣ) / (1 + ρᵣ)
    γᵨ  = 1 / (1 + ρᵣ)
    atol, rtol = 1e-10, 1e-8
    χ(k) = sqrt(β*k + γᵨ*α*k^3)

    function chunk!(x_chunk)
        M = length(x_chunk); out = zeros(M)
        function integrand!(vals, k)
            k == 0 && (fill!(vals, 0.0); return vals)
            c = χ(k)
            amp = k / ((k + c) * (1 + α*k^2 - ρᵣ))
            s_ph, c_ph = sincos(t*(k + c))
            @inbounds @simd for i in 1:M
                s_kx, c_kx = sincos(k * x_chunk[i])
                vals[i] = amp * (c_ph*c_kx + s_ph*s_kx)
            end
            vals
        end
        quadgk!(integrand!, out, 0.0, Inf; atol=atol, rtol=rtol, order=15,
                norm=v->maximum(abs, v))
        out
    end

    N = length(x_grid); I4 = Vector{Float64}(undef, N)
    nchunks = min(Threads.nthreads(), N)
    clen = cld(N, nchunks)
    Threads.@threads :static for c in 1:nchunks
        lo = (c-1)*clen + 1; hi = min(c*clen, N)
        lo <= hi && copyto!(view(I4, lo:hi), chunk!(view(x_grid, lo:hi)))
    end
    return I4
end

x_grid = collect(range(-15.0, 15.0; length=2001))
filter!(x -> abs(x) > 1e-12, x_grid)
t_I4 = 0.37
I4 = cg_I4_profile(x_grid, t_I4)

plot(x_grid, -I4; color="purple", linewidth=3,
     guidefontsize=16, tickfontsize=14,
     xlabel=L"x", ylabel=L"-\mathbb{I}_{4}",
     xlims=(-10,10), ylims=(-2,14), size=(800,400))
```

```matlab
%% Standalone regularized I4 profile at t = 0.37 (Fig 7)
U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0; rho_u = 0.001;
l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
beta      = (1 - rho_r) / (1 + rho_r);
gamma_rho = 1 / (1 + rho_r);
AbsTol = 1e-10; RelTol = 1e-8;
chi = @(k) sqrt(beta*k + gamma_rho*alpha*k.^3);

t = 0.37;
x_grid = linspace(-15, 15, 2001);
x_grid(abs(x_grid) < 1e-12) = [];

% Regularized I4 integrand (pole cancelled): single integral over [0, Inf).
% 'ArrayValued' integrates the whole x_grid in one adaptive quadrature over k.
I4_integrand = @(k) (k == 0) * zeros(1, numel(x_grid)) + (k ~= 0) * ...
    ( k .* cos(t*(k + chi(k)) - k*x_grid) ./ ((k + chi(k)) * (1 + alpha*k^2 - rho_r)) );
I4 = integral(I4_integrand, 0, Inf, 'ArrayValued', true, ...
              'AbsTol', AbsTol, 'RelTol', RelTol);

figure;
plot(x_grid, -I4, 'Color', [0.5 0 0.5], 'LineWidth', 3);
xlabel('x'); ylabel('-I_4');
xlim([-10 10]); ylim([-2 14]);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/cg_I4_fig7.png" alt="Fig 7(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 7(i).</strong> Transient component <em>-&#x1D540;<sub>4</sub></em> at <em>t</em> = 0.37 (Julia).</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig7_overlay.png" alt="Fig 7(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 7(ii).</strong> Julia (line) and MATLAB (markers) overlay, demonstrating the two computations are equivalent.</figcaption>
</figure>
```
This is plotted ($t=0.37 profile$) as Fig 7 in the manuscript showing time evolution of $I_4(x, t)$ from eqn. 4.5(d) for $\rho_r = 0.001$ and $\alpha = 0.1389$. 

## Full IVP profile (Fig. 8)

The capillary–gravity IVP has the transient contribution $\eta_{\mathrm{tr}}$. The follwoing code blocks in Julia/MATLAB reproduce the fig. 8 in the manuscript that plots $\eta$ and $\eta_{tr}$ at $t=367.35$

```julia
using QuadGK, Plots, LaTeXStrings

# Full CG IVP profile: replicates compute_cg_ivp_profile + compute_cg_steady_profile
# using the package's :threaded_vector structure. Each thread processes a spatial
# chunk: IVP η from the combined CPV integrand (eqn 4.5a–d, 3 pole-split quadratures)
# and steady G(x) for η_s (eqn 4.5b via Lamb's formula) are computed together per chunk.
function cg_ivp_profile(x_grid, t)
    U, g, T = 26.7046, 981.0, 72.0
    ρₗ, ρᵤ  = 1.0, 0.001
    l_c = U^2 / g
    α   = T / (ρₗ * U^2 * l_c)
    ρᵣ  = ρᵤ / ρₗ
    β   = (1 - ρᵣ) / (1 + ρᵣ)
    γᵨ  = 1 / (1 + ρᵣ)
    F₀  = 0.01 * T / (ρₗ * U^2 * l_c)
    Δ   = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
    kₗ  = ((1 + ρᵣ) + √Δ) / (2α)
    kₛ  = ((1 + ρᵣ) - √Δ) / (2α)
    ε, atol, rtol = 1e-6, 1e-10, 1e-8
    χ(k) = sqrt(β*k + γᵨ*α*k^3)

    N = length(x_grid)
    η = Vector{Float64}(undef, N)
    G = Vector{Float64}(undef, N)  # Lamb G(x) for η_s

    nchunks = min(Threads.nthreads(), N); clen = cld(N, nchunks)
    Threads.@threads :static for ci in 1:nchunks
        lo = (ci-1)*clen + 1; hi = min(ci*clen, N); lo > hi && continue
        xc = view(x_grid, lo:hi); M = length(xc)

        # IVP η: combined integrand 𝕀(k;x,t), split at kₛ±ε and kₗ±ε
        I1 = zeros(M); I2 = zeros(M); I3 = zeros(M)
        function ivp!(vals, k)
            c = χ(k); invp = 1/(α*(k-kₗ)*(k-kₛ)); invd = (1+ρᵣ)/(1+α*k^2-ρᵣ)
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
η, η_s, η_tr = cg_ivp_profile(x_grid, t_fig8)

plot(x_grid, η .* 1e3; label=L"\eta", color="blue", ls=:dash, linewidth=3,
     guidefontsize=16, tickfontsize=14, legendfontsize=14,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}",
     xlims=(-10,10), ylims=(-4.8, 8.2), yticks=[-4, 0, 4, 8],
     legend=:outerright, size=(800,400))
plot!(x_grid, η_tr .* 1e3; label=L"\eta_{tr}", color="magenta", ls=:dot, linewidth=3)
```

```matlab
%% Full CG IVP profile at t = 367.35 (Fig 8)
U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0; rho_u = 0.001;
l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
gamma_rho = 1 / (1 + rho_r);
discriminant = (1 + rho_r)^2 - 4*alpha*(1 - rho_r);
k_l = ((1 + rho_r) + sqrt(discriminant)) / (2*alpha);
k_s = ((1 + rho_r) - sqrt(discriminant)) / (2*alpha);
F0  = 0.01 * T / (rho_l * U^2 * l_c);
epsilon_pv = 1e-6; AbsTol = 1e-10; RelTol = 1e-8;
chi = @(k) sqrt((1-rho_r)/(1+rho_r)*k + gamma_rho*alpha*k.^3);

x_grid = linspace(-15, 15, 2001); x_grid(abs(x_grid) < 1e-12) = [];
t = 367.35;

% Combined integrand 𝕀(k;x,t) — ArrayValued integrates the whole x_grid at once.
combined = @(k) ...
    2*cos(k*x_grid) ./ (alpha*(k - k_l).*(k - k_s)) ...
    - (1+rho_r)*(k + chi(k)).*cos(k*(t - x_grid) - t*chi(k)) ./ ...
      ((1 - rho_r + alpha*k^2)*alpha.*(k - k_l).*(k - k_s)) ...
    - (1+rho_r)*(k - chi(k)).*cos(k*(t - x_grid) + t*chi(k)) ./ ...
      ((1 - rho_r + alpha*k^2)*alpha.*(k - k_l).*(k - k_s));

I_lo = integral(combined, 0, k_s - epsilon_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
I_mi = integral(combined, k_s + epsilon_pv, k_l - epsilon_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
I_hi = integral(combined, k_l + epsilon_pv, Inf, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_8 = -F0/(2*pi) * (I_lo + I_mi + I_hi);

% Steady η_s (eqn 4.5b) via Lamb G(x) — ArrayValued
G_int = @(k) (cos(k*x_grid)./(k + k_s) - cos(k*x_grid)./(k + k_l)) / (k_l - k_s);
G_x   = integral(G_int, 0, Inf, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_s_8  = F0/(alpha*(k_l - k_s)) .* (-sin(k_s*abs(x_grid)) + sin(k_l*abs(x_grid))) + F0*G_x/(pi*alpha);
eta_tr_8 = eta_8 - eta_s_8;

figure; hold on;
plot(x_grid, eta_8*1e3,    'b--', 'LineWidth', 3);
plot(x_grid, eta_tr_8*1e3, 'm:',  'LineWidth', 3);
xlabel('x'); ylabel('\eta \times 10^3');
legend('\eta', '\eta_{tr}'); xlim([-10 10]);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/cg_ivp_fig8.png" alt="Fig 8(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 8(i).</strong> Capillary–gravity IVP at <em>t</em> = 367.35 (Julia): total displacement <em>&eta;</em> and transient part <em>&eta;<sub>tr</sub></em>.</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig8_overlay.png" alt="Fig 8(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 8(ii).</strong> Julia (lines) and MATLAB (markers) overlay, demonstrating the two computations are equivalent.</figcaption>
</figure>
```

## Comparison with nonlinear simulations (Fig. 10)

The IVP solution is compared against a nonlinear simulation (Basilisk, Navier–Stokes/VOF) at $t_{\dim} = 25$ s. Simulation data are stored in `notebooks/if_25.csv`; valid time indices are $t_{\dim} \in \{1, 3, 7, 15, 25, 60, 145, 300\}$ s.

```julia
using DelimitedFiles, Plots, LaTeXStrings

# Nondimensional time and length scales (no API)
U, g = 26.7046, 981.0
t_c = U / g          # characteristic time [s]
l_c = U^2 / g        # characteristic length [cm]

t_dim = 25
t_sim = t_dim / (100 * t_c)   # nondimensional (data in CGS: 1 cm = l_c)

# IVP η using the same cg_ivp_profile driver as Fig. 8 (defined above)
x_grid = collect(range(-15.0, 15.0; length=2001)); filter!(x -> abs(x) > 1e-12, x_grid)
η_sim, _, _ = cg_ivp_profile(x_grid, t_sim)

# Basilisk simulation data: column 6 = x [cm], column 7 = y [cm]
data  = sortslices(readdlm(joinpath("notebooks", "if_25.csv"), ',', Float64; skipstart=1),
                   dims=1, by=r -> r[6])
x_bsk = data[:, 6] ./ l_c
y_bsk = data[:, 7] ./ l_c

plot(x_grid, η_sim .* 1e3; label=L"\eta", color="blue", ls=:dash, linewidth=3,
     guidefontsize=16, tickfontsize=14, legendfontsize=14,
     xlabel=L"x", ylabel=L"\eta \times 10^{3}",
     xlims=(-6,10), ylims=(-6, 8.2), yticks=[-4, 0, 4, 8],
     legend=:outerright, size=(800,400))
plot!(x_bsk, y_bsk .* 1e3; label="Simulation", color="red", ls=:dot, linewidth=3)
```

```matlab
%% IVP vs nonlinear simulation at t_dim = 25 s (Fig 10)
U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0; rho_u = 0.001;
l_c   = U^2 / g;
t_c   = U / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
gamma_rho = 1 / (1 + rho_r);
discriminant = (1 + rho_r)^2 - 4*alpha*(1 - rho_r);
k_l = ((1 + rho_r) + sqrt(discriminant)) / (2*alpha);
k_s = ((1 + rho_r) - sqrt(discriminant)) / (2*alpha);
F0  = 0.01 * T / (rho_l * U^2 * l_c);
epsilon_pv = 1e-6; AbsTol = 1e-10; RelTol = 1e-8;
chi = @(k) sqrt((1-rho_r)/(1+rho_r)*k + gamma_rho*alpha*k.^3);

x_grid = linspace(-15, 15, 2001); x_grid(abs(x_grid) < 1e-12) = [];
t_dim  = 25;
t      = t_dim / (100 * t_c);

% IVP η — ArrayValued combined-integrand inversion (same as Fig 8)
combined = @(k) ...
    2*cos(k*x_grid) ./ (alpha*(k - k_l).*(k - k_s)) ...
    - (1+rho_r)*(k + chi(k)).*cos(k*(t - x_grid) - t*chi(k)) ./ ...
      ((1 - rho_r + alpha*k^2)*alpha.*(k - k_l).*(k - k_s)) ...
    - (1+rho_r)*(k - chi(k)).*cos(k*(t - x_grid) + t*chi(k)) ./ ...
      ((1 - rho_r + alpha*k^2)*alpha.*(k - k_l).*(k - k_s));

I_lo = integral(combined, 0, k_s - epsilon_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
I_mi = integral(combined, k_s + epsilon_pv, k_l - epsilon_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
I_hi = integral(combined, k_l + epsilon_pv, Inf, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_sim = -F0/(2*pi) * (I_lo + I_mi + I_hi);

% Basilisk simulation data
data  = readmatrix('notebooks/if_25.csv');
data  = sortrows(data, 6);
x_bsk = data(:,6) / l_c;
y_bsk = data(:,7) / l_c;

figure; hold on;
plot(x_grid, eta_sim*1e3, 'b--', 'LineWidth', 3);
plot(x_bsk,  y_bsk*1e3,  'r:',  'LineWidth', 3);
xlabel('x'); ylabel('\eta \times 10^3');
legend('\eta (IVP)', 'Simulation'); xlim([-6 10]);
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/cg_sim_fig10.png" alt="Fig 10(i)" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 10(i).</strong> IVP vs nonlinear simulation (Basilisk, Navier–Stokes/VOF) at <em>t</em><sub>dim</sub> = 25 s (Julia).</figcaption>
</figure>
```

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/fig10_overlay.png" alt="Fig 10(ii) overlay" style="max-width:80%; height:auto;">
  <figcaption style="text-align:center;"><strong>Fig. 10(ii).</strong> Julia (lines), MATLAB (markers), and Basilisk simulation overlay, demonstrating the two computations are equivalent.</figcaption>
</figure>
```
