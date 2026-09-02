%% generate_cg_comparison_data.m
%  Generates MATLAB reference datasets for capillary-gravity Figs 7, 8, 10.
%  Uses the same parameter defaults and integrands as the Julia package.
%
%  Grid: 401 points linspace(-15,15) excluding |x|<1e-12
%  (reduced from 2001 to keep MATLAB IVP integrations tractable)
%
%  Writes to docs/src/assets/:
%    cg_I4_t034.csv          — -I4/(2*pi) at t=0.34   (Fig 7)
%    cg_ivp_eta_t36735.csv   — eta, eta_tr at t=367.35 (Fig 8)
%    cg_ivp_eta_t25s.csv     — eta at t_dim=25s        (Fig 10)
%
%  Usage:
%    matlab -batch "run('/full/path/matlab/generate_cg_comparison_data.m')"

clear; clc;

%% ── Parameters (identical to compute_cg_parameters() defaults) ───────────────
U      = 26.7046;
g      = 981.0;
T      = 72.0;
rho_l  = 1.0;
rho_u  = 0.001;
l_c    = U^2 / g;
t_c    = U / g;
alpha  = T / (rho_l * U^2 * l_c);
rho_r  = rho_u / rho_l;
disc   = (1 + rho_r)^2 - 4*alpha*(1 - rho_r);
k_l    = ((1 + rho_r) + sqrt(disc)) / (2*alpha);
k_s    = ((1 + rho_r) - sqrt(disc)) / (2*alpha);
F0     = 0.01 * T / (rho_l * U^2 * l_c);
AbsTol = 1e-8;   % slightly relaxed for oscillatory IVP integrals
RelTol = 1e-6;
epsilon_pv = 1e-6;
gamma_rho  = 1 / (1 + rho_r);

chi = @(k) sqrt(beta_val(rho_r)*k + gamma_rho*alpha*k.^3);

function b = beta_val(rho_r), b = (1-rho_r)/(1+rho_r); end

beta = (1 - rho_r) / (1 + rho_r);
chi  = @(k) sqrt(beta.*k + gamma_rho.*alpha.*k.^3);

fprintf('alpha=%.6e  k_s=%.6f  k_l=%.6f  F0=%.6e\n', alpha, k_s, k_l, F0);

%% ── Reduced spatial grid (401 pts) ──────────────────────────────────────────
Nx    = 401;
x_raw = linspace(-15, 15, Nx);
x     = x_raw(abs(x_raw) >= 1e-12);
N     = length(x);
fprintf('Grid: %d points\n', N);

outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', 'src', 'assets');

%% ══════════════════════════════════════════════════════════════════════════════
%% Fig 7 — I4 component at t = 0.34
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\nFig 7: I4 at t=0.34 (%d pts) ...\n', N);
t = 0.34;
I4 = zeros(1, N);
for i = 1:N
    xi = x(i);
    % Integrand for I4 from eqn (4.5d): uses (k - chi(k)) and cos(t*(k+chi(k))-k*x)
    f = @(k) -(1+rho_r)/alpha .* (k - chi(k)) .* cos(t.*(k + chi(k)) - k.*xi) ./ ...
              ((1 + alpha.*k.^2 - rho_r) .* (k - k_l) .* (k - k_s));
    I4(i) = integral(f, 0,           k_s-epsilon_pv, 'AbsTol', AbsTol, 'RelTol', RelTol) ...
          + integral(f, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol', AbsTol, 'RelTol', RelTol) ...
          + integral(f, k_l+epsilon_pv, Inf,           'AbsTol', AbsTol, 'RelTol', RelTol);
    if mod(i,80)==0, fprintf('  %d/%d\n',i,N); end
end
writematrix([x', (-(1/(2*pi))*I4)'], fullfile(outdir, 'cg_I4_t034.csv'));
fprintf('  Saved cg_I4_t034.csv\n');

%% ══════════════════════════════════════════════════════════════════════════════
%% Fig 8 — Full IVP (eta, eta_tr) at t = 367.35
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\nFig 8: IVP at t=367.35 (%d pts) ...\n', N);
t = 367.35;
eta_8   = zeros(1, N);
eta_s_8 = zeros(1, N);

for i = 1:N
    xi = x(i);
    combined = @(k) ...
         2.*cos(k.*xi)./(alpha.*(k-k_l).*(k-k_s)) ...
        -(1+rho_r).*(k+chi(k)).*cos(k.*(t-xi)-t.*chi(k)) ./ ...
         ((1-rho_r+alpha.*k.^2).*alpha.*(k-k_l).*(k-k_s)) ...
        -(1+rho_r).*(k-chi(k)).*cos(k.*(t-xi)+t.*chi(k)) ./ ...
         ((1-rho_r+alpha.*k.^2).*alpha.*(k-k_l).*(k-k_s));
    I_tot = integral(combined, 0,             k_s-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
          + integral(combined, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
          + integral(combined, k_l+epsilon_pv, Inf,            'AbsTol',AbsTol,'RelTol',RelTol);
    eta_8(i) = -F0/(2*pi) * I_tot;

    steady_f = @(k) cos(k.*xi)./(alpha.*(k-k_l).*(k-k_s));
    eta_s_8(i) = -F0/pi * ( ...
        integral(steady_f, 0,             k_s-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
      + integral(steady_f, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
      + integral(steady_f, k_l+epsilon_pv, Inf,            'AbsTol',AbsTol,'RelTol',RelTol));

    if mod(i,80)==0, fprintf('  %d/%d\n',i,N); end
end
eta_tr_8 = eta_8 - eta_s_8;
writematrix([x', eta_8', eta_tr_8'], fullfile(outdir, 'cg_ivp_eta_t36735.csv'));
fprintf('  Saved cg_ivp_eta_t36735.csv\n');

%% ══════════════════════════════════════════════════════════════════════════════
%% Fig 10 — IVP at t_dim = 25 s
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\nFig 10: IVP at t_dim=25s (%d pts) ...\n', N);
t = 25 / (100 * t_c);   % same factor as Julia: t_dim/(100*t_c)
eta_sim = zeros(1, N);

for i = 1:N
    xi = x(i);
    combined = @(k) ...
         2.*cos(k.*xi)./(alpha.*(k-k_l).*(k-k_s)) ...
        -(1+rho_r).*(k+chi(k)).*cos(k.*(t-xi)-t.*chi(k)) ./ ...
         ((1-rho_r+alpha.*k.^2).*alpha.*(k-k_l).*(k-k_s)) ...
        -(1+rho_r).*(k-chi(k)).*cos(k.*(t-xi)+t.*chi(k)) ./ ...
         ((1-rho_r+alpha.*k.^2).*alpha.*(k-k_l).*(k-k_s));
    eta_sim(i) = -F0/(2*pi) * ( ...
        integral(combined, 0,             k_s-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
      + integral(combined, k_s+epsilon_pv, k_l-epsilon_pv, 'AbsTol',AbsTol,'RelTol',RelTol) ...
      + integral(combined, k_l+epsilon_pv, Inf,            'AbsTol',AbsTol,'RelTol',RelTol));
    if mod(i,80)==0, fprintf('  %d/%d\n',i,N); end
end
writematrix([x', eta_sim'], fullfile(outdir, 'cg_ivp_eta_t25s.csv'));
fprintf('  Saved cg_ivp_eta_t25s.csv\n');

fprintf('\nAll done.\n');
