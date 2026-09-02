%% generate_steady_state_comparison_data.m
%  Computes the capillary-gravity steady-state profiles using parameters and
%  formulas identical to the Julia package defaults.
%
%  Local term uses the exponentially-decaying form (eqn. 3.11 in manuscript):
%    eta_local = F0*(k_l+k_s)/(pi*alpha) *
%                integral_0^inf  y*exp(-|x|*y) / ((y^2+k_l^2)*(y^2+k_s^2)) dy
%
%  This is mathematically equivalent to F0*G(x)/(pi*alpha) but is well-
%  conditioned for MATLAB's adaptive quadrature (exponential decay vs
%  oscillatory tails).  The Julia package's cg_Gx_integral uses the cosine
%  form internally; both expressions evaluate to identical values (proved in
%  docs/src/capillary_gravity_rayleigh_dissipation.md).
%
%  Far-field:
%    No Rayleigh: F0/(alpha*(k_l-k_s)) * (-sin(k_s*|x|) + sin(k_l*|x|))
%    Rayleigh:    -2*F0/(alpha*(k_l-k_s)) * sin(k_l*x)  for x<0
%                 -2*F0/(alpha*(k_l-k_s)) * sin(k_s*x)  for x>0
%
%  Grid: 2001 points linspace(-15,15), then remove the single point closest
%  to zero (matching make_cg_xgrid which filters |x| < 1e-12).
%
%  Writes to docs/src/assets/:
%    ssl_no_rayleigh_farfield.csv    ssl_no_rayleigh_local.csv
%    ssl_rayleigh_farfield.csv       ssl_rayleigh_local.csv
%
%  Usage (from repo root or anywhere):
%    matlab -batch "run('/full/path/matlab/generate_steady_state_comparison_data.m')"

clear; clc;

%% ── Parameters (matching compute_cg_parameters() defaults) ──────────────────
U     = 26.7046;
g     = 981.0;
T     = 72.0;
rho_l = 1.0;
rho_u = 0.001;

l_c   = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
disc  = (1 + rho_r)^2 - 4*alpha*(1 - rho_r);
k_l   = ((1 + rho_r) + sqrt(disc)) / (2*alpha);
k_s   = ((1 + rho_r) - sqrt(disc)) / (2*alpha);
F0    = 0.01 * T / (rho_l * U^2 * l_c);

AbsTol = 1e-10;
RelTol = 1e-8;

fprintf('alpha = %.15e\n', alpha);
fprintf('k_s   = %.15e\n', k_s);
fprintf('k_l   = %.15e\n', k_l);
fprintf('F0    = %.15e\n', F0);

%% ── Grid matching make_cg_xgrid(p; Nx=2001, xlim=(-15,15)) ──────────────────
Nx    = 2001;
x_raw = linspace(-15, 15, Nx);
x     = x_raw(abs(x_raw) >= 1e-12);   % removes the x=0 point
N     = length(x);
fprintf('Grid: %d points\n', N);

%% ── Local term (exponentially-decaying integrand, well-conditioned) ──────────
fprintf('Computing local term (%d points) ...\n', N);
eta_local = zeros(1, N);
prefactor = F0 * (k_l + k_s) / (pi * alpha);
for i = 1:N
    xi = x(i);
    integrand = @(y) y .* exp(-abs(xi) .* y) ./ ((y.^2 + k_l^2) .* (y.^2 + k_s^2));
    eta_local(i) = prefactor * integral(integrand, 0, Inf, 'AbsTol', AbsTol, 'RelTol', RelTol);
    if mod(i, 400) == 0, fprintf('  %d/%d\n', i, N); end
end
fprintf('  Done.\n');

%% ── Far-field: no Rayleigh (eqn. 3.11) ──────────────────────────────────────
denom     = alpha * (k_l - k_s);
eta_far_nr = F0 / denom .* (-sin(k_s * abs(x)) + sin(k_l * abs(x)));

%% ── Far-field: Rayleigh dissipation (eqn. 3.12) ─────────────────────────────
eta_far_r          = zeros(1, N);
eta_far_r(x > 0)   = -2*F0/denom * sin(k_s * x(x > 0));
eta_far_r(x < 0)   = -2*F0/denom * sin(k_l * x(x < 0));

%% ── Write CSVs ───────────────────────────────────────────────────────────────
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
outdir    = fullfile(repoRoot, 'docs', 'src', 'assets');

writematrix([x', eta_far_nr'],  fullfile(outdir, 'ssl_no_rayleigh_farfield.csv'));
writematrix([x', eta_local'],   fullfile(outdir, 'ssl_no_rayleigh_local.csv'));
writematrix([x', eta_far_r'],   fullfile(outdir, 'ssl_rayleigh_farfield.csv'));
writematrix([x', eta_local'],   fullfile(outdir, 'ssl_rayleigh_local.csv'));   % same local

fprintf('\nSaved four CSV files to: %s\n', outdir);
fprintf('Done.\n');
