%% Fig. 6(ii) — Julia vs MATLAB overlay, pure-gravity IVP at t = 183.68
% Workflow (as in Fig. 5a(ii) and Fig. 5b(ii)):
%   1. Load Julia profiles: x, eta, eta_CPV, eta_s, eta_tr.
%   2. Recompute all four profiles independently in MATLAB.
%   3. Plot Julia lines and overlay downsampled MATLAB markers.

clear; close all; clc;

%% ─── Load Julia reference profiles ───────────────────────────────────────────
julia = readmatrix('julia_fig6_profiles.csv');
x_j       = julia(:,1);
eta_j     = julia(:,2);
eta_cpv_j = julia(:,3);
eta_s_j   = julia(:,4);
eta_tr_j  = julia(:,5);

%% ─── MATLAB computation: eqns. (4.3), (4.4b–e), and direct CPV ──────────────
U = 26.7046; g = 981.0;
rho_l = 1.0; rho_u = 0.001;
rho_r = rho_u / rho_l;
beta = (1.0 - rho_r) / (1.0 + rho_r);
sqrt_beta = sqrt(beta);
F0 = 0.01 * 72.0 / (rho_l * U^2 * (U^2 / g));

eps_pv = 1e-6;
AbsTol = 1e-10;
RelTol = 1e-8;
K_MAX = 100.0;
t = 183.68;

x = linspace(-12, 12, 2001);
x(abs(x) < 1e-6) = [];
a = t - x;
s = sign(a);

assert(numel(x_j) == numel(x) && max(abs(x_j(:) - x(:))) < 1e-12, ...
       'Julia and MATLAB Figure 6 grids do not match.');

% Eqn. (4.3): steady contribution eta_s = F0*T0.
I0 = integral(@(y) exp(-y .* abs(x)) .* y ./ (beta^2 + y.^2), 0, Inf, ...
              'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
eta_s_m = F0 * (-pi * sin(beta * abs(x)) + I0) / (pi * (1 + rho_r));

% Eqns. (4.4b–e): eta_tr = F0*(T1 + T2 + T3 + T4).
T1 = -s .* sin(beta * x) / (1 + rho_r);
T2fun = @(v) v.^2 .* exp(-s .* 2 .* v.^2 .* a + s .* v * t * sqrt_beta) .* ...
    (sqrt_beta * cos(v * t * sqrt_beta) + ...
     s .* (2 * v - sqrt_beta) * sin(v * t * sqrt_beta)) ./ ...
    (beta + (2 * v - sqrt_beta).^2);
T2 = -4 / (pi * (1 + rho_r) * beta) * ...
     integral(T2fun, 0, K_MAX, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);

b = 0.5 * t * sqrt_beta;
X = b * sqrt(2 ./ (pi * abs(a)));
T3 = 1 / (pi * (1 + rho_r) * sqrt_beta) .* ...
     (1 + t ./ (2 * a)) .* sqrt(pi ./ (2 * abs(a))) .* ...
     (cos(b^2 ./ abs(a)) .* (0.5 - s .* fresnelc(X)) + ...
      sin(b^2 ./ abs(a)) .* (0.5 - s .* fresnels(X)));

T4fun = @(v) cos(v.^2 .* a + v * t * sqrt_beta) ./ (v + sqrt_beta);
T4 = -1 / (pi * (1 + rho_r)) * ...
     integral(T4fun, 0, K_MAX, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);

eta_tr_m = F0 * (T1 + T2 + T3 + T4);
eta_m = eta_s_m + eta_tr_m;

% Direct CPV check: split symmetrically about the removable pole k = beta.
G = @(k) cos(k * x) ./ (pi * (1 + rho_r) * (k - beta)) ...
       - k .* cos(k * (t - x) - t * sqrt(beta * k)) ./ ...
         (2 * pi * (1 - rho_r) * (k - sqrt(beta * k))) ...
       - k .* cos(k * (t - x) + t * sqrt(beta * k)) ./ ...
         (2 * pi * (1 - rho_r) * (k + sqrt(beta * k)));
Ilo = integral(G, 0, beta - eps_pv, 'ArrayValued', true, 'AbsTol', AbsTol, 'RelTol', RelTol);
Ihi = integral(G, beta + eps_pv, K_MAX, 'ArrayValued', true, ...
               'AbsTol', AbsTol, 'RelTol', RelTol);
eta_cpv_m = F0 * (Ilo + Ihi);

%% ─── Julia–MATLAB agreement ─────────────────────────────────────────────────
err_eta = max(abs(eta_j(:) - eta_m(:)));
err_cpv = max(abs(eta_cpv_j(:) - eta_cpv_m(:)));
err_s   = max(abs(eta_s_j(:) - eta_s_m(:)));
err_tr  = max(abs(eta_tr_j(:) - eta_tr_m(:)));
fprintf('max |eta_Julia - eta_MATLAB|         = %.6e\n', err_eta);
fprintf('max |eta_CPV_Julia - eta_CPV_MATLAB| = %.6e\n', err_cpv);
fprintf('max |eta_s_Julia - eta_s_MATLAB|     = %.6e\n', err_s);
fprintf('max |eta_tr_Julia - eta_tr_MATLAB|   = %.6e\n', err_tr);

%% ─── Plot: Julia lines first, MATLAB markers second ──────────────────────────
marker_idx = 1:50:numel(x);
c_eta = [0.00 0.35 0.95];
c_cpv = [0.85 0.10 0.10];
c_s = [0.00 0.00 0.00];
c_tr = [0.85 0.00 0.85];

figure('Color', 'w', 'Position', [100 100 1120 460]);
hold on; box on;

h_eta_j = plot(x_j, eta_j * 1e3, '-', 'Color', c_eta, 'LineWidth', 3);
h_cpv_j = plot(x_j, eta_cpv_j * 1e3, '-.', 'Color', c_cpv, 'LineWidth', 3);
h_s_j = plot(x_j, eta_s_j * 1e3, '-', 'Color', c_s, 'LineWidth', 3);
h_tr_j = plot(x_j, eta_tr_j * 1e3, ':', 'Color', c_tr, 'LineWidth', 3);

h_eta_m = scatter(x(marker_idx), eta_m(marker_idx) * 1e3, 70, c_eta, '*', ...
                  'LineWidth', 1.6);
h_cpv_m = scatter(x(marker_idx), eta_cpv_m(marker_idx) * 1e3, 38, c_cpv, 's', ...
                  'MarkerFaceColor', 'w', 'LineWidth', 1.2);
h_s_m = scatter(x(marker_idx), eta_s_m(marker_idx) * 1e3, 38, c_s, 'd', ...
                'MarkerFaceColor', 'w', 'LineWidth', 1.2);
h_tr_m = scatter(x(marker_idx), eta_tr_m(marker_idx) * 1e3, 38, c_tr, '^', ...
                 'MarkerFaceColor', 'w', 'LineWidth', 1.2);

xlim([-12 12]); ylim([-4.8 8.2]); yticks([-4 0 4 8]);
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'Box', 'on');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\eta \times 10^{3}$', 'Interpreter', 'latex', 'FontSize', 22);
legend([h_eta_j h_cpv_j h_s_j h_tr_j h_eta_m h_cpv_m h_s_m h_tr_m], ...
       {'$\eta$ (Julia)', '$\eta_{\mathrm{CPV}}$ (Julia)', ...
        '$\eta_s$ (Julia)', '$\eta_{tr}$ (Julia)', ...
        '$\eta$ (MATLAB)', '$\eta_{\mathrm{CPV}}$ (MATLAB)', ...
        '$\eta_s$ (MATLAB)', '$\eta_{tr}$ (MATLAB)'}, ...
       'Interpreter', 'latex', 'FontSize', 16, 'NumColumns', 2, ...
       'Location', 'eastoutside', 'Box', 'off');

%% ─── Save publication and editable assets ────────────────────────────────────
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, '../fig6_comparison.png', 'Resolution', 300);
print(gcf, 'fig6_comparison.svg', '-dsvg');
savefig(gcf, 'fig6_comparison.fig');
disp('Saved fig6_comparison.png, fig6_comparison.svg, and fig6_comparison.fig');
