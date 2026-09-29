%% Fig. 5a(ii) — Julia vs MATLAB overlay, steady-state WITHOUT Rayleigh dissipation (eqn. 3.11)
% Plots eta_s^{local}(x) and eta_s^{far-field}(x) [notation from the manuscript/docs text]
%
% Workflow:
%   1. Load the Julia-computed (x, eta) datasets (saved as CSV alongside this script).
%   2. Plot the Julia curves first.
%   3. Compute the same quantities independently in MATLAB and overlay with `hold on`.
%   4. All labels/legend use LaTeX interpreter — edit interactively in the MATLAB GUI
%      (Edit Plot mode) if further tweaks are needed, then File > Export/Save as needed.

clear; close all; clc;

%% ─── Load Julia reference dataset ───
julia_local    = readmatrix('julia_ssl_no_rayleigh_local.csv');
julia_farfield = readmatrix('julia_ssl_no_rayleigh_farfield.csv');

x_j        = julia_local(:,1);
eta_loc_j  = julia_local(:,2);
eta_ff_j   = julia_farfield(:,2);

%% ─── MATLAB computation (independent) — eqn. (3.11) ───
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
eta_far_m = F0/(alpha*(k_l - k_s)) .* (-sin(k_s*abs(x)) + sin(k_l*abs(x)));

% Local: exponentially-decaying integral
eta_loc_m = zeros(size(x));
for i = 1:length(x)
    eta_loc_m(i) = F0*(k_l + k_s)/(pi*alpha) * ...
        integral(@(y) y.*exp(-abs(x(i)).*y)./((y.^2+k_l^2).*(y.^2+k_s^2)), 0, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-8);
end

%% ─── Downsample MATLAB curve to markers (denser sampling for far-field) ───
[pk_val, pk_idx]   = findpeaks(eta_far_m);
[tr_val, tr_idx]   = findpeaks(-eta_far_m);
zero_idx   = find(diff(sign(eta_far_m)) ~= 0);
stride_idx = 1:6:numel(x);                          % regular dense sampling
marker_idx = sort(unique([pk_idx, tr_idx, zero_idx, stride_idx]));

%% ─── Plot: Julia lines first, then MATLAB overlay via hold on ───
figure('Color', 'w');
hold on; box on;

% Julia curves (solid lines)
plot(x_j, eta_loc_j*1e3, '-', 'Color', [0.85 0 0], 'LineWidth', 3);
plot(x_j, eta_ff_j*1e3,  '-', 'Color', [0 0 0.85], 'LineWidth', 3);

% MATLAB overlay (scatter markers, high-contrast colors distinct from the Julia lines)
scatter(x(marker_idx), eta_far_m(marker_idx)*1e3, 45, [0.95 0.60 0.05], 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9);
scatter(x(1:25:end), eta_loc_m(1:25:end)*1e3, 45, [0 0.55 0.30], 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9);

xlim([-10 10]);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\eta \times 10^{3}$', 'Interpreter', 'latex', 'FontSize', 18);
legend({'$\eta_s^{\mathrm{local}}$ (Julia)', '$\eta_s^{\mathrm{far\mbox{-}field}}$ (Julia)', ...
        '$\eta_s^{\mathrm{far\mbox{-}field}}$ (MATLAB)', '$\eta_s^{\mathrm{local}}$ (MATLAB)'}, ...
       'Interpreter', 'latex', 'FontSize', 13, 'Location', 'eastoutside', 'Box', 'off');

%% ─── Save figure ───
% Edit interactively in the GUI (Edit Plot / Property Inspector) if needed, then re-run
% the two lines below, or use File > Export Setup for further control.
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, '../fig5a_overlay.png', 'Resolution', 300);
savefig(gcf, 'fig5a_overlay.fig');   % editable MATLAB figure for later GUI tweaks
disp('Saved fig5a_overlay.png and fig5a_overlay.fig');
