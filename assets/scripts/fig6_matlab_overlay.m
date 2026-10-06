%% Fig. 7(ii) — Julia vs MATLAB overlay for -I4 at t = 0.37
clear; close all; clc;

julia = readmatrix('julia_fig7_i4.csv');
x_j = julia(:,1);
y_j = julia(:,2);

U = 26.7046; g = 981.0; T = 72.0;
rho_l = 1.0; rho_u = 0.001;
l_c = U^2 / g;
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;
beta = (1.0 - rho_r) / (1.0 + rho_r);
gamma_rho = 1.0 / (1.0 + rho_r);
t = 0.37;
AbsTol = 1e-10; RelTol = 1e-8;

x = linspace(-15, 15, 2001);
x(abs(x) < 1e-12) = [];
assert(numel(x) == numel(x_j) && max(abs(x(:) - x_j(:))) < 1e-12, ...
       'Julia and MATLAB grids do not match.');
chi = @(k) sqrt(beta*k + gamma_rho*alpha*k.^3);
I4fun = @(k) (k == 0) * zeros(size(x)) + (k ~= 0) * ...
    (k .* cos(t*(k + chi(k)) - k*x) ./ ((k + chi(k)) .* (1 + alpha*k.^2 - rho_r)));
y_m = -integral(I4fun, 0, Inf, 'ArrayValued', true, ...
                'AbsTol', AbsTol, 'RelTol', RelTol);

fprintf('max |-I4_Julia - I4_MATLAB| = %.6e\n', max(abs(y_j(:) - y_m(:))));
idx = 1:25:numel(x);
figure('Color','w', 'Position', [100 100 460 270]); hold on; box on;
h_j = plot(x_j, y_j, '-', 'Color', [0.5 0 0.5], 'LineWidth', 2);
h_m = plot(x(idx), y_m(idx), 'o', 'LineStyle', 'none', ...
           'Color', [0.95 0.60 0.05], 'MarkerSize', 6, ...
           'MarkerFaceColor', 'none', 'LineWidth', 0.7);
xlim([-10 10]); ylim([-2 14]); yticks([0 3 6 9 12]);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('−𝕀₄', 'Interpreter', 'none', 'FontName', 'DejaVu Serif', 'FontSize', 20);
legend([h_j h_m], {'Julia', 'MATLAB'}, 'Interpreter','latex', ...
       'FontSize', 16, 'Location','northeast', 'Box','off');
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, '../fig7_overlay.png', 'Resolution', 300);
savefig(gcf, 'fig7_overlay.fig');
disp('Saved fig7_overlay.png and fig7_overlay.fig');
