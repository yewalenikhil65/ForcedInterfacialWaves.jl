%% Fig10_jfm_vinod.m
%
% Two-fluid capillary-gravity IVP solution.
%
% The time-dependent IVP surface elevation is obtained by a Cauchy
% principal value (CPV) integral around the two real poles k_s (gravity
% root) and k_l (capillary root) of the two-fluid dispersion relation.
% At the final time the classical steady-state solution is also
% evaluated for comparison.
%
% Optional: if simulation snapshots are available at
%   fullfile(pwd, 'matlab', 'fig10_simulation_data', 'if_<time_index>.csv')
% they are overlaid on the theoretical curves. Snapshots are skipped
% silently (with a warning) if not found.

clear;
clc;
close all;

%% Spatial resolution
Nx_plot = 2001;

%% Dimensional physical parameters (CGS units)
U     = 26.7046;    % Uniform base flow speed [cm/s]
g     = 981.0;       % Gravitational acceleration [cm/s^2]
T     = 72.0;         % Surface tension [dyn/cm]
rho_l = 1.0;          % Lower-fluid density [g/cm^3]
rho_u = 0.001;        % Upper-fluid density [g/cm^3]
L     = 75.5996;      % Dimensional domain length [cm]

%% Characteristic scales
l_c = U^2 / g;
t_c = U / g;
F_c = rho_l * U^2 * l_c;

%% Nondimensional parameters
alpha = T / (rho_l * U^2 * l_c);
rho_r = rho_u / rho_l;

beta      = (1.0 - rho_r) / (1.0 + rho_r);
gamma_rho = 1.0 / (1.0 + rho_r);

%% Gravity and capillary wave roots
discriminant = (1.0 + rho_r)^2 - 4.0 * alpha * (1.0 - rho_r);

if discriminant <= 0.0
    error('The steady capillary-gravity roots are not distinct positive real numbers.');
end

k_l = ((1.0 + rho_r) + sqrt(discriminant)) / (2.0 * alpha);
k_s = ((1.0 + rho_r) - sqrt(discriminant)) / (2.0 * alpha);

%% Dimensional wavenumbers and wavelengths
k_l_dim = k_l / l_c;
k_s_dim = k_s / l_c;

lambda_c = 2.0 * pi / k_l_dim;  % Capillary wavelength
lambda_g = 2.0 * pi / k_s_dim;  % Gravity wavelength

fprintf('alpha    = %.2e\n', alpha);
fprintf('rho_r    = %.2e\n', rho_r);
fprintf('k_s      = %.2e\n', k_s);
fprintf('k_l      = %.2e\n', k_l);
fprintf('lambda_c = %.2e cm\n', lambda_c);
fprintf('lambda_g = %.2e cm\n', lambda_g);

%% Nondimensional forcing amplitude
F0_dim = 0.01 * T;
F0 = F0_dim / F_c;

%% Nondimensional spatial grid
x_grid = linspace(-L / 2.0, L / 2.0, Nx_plot) / l_c;

% Remove x = 0 from the grid
x_grid(abs(x_grid) < eps) = [];

Nx = numel(x_grid);

%% Quadrature parameters
epsilon_pv = 1.0e-6;
AbsTol     = 1.0e-10;
RelTol     = 1.0e-8;

k_max        = Inf;
k_max_steady = Inf;

%% Optional simulation data directory (leave empty to disable overlay)
simulation_data_directory = fullfile(pwd, 'matlab', 'fig10_simulation_data');

%% Archived two-column (x, eta) CSV comparison at time_index = 300
comparison_file = fullfile(pwd, 'docs', 'src', 'assets', ...
    'matlab_cg_ivp_t300.csv');
have_comparison = exist(comparison_file, 'file') == 2;

if have_comparison
    comparison_matrix = readmatrix(comparison_file);
    if size(comparison_matrix, 2) < 2
        error('Expected at least two columns in comparison CSV: %s', comparison_file);
    end
    x_comparison   = comparison_matrix(:, 1);
    eta_comparison = comparison_matrix(:, 2);
    fprintf('CSV comparison: %s\n', comparison_file);
else
    warning('CSV comparison not found (skipping overlay): %s', comparison_file);
end

%% Output directory
output_directory = fullfile(pwd, 'output', 'figure10');

if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

%% Simulation times, stored as hundredths of a second
time_indices = [1, 3, 7, 15, 25, 60, 145, 300];

for it = 1:numel(time_indices)

    time_index = time_indices(it);

    % Dimensional and nondimensional times
    t_dim = time_index / 100.0;
    t = t_dim / t_c;

    %% Read optional simulation snapshot
    have_simulation = false;
    simulation_file = fullfile(simulation_data_directory, ...
        sprintf('if_%d.csv', time_index));

    if exist(simulation_file, 'file')
        simulation_matrix = readmatrix(simulation_file, 'NumHeaderLines', 1);

        % Columns 6 and 7 contain x and eta
        simulation_matrix = sortrows(simulation_matrix(:, [6, 7]), 1);

        x_sim   = simulation_matrix(:, 1) / l_c;
        eta_sim = simulation_matrix(:, 2) / l_c;

        have_simulation = true;
    else
        warning('Simulation file not found (skipping overlay): %s', simulation_file);
    end

    %% Allocate theoretical solutions
    eta_ivp = zeros(size(x_grid));

    if time_index == 300
        eta_steady = nan(size(x_grid));
    else
        eta_steady = [];
    end

    %% Evaluate the theoretical solution at each grid point
    for ix = 1:Nx

        x = x_grid(ix);

        % Two-fluid dispersion function
        chi = @(k) sqrt(beta .* k + gamma_rho * alpha .* k.^3);

        % Nominally time-independent contribution
        integrand_eta_s = @(k) 2.0 .* cos(k .* x) ./ (alpha .* (k - k_l) .* (k - k_s));

        % First time-dependent contribution
        integrand_I3 = @(k) ...
            -((1.0 + rho_r) .* (k + chi(k))) ./ (1.0 - rho_r + alpha .* k.^2) .* ...
            cos(k .* (t - x) - t .* chi(k)) ./ (alpha .* (k - k_l) .* (k - k_s));

        % Second time-dependent contribution
        integrand_I4 = @(k) ...
            -((1.0 + rho_r) .* (k - chi(k))) ./ (1.0 - rho_r + alpha .* k.^2) .* ...
            cos(k .* (t - x) + t .* chi(k)) ./ (alpha .* (k - k_l) .* (k - k_s));

        % Combine the terms so removable singularities at k_s, k_l cancel
        % before quadrature.
        total_integrand = @(k) integrand_eta_s(k) + integrand_I3(k) + integrand_I4(k);

        % Integration below k_s
        integral_1 = integral(total_integrand, 0.0, k_s - epsilon_pv, ...
            'AbsTol', AbsTol, 'RelTol', RelTol);

        % Integration between k_s and k_l
        integral_2 = integral(total_integrand, k_s + epsilon_pv, k_l - epsilon_pv, ...
            'AbsTol', AbsTol, 'RelTol', RelTol);

        % Integration above k_l
        integral_3 = integral(total_integrand, k_l + epsilon_pv, k_max, ...
            'AbsTol', AbsTol, 'RelTol', RelTol);

        % Time-dependent IVP solution
        eta_ivp(ix) = -F0 / (2.0 * pi) * (integral_1 + integral_2 + integral_3);

        %% Classical steady solution (evaluated only at the final time)
        if time_index == 300

            G_x_combined_integrand = @(k) cos(k .* x) ./ (k + k_s) - cos(k .* x) ./ (k + k_l);

            G_x = 1.0 / (k_l - k_s) * integral(G_x_combined_integrand, 0.0, k_max_steady, ...
                'AbsTol', AbsTol, 'RelTol', RelTol);

            if x > 0.0
                eta_steady(ix) = -2.0 * F0 / (alpha * (k_l - k_s)) * sin(k_s * x) ...
                    + F0 / (pi * alpha) * G_x;
            else
                eta_steady(ix) = -2.0 * F0 / (alpha * (k_l - k_s)) * sin(k_l * x) ...
                    + F0 / (pi * alpha) * G_x;
            end
        end
    end

    %% Publication-quality figure settings
    figure_size_in = [3.45, 2.30];  % width x height in inches

    axis_line_width  = 0.9;
    curve_line_width = 1.8;

    axis_font_size  = 10;
    label_font_size = 12;
    legend_font_size = 10;

    fig = figure('Visible', 'on', 'Color', 'w', 'Renderer', 'painters');
    set(fig, 'Units', 'inches', 'Position', [1, 1, figure_size_in]);

    layout = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax = nexttile(layout);
    hold(ax, 'on');

    set(ax, 'TickDir', 'out', 'LineWidth', axis_line_width, ...
        'FontSize', axis_font_size, 'XAxisLocation', 'bottom', ...
        'YAxisLocation', 'left');

    %% Plot steady, simulation (if available), and IVP curves
    if time_index == 300

        % 1. Steady solution
        h1 = plot(ax, x_grid, 1.0e3 * eta_steady, 'k-', 'LineWidth', curve_line_width);

        handles = h1;
        labels  = {'Steady'};

        % 2. Simulation overlay
        if have_simulation
            h2 = plot(ax, x_sim, 1.0e3 * eta_sim, 'r:', 'LineWidth', curve_line_width);
            handles = [handles, h2];
            labels  = [labels, {'Simulation'}];
        end

        % 3. Archived CSV comparison
        if have_comparison
            h3 = plot(ax, x_comparison, 1.0e3 * eta_comparison, 'r:', ...
                'LineWidth', curve_line_width);
            handles = [handles, h3];
            labels  = [labels, {'CSV reference'}];
        end

        % 4. Theoretical IVP solution
        h4 = plot(ax, x_grid, 1.0e3 * eta_ivp, 'b--', 'LineWidth', curve_line_width);
        handles = [handles, h4];
        labels  = [labels, {'Theory'}];

        lgd = legend(ax, handles, labels, 'Interpreter', 'latex', ...
            'FontSize', legend_font_size, 'Location', 'northeast');

    else

        handles = [];
        labels  = {};

        % Simulation overlay
        if have_simulation
            h2 = plot(ax, x_sim, 1.0e3 * eta_sim, 'r:', 'LineWidth', curve_line_width);
            handles = [handles, h2];
            labels  = [labels, {'Simulation'}];
        end

        % Theoretical IVP solution
        h3 = plot(ax, x_grid, 1.0e3 * eta_ivp, 'b-.', 'LineWidth', curve_line_width);
        handles = [handles, h3];
        labels  = [labels, {'Theory'}];

        lgd = legend(ax, handles, labels, 'Interpreter', 'latex', ...
            'FontSize', legend_font_size, 'Location', 'northeast');
    end

    % Shorter line samples inside legend
    lgd.ItemTokenSize = [12, 9];

    %% Axis labels and limits
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', label_font_size);
    ylabel(ax, '$10^3\eta$', 'Interpreter', 'latex', 'FontSize', label_font_size);

    xlim(ax, [-5, 10]);
    ylim(ax, [-6, 13]);
    yticks(ax, [-4, 0, 4, 8]);

    %% Draw the top and right frame lines (box off, manual frame on)
    set(ax, 'Box', 'off');
    delete(findall(ax, 'Tag', 'TopRightFrame'));

    x_limits = xlim(ax);
    y_limits = ylim(ax);

    line(ax, [x_limits(1), x_limits(2)], [y_limits(2), y_limits(2)], ...
        'Color', 'k', 'LineWidth', axis_line_width, 'Clipping', 'off', ...
        'HandleVisibility', 'off', 'Tag', 'TopRightFrame');

    line(ax, [x_limits(2), x_limits(2)], [y_limits(1), y_limits(2)], ...
        'Color', 'k', 'LineWidth', axis_line_width, 'Clipping', 'off', ...
        'HandleVisibility', 'off', 'Tag', 'TopRightFrame');

    %% Save figure
    output_name = sprintf('figure10_t_%04d', time_index);

    exportgraphics(fig, fullfile(output_directory, [output_name, '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');

    exportgraphics(fig, fullfile(output_directory, [output_name, '.png']), ...
        'Resolution', 900, 'BackgroundColor', 'white');
end

fprintf('Figures saved in: %s\n', output_directory);
