%% === Helmsdale Sensitivity Study ===
clear; clc; close all;

% 1. Shared Configuration
Lx = 14000; Ly = 6000; Nx = 141; Ny = 61;
x = linspace(0, Lx, Nx); y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
dx = x(2)-x(1); dy = y(2)-y(1);
years = 150000; sec_per_year = 365.25 * 24 * 3600; dt_mult = 0.1;

% Define 4 Scenarios: [Qr_granite, KD_fault, Scenario_Name]
scenarios = {
    6.53e-6, 2.5e-7, 'Baseline';
    4.00e-6, 2.5e-7, 'Low Radiogenic Heat (Qr)';
    6.53e-6, 1.0e-8, 'Low Fault Permeability (KD)';
    6.53e-6, 1.0e-6, 'High Fault Permeability (KD)'
};

% Storage for results
T_results = cell(4, 1);
Depth_100C_Fault = zeros(4, 1);
idx_fault = round(5100 / dx); % Fault location index

fprintf('Starting Sensitivity Study (4 Scenarios)...\n');

for s = 1:4
    fprintf('Running Scenario %d: %s...\n', s, scenarios{s, 3});
    
    % --- Base Matrices ---
    K_mat = ones(Ny, Nx)*2.5; Rho_mat = ones(Ny, Nx)*2700;
    Cp_mat = ones(Ny, Nx)*900; Qr_mat = ones(Ny, Nx)*1.0e-6; 
    KD_mat = ones(Ny, Nx)*1.0e-11; Type_mat = zeros(Ny, Nx);
    
    % --- Apply Geometry & Scenario Parameters ---
    % Granite
    mask_granite = (X < 5000); 
    K_mat(mask_granite) = 2.67; Rho_mat(mask_granite) = 2630; Cp_mat(mask_granite) = 836;
    Qr_mat(mask_granite) = scenarios{s, 1}; % <--- INJECT Qr
    Type_mat(mask_granite) = 1;
    
    % Fault
    fault_center_x = 5000 + Y / tand(60); mask_fault = abs(X - fault_center_x) < 250;
    K_mat(mask_fault) = 2.70; Rho_mat(mask_fault) = 2299; Cp_mat(mask_fault) = 1031;
    KD_mat(mask_fault) = scenarios{s, 2};   % <--- INJECT KD
    Type_mat(mask_fault) = 3;
    
    % Basin
    mask_basin = (X > fault_center_x + 250) & (Y < 2500);
    K_mat(mask_basin) = 1.78; Rho_mat(mask_basin) = 2073; Cp_mat(mask_basin) = 1361;
    Qr_mat(mask_basin) = 0.5e-6; Type_mat(mask_basin) = 2;
    
    % --- Fluid Flow ---
    Vy = -2.0e-2 * KD_mat; Vy(Type_mat ~= 3) = 0;
    Vy = smoothdata(Vy, 1, 'gaussian', 5);
    
    % --- Solver Init ---
    D_mat = K_mat ./ (Rho_mat .* Cp_mat);
    dt_diff = dt_mult * min(dx^2, dy^2) / max(D_mat(:));
    max_v = max(abs(Vy(:)));
    if max_v > 0, dt = min(dt_diff, dt_mult * dy / max_v); else, dt = dt_diff; end
    nt = ceil(years * sec_per_year / dt);
    
    T = 10 + 0.030 * Y; T_new = T;
    Const_Adv = (1000 * 4200) ./ (Rho_mat .* Cp_mat);
    
    % --- Main Loop ---
    for n = 1:nt
        k_right = 0.5*(K_mat + circshift(K_mat, [0,-1])); k_left = 0.5*(K_mat + circshift(K_mat, [0, 1]));
        dT2_dx2 = (k_right.*(circshift(T,[0,-1])-T) - k_left.*(T-circshift(T,[0,1]))) / dx^2;
        k_down = 0.5*(K_mat + circshift(K_mat, [-1,0])); k_up = 0.5*(K_mat + circshift(K_mat, [1, 0]));
        dT2_dy2 = (k_down.*(circshift(T,[-1,0])-T) - k_up.*(T-circshift(T,[1,0]))) / dy^2;
        Diff_Term = (1./(Rho_mat.*Cp_mat)) .* (dT2_dx2 + dT2_dy2);
        
        dT_dy_up = zeros(size(T)); dT_dy_up(Vy < 0) = (circshift(T(Vy < 0), [-1, 0]) - T(Vy < 0)) / dy; 
        Adv_Term = Const_Adv .* (Vy .* dT_dy_up);
        Source_Term = Qr_mat ./ (Rho_mat .* Cp_mat);
        
        T_new = T + dt * (Diff_Term - Adv_Term + Source_Term);
        T_new(1, :) = 10; T_new(end, :) = T_new(end-1, :) + 0.035 * dy; 
        T_new(:, 1) = T_new(:, 2); T_new(:, end) = T_new(:, end-1);    
        T = T_new;
    end
    
    % Store Results
    T_results{s} = T;
    
    % Find depth of 100C at the fault
    fault_temp_profile = T(:, idx_fault);
    depth_idx = find(fault_temp_profile >= 100, 1, 'first');
    if isempty(depth_idx)
        Depth_100C_Fault(s) = NaN; % Did not reach 100C
    else
        Depth_100C_Fault(s) = y(depth_idx) / 1000; % Convert to km
    end
    fprintf('  -> 100C Depth at Fault: %.2f km\n', Depth_100C_Fault(s));
end

%% === Visualization: Isotherm Shift Comparison ===
figure('Position', [150, 150, 900, 600], 'Color', 'w'); hold on;
% Background geological context (from Baseline)
imagesc(x/1000, y/1000, Type_mat); colormap([0.9 0.9 0.9; 1 0.8 0.8; 0.8 0.9 1; 1 0.8 1]); caxis([0 3]);

% Plot 100C Isotherms for all scenarios
colors = {'k', 'b', 'r', 'g'};
line_styles = {'-', '--', '-.', ':'};
lines = zeros(4,1);

for s = 1:4
    [C, h] = contour(x/1000, y/1000, T_results{s}, [100 100], ...
        'LineColor', colors{s}, 'LineStyle', line_styles{s}, 'LineWidth', 2.5);
    lines(s) = h;
end

% Formatting
axis ij; pbaspect([2.5 1 1]); 
xlim([0 Lx/1000]); ylim([0 Ly/1000]);
xlabel('Distance (km)', 'FontSize', 12); ylabel('Depth (km)', 'FontSize', 12);
title('Sensitivity Study: Shift of the 100\circC Isotherm', 'FontSize', 14);

% Custom Legend
legend(lines, ...
    sprintf('Baseline (100\\circC at %.2f km)', Depth_100C_Fault(1)), ...
    sprintf('Low Q_r (100\\circC at %.2f km)', Depth_100C_Fault(2)), ...
    sprintf('Low K_D (100\\circC at %.2f km)', Depth_100C_Fault(3)), ...
    sprintf('High K_D (100\\circC at %.2f km)', Depth_100C_Fault(4)), ...
    'Location', 'SouthEast', 'FontSize', 11);
grid on;