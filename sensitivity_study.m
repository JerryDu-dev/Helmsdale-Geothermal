%% === Helmsdale Systematic Parameter Sweep (Revised & Calibrated) ===
% This script evaluates the sensitivity of the 100C isotherm depth 
% to Fault Permeability (KD) and Granite Radiogenic Heat (Qr).
% Optimized for the 22 K/km Basal Gradient scenario.

clear; clc; close all;

% 1. Shared Numerical Configuration
Lx = 14000; Ly = 6000; Nx = 141; Ny = 61;
x = linspace(0, Lx, Nx); y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y); dx = x(2)-x(1); dy = y(2)-y(1);

% %% === CRITICAL CHANGE 1: Increase simulation time for steady state ===
years = 150000; % Increased to 300 kyr to ensure full thermal equilibrium
sec_per_year = 365.25 * 24 * 3600;
dt_mult = 0.1;

% --- Define Parameter Space ---
KD_vec = logspace(-9, -6, 7); 
depth_KD_fault = zeros(size(KD_vec));

Qr_vec = linspace(3e-6, 9e-6, 7);
depth_Qr_granite = zeros(size(Qr_vec));

KD_base = 2.5e-7;
Qr_base = 6.53e-6;

fprintf('Starting Systematic Parameter Sweep (22 K/km Baseline)...\n');

%% --- SWEEP 1: Fault Permeability (KD) Variation ---
fprintf('\n--- Sweep 1: KD Variation (Measuring at Fault x=5.1km) ---\n');
for i = 1:length(KD_vec)
    fprintf('Running KD = %.1e... ', KD_vec(i));
    [d_fault, ~] = run_model(X, Y, dx, dy, Ny, Nx, years, sec_per_year, dt_mult, y, Qr_base, KD_vec(i));
    depth_KD_fault(i) = d_fault;
    fprintf('Depth: %.2f km\n', depth_KD_fault(i));
end

%% --- SWEEP 2: Radiogenic Heat (Qr) Variation ---
fprintf('\n--- Sweep 2: Qr Variation (Measuring in Granite x=2.0km) ---\n');
for i = 1:length(Qr_vec)
    fprintf('Running Qr = %.2e... ', Qr_vec(i));
    [~, d_granite] = run_model(X, Y, dx, dy, Ny, Nx, years, sec_per_year, dt_mult, y, Qr_vec(i), KD_base);
    depth_Qr_granite(i) = d_granite;
    fprintf('Depth: %.2f km\n', depth_Qr_granite(i));
end

%% === Visualization: Professional Sensitivity Curves ===
figure('Position', [100, 100, 1000, 450], 'Color', 'w');

subplot(1, 2, 1);
semilogx(KD_vec, depth_KD_fault, '-ko', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
set(gca, 'YDir', 'reverse'); 
grid on;
xlabel('Fault Permeability Factor K_D (m^2/Pa\cdot s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('100\circC Depth at Fault (km)', 'FontSize', 12, 'FontWeight', 'bold');
title('A. Sensitivity to Fault Permeability', 'FontSize', 14);
ylim([1, 4.5]); % Consistent Y-axis

subplot(1, 2, 2);
plot(Qr_vec*1e6, depth_Qr_granite, '-ko', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
set(gca, 'YDir', 'reverse'); 
grid on;
xlabel('Granite Radiogenic Heat Q_r (\muW/m^3)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('100\circC Depth in Granite (km)', 'FontSize', 12, 'FontWeight', 'bold');
title('B. Sensitivity to Radiogenic Heat', 'FontSize', 14);
ylim([1, 4.5]); 

%% === Internal Function: Core FDM Solver ===
function [depth_fault, depth_granite] = run_model(X, Y, dx, dy, Ny, Nx, years, sec_per_year, dt_mult, y_vec, Qr_val, KD_val)
    % Property Matrices initialization (Basement defaults)
    K_mat = ones(Ny, Nx)*2.5; Rho_mat = ones(Ny, Nx)*2700;
    Cp_mat = ones(Ny, Nx)*900; Qr_mat = ones(Ny, Nx)*1.0e-6; 
    KD_mat = ones(Ny, Nx)*1.0e-11; Type_mat = zeros(Ny, Nx);
    
    % 1. Granite
    mask_granite = (X < 5000); 
    K_mat(mask_granite) = 2.67; Rho_mat(mask_granite) = 2630; Cp_mat(mask_granite) = 836;
    Qr_mat(mask_granite) = Qr_val; Type_mat(mask_granite) = 1;
    
    % 2. Fault
    fault_center_x = 5000 + Y / tand(60); 
    mask_fault = abs(X - fault_center_x) < 250;
    K_mat(mask_fault) = 2.70; Rho_mat(mask_fault) = 2299; Cp_mat(mask_fault) = 1031;
    KD_mat(mask_fault) = KD_val; Type_mat(mask_fault) = 3;   
    
    % 3. Basin
    mask_basin = (X > fault_center_x + 250) & (Y < 2500);
    K_mat(mask_basin) = 1.78; Rho_mat(mask_basin) = 2073; Cp_mat(mask_basin) = 1361;
    Qr_mat(mask_basin) = 0.5e-6; Type_mat(mask_basin) = 2;
    
    % Flow Proxy
    Vy = -2.0e-2 * KD_mat; 
    Vy(Type_mat ~= 3) = 0; 
    Vy = smoothdata(Vy, 1, 'gaussian', 5);
    
    % Solver Init
    D_mat = K_mat ./ (Rho_mat .* Cp_mat);
    dt_diff = dt_mult * min(dx^2, dy^2) / max(D_mat(:));
    max_v = max(abs(Vy(:)));
    if max_v > 0, dt = min(dt_diff, dt_mult * dy / max_v); else, dt = dt_diff; end
    nt = ceil(years * sec_per_year / dt);
    
    % %% === CRITICAL CHANGE 2: Align Initial Condition ===
    T_surf = 10;
    DT_dy_init = 0.022; % Changed from 0.030 to 0.022
    T = T_surf + DT_dy_init * Y; 
    T_new = T;
    
    Const_Adv = (1000 * 4200) ./ (Rho_mat .* Cp_mat);
    
    for n = 1:nt
        k_right = 0.5*(K_mat + circshift(K_mat, [0,-1])); 
        k_left  = 0.5*(K_mat + circshift(K_mat, [0, 1]));
        dT2_dx2 = (k_right.*(circshift(T,[0,-1])-T) - k_left.*(T-circshift(T,[0,1]))) / dx^2;
        k_down  = 0.5*(K_mat + circshift(K_mat, [-1,0])); 
        k_up    = 0.5*(K_mat + circshift(K_mat, [1, 0]));
        dT2_dy2 = (k_down.*(circshift(T,[-1,0])-T) - k_up.*(T-circshift(T,[1,0]))) / dy^2;
        Diff_Term = (1./(Rho_mat.*Cp_mat)) .* (dT2_dx2 + dT2_dy2);
        
        dT_dy_up = zeros(size(T));
        T_shift = circshift(T, [-1, 0]); 
        mask_up = Vy < 0; 
        dT_dy_up(mask_up) = (T_shift(mask_up) - T(mask_up)) / dy; 
        Adv_Term = Const_Adv .* (Vy .* dT_dy_up);
        
        Source_Term = Qr_mat ./ (Rho_mat .* Cp_mat);
        T_new = T + dt * (Diff_Term - Adv_Term + Source_Term);
        
        T_new(1, :) = T_surf; 
        % %% === CRITICAL CHANGE 3: Align Basal Boundary ===
        T_new(end, :) = T_new(end-1, :) + 0.022 * dy; % Changed from 0.035 to 0.022
        T_new(:, 1) = T_new(:, 2); T_new(:, end) = T_new(:, end-1);
        T = T_new;
    end
    
    % Extract depths
    idx_f = round(5100 / dx); prof_f = T(:, idx_f);
    depth_fault = interp1(prof_f, y_vec, 100, 'linear') / 1000;
    
    idx_g = round(2000 / dx); prof_g = T(:, idx_g);
    depth_granite = interp1(prof_g, y_vec, 100, 'linear') / 1000;
end