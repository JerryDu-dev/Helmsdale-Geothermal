clear; clc; close all;

%% === 1. Model Configuration  ===
% Spatial Discretization
Lx = 14000;  % Width: 14 km
Ly = 6000;   % Depth: 6 km
Nx = 141;    % Grid size X
Ny = 61;     % Grid size Y

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
dx = x(2)-x(1); 
dy = y(2)-y(1);

% Time Stepping
years = 150000; % 150 kyr to reach quasi-steady state
sec_per_year = 365.25 * 24 * 3600;
dt_mult = 0.1;  % Safety factor for CFL

%% === 2. Material Properties from Data  ===
% Data Source: Helmsdale_Data.xlsx (Calculated Mean Values)

% --- A. Define Property Matrices (Default = Basement) ---
K_mat   = ones(Ny, Nx) * 2.5;  
Rho_mat = ones(Ny, Nx) * 2700;
Cp_mat  = ones(Ny, Nx) * 900;
Qr_mat  = ones(Ny, Nx) * 1.0e-6; 
KD_mat  = ones(Ny, Nx) * 1.0e-11; 
Type_mat = zeros(Ny, Nx);         % ID for plotting

% --- B. Define Geometry & Assign Data Values ---

% 1. Granite Intrusion (West) - He1/He2
mask_granite = (X < 5000); 
K_mat(mask_granite)   = 2.67;     % Data: Granite.csv Mean kT
Rho_mat(mask_granite) = 2630;     
Cp_mat(mask_granite)  = 836;      
Qr_mat(mask_granite)  = 6.53e-6;  % Data: High Radiogenic Heat
KD_mat(mask_granite)  = 1.0e-10;  
Type_mat(mask_granite) = 1;

% 2. Fault Zone (Center) - Fz
% Geometry: Dipping Normal Fault (~60 degrees)
fault_surf_x = 5000;
dip_angle = 60; 
fault_center_x = fault_surf_x + Y / tand(dip_angle);
half_width = 250; 
mask_fault = abs(X - fault_center_x) < half_width;

K_mat(mask_fault)   = 2.70;      % Data: Fault Zone.csv
Rho_mat(mask_fault) = 2299;
Cp_mat(mask_fault)  = 1031;
Qr_mat(mask_fault)  = 1.0e-6;    
KD_mat(mask_fault)  = 2.5e-7;    % Data: High Permeability
Type_mat(mask_fault) = 3;

% 3. Sedimentary Basin (East) - Ms
% Geometry: Hanging wall (Right of fault), depth ~2.5km
mask_basin = (X > fault_center_x + half_width) & (Y < 2500);
K_mat(mask_basin)   = 1.78;      % Data: Low k -> Thermal Blanket
Rho_mat(mask_basin) = 2073;
Cp_mat(mask_basin)  = 1361;
Qr_mat(mask_basin)  = 0.5e-6;    
KD_mat(mask_basin)  = 1.0e-8;    
Type_mat(mask_basin) = 2;

%% === 3. Fluid Flow Setup ===
Vx = zeros(Ny, Nx);
Vy = zeros(Ny, Nx);

% Flow driven by KD scaling (Darcy approximation)
V_base_scale = -2.0e-2; 
Vy = V_base_scale * KD_mat; 
Vy(Type_mat ~= 3) = 0; % Constrain to Fault
Vy = smoothdata(Vy, 1, 'gaussian', 5);

%% === 4. Solver Initialization  ===
D_mat = K_mat ./ (Rho_mat .* Cp_mat);

% CFL Check
dt_diff = dt_mult * min(dx^2, dy^2) / max(D_mat(:));
max_v = max(abs(Vy(:)));
if max_v > 0
    dt_adv = dt_mult * dy / max_v;
    dt = min(dt_diff, dt_adv);
else
    dt = dt_diff;
end
nt = ceil(years * sec_per_year / dt);

% Initial Condition
T_surf = 10;
DT_dy_init = 0.022; 
T = T_surf + DT_dy_init * Y;
T_new = T;

% Fluid properties
rho_f = 1000; Cp_f = 4200;
Const_Adv = (rho_f * Cp_f) ./ (Rho_mat .* Cp_mat);

%% === 5. Main Loop (FDM Solver) ===
fprintf('Starting Simulation... %d years, %d steps\n', years, nt);

for n = 1:nt
    % 1. Diffusion (Variable k)
    k_right = 0.5*(K_mat + circshift(K_mat, [0,-1]));
    k_left  = 0.5*(K_mat + circshift(K_mat, [0, 1]));
    dT2_dx2 = (k_right.*(circshift(T,[0,-1])-T) - k_left.*(T-circshift(T,[0,1]))) / dx^2;
    
    k_down  = 0.5*(K_mat + circshift(K_mat, [-1,0]));
    k_up    = 0.5*(K_mat + circshift(K_mat, [1, 0]));
    dT2_dy2 = (k_down.*(circshift(T,[-1,0])-T) - k_up.*(T-circshift(T,[1,0]))) / dy^2;
    
    Diff_Term = (1./(Rho_mat.*Cp_mat)) .* (dT2_dx2 + dT2_dy2);
    
    % 2. Advection (Upwind)
    dT_dy_up = zeros(size(T));
    T_shift = circshift(T, [-1, 0]); 
    mask = Vy < 0;                   
    dT_dy_up(mask) = (T_shift(mask) - T(mask)) / dy; 
    Adv_Term = Const_Adv .* (Vy .* dT_dy_up);
    
    % 3. Source
    Source_Term = Qr_mat ./ (Rho_mat .* Cp_mat);
    
    % 4. Update
    T_new = T + dt * (Diff_Term - Adv_Term + Source_Term);
    
    % 5. Boundaries
    T_new(1, :) = T_surf;               
    T_new(end, :) = T_new(end-1, :) + 0.022 * dy; 
    T_new(:, 1) = T_new(:, 2);          
    T_new(:, end) = T_new(:, end-1);    
    T = T_new;
    
    if mod(n, floor(nt/10)) == 0
        fprintf('Progress: %.0f%%\n', n/nt*100);
    end
end

%% === 6. Visualization and Verification ===
% Create a figure with a white background
figure('Position', [100, 100, 1000, 1000], 'Color', 'w');

% Define drill hole He1 observation data (used in both subplots)
Real_Depth_m = [51, 150, 249, 351, 450, 550, 649, 750, 850];
Real_Temp_C  = [15.1, 15.9, 21.5, 21.0, 25.0, 28.3, 31.6, 34.2, 38.5];
He1_x_pos_km = 2.0; % The drill hole is located at x = 2 km in the granite zone

% ---------------------------------------------------------
% Subplot 1: 2D Thermal Structure
% ---------------------------------------------------------
subplot(2, 1, 1); hold on;

% Plot heatmap
imagesc(x/1000, y/1000, T);
colormap(jet); 
caxis([10, 160]);

% Overlay geological boundaries
contour(x/1000, y/1000, Type_mat, [0.5 0.5], 'k:', 'LineWidth', 1); 
contour(x/1000, y/1000, Type_mat, [1.5 1.5], 'k-', 'LineWidth', 2); 
contour(x/1000, y/1000, Type_mat, [2.5 2.5], 'm--', 'LineWidth', 2); 

% Overlay key economic isotherms
[C, h] = contour(x/1000, y/1000, T, [50, 70, 100, 120], 'w-', 'LineWidth', 1.5);
clabel(C, h, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

% Add geological labels
text(1, 1.5, 'Granite (High Qr)', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 11);
text(9, 1, 'Sediments (Blanket)', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 11);
text(5.5, 3, 'Fault (Upflow)', 'Color', 'm', 'Rotation', -90, 'FontSize', 10);

% MARKER: Add Drill Hole He1 location to the 2D map
plot(He1_x_pos_km * ones(size(Real_Depth_m)), Real_Depth_m/1000, 'ko', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 6);
text(He1_x_pos_km + 0.3, 0.5, 'He1 Borehole', 'Color', 'y', 'FontWeight', 'bold', 'FontSize', 11);

% Format axes
axis ij; 
pbaspect([2.5 1 1]); 
xlim([0 Lx/1000]); ylim([0 Ly/1000]);
xlabel('Distance (km)'); ylabel('Depth (km)');
%title('A. 2D Thermal Structure & Isotherms', 'FontSize', 14);

% Add colorbar
cb = colorbar; 
ylabel(cb, 'Temperature (\circC)', 'FontSize', 11);

% ---------------------------------------------------------
% Subplot 2: Geothermal Gradient Comparison
% ---------------------------------------------------------
subplot(2, 1, 2); hold on;

% Extract 1D profiles from the 2D temperature matrix
idx_granite = round(2000 / dx);  % x = 2 km
idx_basin   = round(10000 / dx); % x = 10 km
idx_fault   = round(5100 / dx);  % x = 5.1 km

% Plot modeled profiles
plot(T(:, idx_granite), y/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Model: Granite Zone');
plot(T(:, idx_basin),   y/1000, 'b--', 'LineWidth', 2, 'DisplayName', 'Model: Basin Zone');
plot(T(:, idx_fault),   y/1000, 'g-.', 'LineWidth', 2, 'DisplayName', 'Model: Fault Zone');

% Plot observational data for validation
plot(Real_Temp_C, Real_Depth_m/1000, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, ...
    'DisplayName', 'Data: Drill Hole He1');

% Format axes
set(gca, 'YDir', 'reverse');
grid on;
legend('Location', 'SouthWest', 'FontSize', 10);
xlabel('Temperature (\circC)'); ylabel('Depth (km)');
%title('B. Geothermal Gradient Comparison & Validation', 'FontSize', 14);
ylim([0, 4]); 

% Save the final high-resolution figure
% print('Helmsdale_Final_Results', '-dpng', '-r300');